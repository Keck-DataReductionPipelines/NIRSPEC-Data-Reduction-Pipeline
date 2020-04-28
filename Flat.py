import os
import numpy as np
from scipy.signal import argrelextrema
# from astropy.io import fits
import logging
# import traceback

import config
import GratingEq
from DrpException import DrpException
import FlatOrder
import nirspec_lib
import image_lib
from numpy.testing.utils import measure


class Flat:

    def __init__(self, fn, baseFns, header, data, logDir=None):

        self.fn = fn
        self.baseFns = baseFns
        self.header = header
        self.logDir = logDir

        self.baseName = self.getBaseName()

        if logDir is None:
            self.logger = logging.getLogger('obj')
        else:
            self.logger = logging.getLogger('flat')
            self.initLogger()

        self.gratingEq = GratingEq.GratingEq(self.logger)

        names = ', '.join(str(x) for x in self.baseFns)
        self.logger.info('creating {} from {}'.format(self.fn, names))

        self.flatImg = data

        self.filterName = self.header['filname']
        self.slit = self.header['slitname']
        self.echelleAngle = self.header['echlpos']
        self.disperserAngle = self.header['disppos']

        self.topEdgeImg = None          # top edge image (shift subtract)
        self.botEdgeImg = None          # bottom edge image (shift subtract)
        self.topEdgeProfile = None      # top edge profiles
        self.botEdgeProfile = None      # bottom edge profiles
        self.topEdgePeaks = None        # filtered top edge profile peaks
        self.botEdgePeaks = None        # filtered bottom edge profile peaks

        self.nOrdersExpected = 0
        self.nOrdersFound = 0

        self.flatOrders = []

        try:
            self.reduce()
        except Exception as e:
            self.logger.error('flat reduction failed: ' + e.message)
#             traceback.print_tb()
            raise

    def getShape(self):
        return((self.objHeader['NAXIS1'], self.objHeader['NAXIS2']))

    def reduce(self):
        """
        """
        self.logger.info('reducing flat {}'.format(self.baseName))

        self.findEdgeProfilePeaks()

        self.nOrdersExpected = 0
        firstOrderFound = False

        for orderNum in range(
                config.get_starting_order(self.filterName), 0, -1):

            self.logger.info('***** flat order {} *****'.format(orderNum))

            flatOrder = FlatOrder.FlatOrder(
                self.baseName, orderNum, self.logger)

            # get expected location of order on detector
            flatOrder.topCalc, flatOrder.botCalc, flatOrder.gratingEqWaveScale = self.gratingEq.evaluate(
                orderNum, self.filterName, self.slit, self.echelleAngle, self.disperserAngle)

            self.logger.info(
                'predicted top edge location = {:.0f} pixels'.format(
                    flatOrder.topCalc))
            self.logger.info(
                'predicted bot edge location = {:.0f} pixels'.format(
                    flatOrder.botCalc))

            # determine if order is expected to be on the detector
            # if this order is off but previous order(s) was/were on then no more orders
            # because orders are contiguous
            if not self.gratingEq.is_on_detector(
                    flatOrder.topCalc, flatOrder.botCalc):
                self.logger.info(
                    'order {} is not on the detector'.format(orderNum))
                if firstOrderFound:
                    break

            else:
                firstOrderFound = True
                self.nOrdersExpected += 1

                # determine top and bottom LHS of order by edge detection
                if config.params['sowc'] is True:
                    self.findOrderSowc(flatOrder)
                else:
                    self.findOrder(flatOrder)

                if flatOrder.topMeas is None and flatOrder.botMeas is None:
                    continue

                # find spatial trace from edge traces
                try:
                    self.findSpatialTrace(flatOrder)
                except DrpException as e:
                    self.logger.info(
                        'failed to find spatial trace: {}'.format(
                            e.message))
                    flatOrder.valid = False
                    continue

                if flatOrder.spatialTraceFitResidual > config.params['max_spatial_trace_res']:
                    self.logger.info('spatial trace fit residual too large, limit = {}'.format(
                        config.params['max_spatial_trace_res']))
                    flatOrder.valid = False
                    continue

                try:
                    self.cutOutOrder(flatOrder)
                except DrpException as e:
                    self.logger.warning('failed to extract flat order {}: {}'.format(
                        str(orderNum), e.message))
                    flatOrder.valid = False
                    continue

                try:
                    flatOrder.reduce()
                except DrpException as e:
                    self.logger.warning('failed to reduce flat order{}: {}'.format(
                        str(orderNum), e.message))
                    flatOrder.valid = False
                    continue

                flatOrder.valid = True
                self.logger.debug('flat order {} validated'.format(orderNum))
                self.flatOrders.append(flatOrder)

        self.logger.info('flat reduction compete')
        self.logger.info('n orders expected = {}'.format(self.nOrdersExpected))
        self.nOrdersFound = len(
            [p for p in self.flatOrders if p.valid is True])
        self.logger.info('n orders found = {}'.format(self.nOrdersFound))
        return

    def getBaseName(self):
        return self.fn[self.fn.rfind('/') + 1:self.fn.rfind('.')]

    def findEdgeProfilePeaks(self):

        # make top and bottom edge profile images
        rolled = np.roll(self.flatImg, 5, axis=0)
        self.topEdgeImg = rolled - self.flatImg
        self.botEdgeImg = self.flatImg - rolled

        self.topEdgeProfile = np.median(self.topEdgeImg[:, 40:50], axis=1)
        self.botEdgeProfile = np.median(self.botEdgeImg[:, 40:50], axis=1)

        self.topEdgePeaks = self.findPeaks(self.topEdgeProfile)
        self.botEdgePeaks = self.findPeaks(self.botEdgeProfile)

        return

    def findPeaks(self, edgeProfile):

        peak_rows = argrelextrema(edgeProfile, np.greater, order=35)[0]
        peak_intensities = edgeProfile[peak_rows]
        tall_peaks_i = np.where(
            peak_intensities > (
                np.amax(peak_intensities) * 0.10))

        return(peak_rows[tall_peaks_i[0]])

    def findOrderSowc(self, flatOrder):

        flatOrder.topMeas = None
        flatOrder.botMeas = None
        w = flatOrder.topCalc - flatOrder.botCalc

        maxDelta = config.get_max_edge_location_error(
            self.filterName, self.slit)

        flatOrder.topMeas = self.findEdge(flatOrder.topCalc, maxDelta, 'top')

        if flatOrder.topMeas is not None:
            flatOrder.botMeas = self.findEdge(
                flatOrder.topMeas - w, maxDelta, 'bot')
        else:
            flatOrder.botMeas = self.findEdge(
                flatOrder.botCalc, maxDelta, 'bot')

        if flatOrder.topMeas is not None and flatOrder.botMeas is None:
            flatOrder.botMeas = self.findEdge(
                flatOrder.topMeas - w, maxDelta, 'bot')
        elif flatOrder.topMeas is None and flatOrder.botMeas is not None:
            flatOrder.topMeas = self.findEdge(
                flatOrder.botMeas + w, maxDelta, 'top')

        if flatOrder.topMeas is None:
            self.logger.info('top edge not found')
        else:
            self.logger.info(
                'measured top edge location = {:.0f} pixels'.format(
                    flatOrder.topMeas))
            self.logger.info('   top edge location delta = {:.0f} pixels'.format(
                flatOrder.topCalc - flatOrder.topMeas))
        if flatOrder.botMeas is None:
            self.logger.info('bottom edge not found')
        else:
            self.logger.info(
                'measured bot edge location = {:.0f} pixels'.format(
                    flatOrder.botMeas))
            self.logger.info('   bot edge location delta = {:.0f} pixels'.format(
                flatOrder.botCalc - flatOrder.botMeas))
        return

    def findEdge(self, calc, maxDelta, topOrBot):
        if topOrBot == 'bot':
            meas = min((abs(calc - i), i) for i in self.botEdgePeaks)[1]
        else:
            meas = min((abs(calc - i), i) for i in self.topEdgePeaks)[1]

        if meas is None or abs(meas - calc) > maxDelta:
            self.logger.debug('{} edge not found at {:.0f} +/- {:.0f}'.format(
                topOrBot, calc, maxDelta))
            return None
        else:
            return meas

    def findOrder(self, flatOrder):

        flatOrder.topMeas = min((abs(flatOrder.topCalc - i), i)
                                for i in self.topEdgePeaks)[1]
        flatOrder.botMeas = min((abs(flatOrder.botCalc - i), i)
                                for i in self.botEdgePeaks)[1]

        max_delta = config.get_max_edge_location_error(
            self.filterName, self.slit)

        if flatOrder.topMeas is None or abs(
                flatOrder.topMeas - flatOrder.topCalc) > max_delta:
            self.logger.info(
                'measured top edge location too far from expected location')
            self.logger.info('\tcalc={:.0f}, meas={:.0f}, delta={:.0f}, max delta={:.0f}'.format(
                flatOrder.topCalc, flatOrder.topMeas,
                abs(flatOrder.topMeas - flatOrder.topCalc), max_delta))
            flatOrder.topMeas = None
            topStr = 'not found'
        else:
            topStr = str(flatOrder.topMeas)

        if flatOrder.botMeas is None or abs(
                flatOrder.botMeas - flatOrder.botCalc) > max_delta:
            self.logger.info(
                'measured bottom edge location too far from expected location')
            self.logger.info('\tcalc={:.0f}, meas={:.0f}, delta={:.0f}, max delta={:.0f}'.format(
                flatOrder.botCalc, flatOrder.botMeas,
                abs(flatOrder.botMeas - flatOrder.botCalc), max_delta))
            flatOrder.botMeas = None
            botStr = 'not found'
        else:
            botStr = str(flatOrder.botMeas)

        self.logger.info(
            'measured y location:  top = ' +
            topStr +
            ', bottom = ' +
            botStr)

        return

    def findSpatialTrace(self, flatOrder):

        if flatOrder.topMeas is not None:
            self.logger.debug('tracing top of order')
            flatOrder.topEdgeTrace = nirspec_lib.trace_order_edge(
                self.topEdgeImg, flatOrder.topMeas)

        if flatOrder.botMeas is not None:
            self.logger.debug('tracing bottom of order')
            flatOrder.botEdgeTrace = nirspec_lib.trace_order_edge(
                self.botEdgeImg, flatOrder.botMeas)

        if flatOrder.topEdgeTrace is None and flatOrder.botEdgeTrace is None:
            raise DrpException('could not trace top or bottom edge')

        if flatOrder.topEdgeTrace is not None and flatOrder.botEdgeTrace is not None:
            self.logger.info('using top and bottom trace')
            flatOrder.avgEdgeTrace = (
                flatOrder.topEdgeTrace + flatOrder.botEdgeTrace) / 2.0

        elif flatOrder.botEdgeTrace is None:
            self.logger.info('using top trace only')
            flatOrder.avgEdgeTrace = flatOrder.topEdgeTrace - \
                ((flatOrder.topMeas - flatOrder.botCalc) / 2.0) + 1.0

        else:
            self.logger.info('using bottom trace only')
            flatOrder.avgEdgeTrace = flatOrder.botEdgeTrace + \
                ((flatOrder.topCalc - flatOrder.botMeas) / 2.0) + 1.0

        # apply long slit edge margin correction to raw traces
        if '24' in self.slit:
            self.logger.info('applying long slit edge margins of {} pixels'.format(
                config.params['long_slit_edge_margin']))
            if flatOrder.topEdgeTrace is not None:
                flatOrder.topEdgeTrace -= config.params['long_slit_edge_margin']
            if flatOrder.botEdgeTrace is not None:
                flatOrder.botEdgeTrace += config.params['long_slit_edge_margin']

        # if bottom edge trace successful, use to refine LHS bottom location
        if flatOrder.botEdgeTrace is not None:
            flatOrder.botMeas = flatOrder.botEdgeTrace[1]

        # smooth spatial trace
        flatOrder.smoothedSpatialTrace, flatOrder.spatialTraceMask = \
            nirspec_lib.smooth_spatial_trace(flatOrder.avgEdgeTrace)

        self.logger.info('spatial trace smoothed, ' +
                         str(self.flatImg.shape[1] - np.count_nonzero(flatOrder.spatialTraceMask)) +
                         ' points ignored')

        flatOrder.spatialTraceFitResidual = np.sqrt(
            np.mean(np.square(flatOrder.avgEdgeTrace - flatOrder.smoothedSpatialTrace)))
        self.logger.info('spatial trace smoothing rms fit residual = {:.2f}'.format(
            flatOrder.spatialTraceFitResidual))

        return

    def cutOutOrder(self, flatOrder):

        # determine cutout padding
        flatOrder.cutoutPadding = config.get_cutout_padding(
            self.filterName, self.slit)

        # add extra padding for orders with large tilt
        tilt = abs(flatOrder.avgEdgeTrace[0] - flatOrder.avgEdgeTrace[-1])
        if tilt > config.params['large_tilt_threshold']:
            self.logger.info('large order tilt detected, tilt = ' + str(round(tilt, 1)) +
                             ' threshold = ' + str(config.params['large_tilt_threshold']) +
                             ' extra padding = ' + str(config.params['large_tilt_extra_padding']))
            flatOrder.cutoutPadding += config.params['large_tilt_extra_padding']
        self.logger.debug('cutout padding = ' +
                          str(round(flatOrder.cutoutPadding, 0)))

        # determine highest point of top trace (ignore edge)
        if flatOrder.topEdgeTrace is None:
            flatOrder.topEdgeTrace = flatOrder.botEdgeTrace + \
                (flatOrder.topCalc - flatOrder.botCalc) - 5

        flatOrder.highestPoint = np.amax(
            flatOrder.topEdgeTrace[0:-config.params['overscan_width']])

        if flatOrder.botEdgeTrace is None:
            flatOrder.botEdgeTrace = flatOrder.topEdgeTrace - \
                (flatOrder.topCalc - flatOrder.botCalc) + 5

        flatOrder.lowestPoint = np.amin(
            flatOrder.botEdgeTrace[0:-config.params['overscan_width']])

        flatOrder.cutout = np.array(image_lib.cut_out(
            self.flatImg, flatOrder.highestPoint, flatOrder.lowestPoint,
            flatOrder.cutoutPadding))

        if float(flatOrder.lowestPoint) > float(flatOrder.cutoutPadding):
            flatOrder.onOrderMask, flatOrder.offOrderMask = get_masks(
                flatOrder.cutout.shape,
                flatOrder.topEdgeTrace - flatOrder.lowestPoint + flatOrder.cutoutPadding,
                flatOrder.botEdgeTrace - flatOrder.lowestPoint + flatOrder.cutoutPadding)
        else:
            flatOrder.onOrderMask, flatOrder.offOrderMask = get_masks(
                flatOrder.cutout.shape, flatOrder.topEdgeTrace, flatOrder.botEdgeTrace)

        flatOrder.cutout = np.ma.masked_array(
            flatOrder.cutout, mask=flatOrder.offOrderMask)

        return

    def initLogger(self):

        self.logger.handlers = []
        if config.params['debug']:
            self.logger.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s ' +
                                          '%(levelname)s - %(filename)s:%(lineno)s - %(message)s')
        else:
            self.logger.setLevel(logging.INFO)
            formatter = logging.Formatter(
                '%(asctime)s %(levelname)s - %(message)s')

        fn = self.logDir + '/' + self.baseName + '.log'

        if os.path.exists(fn):
            os.rename(fn, fn + '.prev')

        fh = logging.FileHandler(filename=fn)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

        if config.params['verbose'] is True:
            if config.params['debug']:
                sformatter = logging.Formatter(
                    '%(asctime)s %(levelname)s - %(filename)s:%(lineno)s - %(message)s')
            else:
                sformatter = logging.Formatter(
                    '%(asctime)s %(levelname)s - %(message)s')
            sh = logging.StreamHandler()
            sh.setLevel(logging.DEBUG)
            sh.setFormatter(sformatter)
            self.logger.addHandler(sh)


def get_masks(shape, top_trace, bot_trace):

    y, x = np.indices(shape, dtype=np.float32)

    off_top = y > top_trace
    off_bot = y < bot_trace
    off_order = off_top | off_bot

    belowtop = y < top_trace
    abovebot = y > bot_trace
    on_order = belowtop & abovebot

    return on_order, off_order
