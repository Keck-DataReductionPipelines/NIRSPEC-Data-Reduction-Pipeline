import numpy as np
from scipy.signal import argrelextrema
from astropy.io import fits
import logging
import traceback

import config
import grating_eq
from DrpException import DrpException
import FlatOrder
import nirspec_lib
import image_lib

logger = logging.getLogger('obj')

class Flat:
    
#     def __init__(self, filename):
#         
#         self.filename = filename
#         
#         self.baseName = self.getBaseName()
#         logger.info('constructing Flat object for {}'.format(self.baseName))
# 
#         self.header = fits.getheader(filename)
#         self.flatImg = fits.getdata(filename)

    def __init__(self, filename, data):
        
        self.filename = filename
        
        self.baseName = self.getBaseName()
        logger.info('constructing Flat object for {}'.format(self.baseName))

        self.header = fits.getheader(filename)
#         self.flatImg = fits.getdata(filename)
        self.flatImg = data

        self.filterName = self.header['filname']
        self.slit = self.header['slitname']
        self.echelleAngle = self.header['echlpos']
        self.disperserAngle = self.header['disppos']
        
        self.topEdgeImg = None          # top edge image (shift subtract)
        self.botEdgeImg = None          # bottom edge image (shift subtract)
        self.topEdgeProfile = None     # top edge profiles
        self.botEdgeProfile = None     # bottom edge profiles
        self.topEdgePeaks = None        # filtered top edge profile peaks
        self.botEdgePeaks = None        # filtered bottom edge profile peaks
        
        
        self.flatOrders = []
        
        try:
            self.reduce()
        except Exception as e:
            logger.error('flat reduction failed: ' + e.message)
#             traceback.print_tb()
            raise
        
    def getShape(self):
        return((self.objHeader['NAXIS1'], self.objHeader['NAXIS2']))   
    
    def reduce(self):
        """
        """
        logger.info('reducing flat {}'.format(self.baseName))
        
        self.findEdgeProfilePeaks()
        
        nOrdersExpected = 0      
        firstOrderFound = False
        
        for orderNum in range(config.get_starting_order(self.filterName), 0, -1):
            
            logger.info('***** flat order {} *****'.format(orderNum))

            flatOrder = FlatOrder.FlatOrder(orderNum)
            
            # get expected location of order on detector
            flatOrder.topCalc, flatOrder.botCalc, flatOrder.waveScaleCalc = grating_eq.evaluate(
                    orderNum, self.filterName, self.slit, self.echelleAngle, self.disperserAngle)
            
            logger.info('predicted y location: top = ' + '{:.0f}'.format(flatOrder.topCalc) + 
                    ', bottom = ' + '{:.0f}'.format(flatOrder.botCalc))
            
            # determine if order is expected to be on the detector
            # if this order is off but previous order(s) was/were on then no more orders
            # because orders are contiguous
            if not grating_eq.is_on_detector(flatOrder.topCalc, flatOrder.botCalc):
                logger.info('order {} is not on the detector'.format(orderNum))
                if firstOrderFound:
                    break
                
            else:
                firstOrderFound = True
                nOrdersExpected += 1
                
                # determine top and bottom LHS of order by edge detection
                self.findOrder(flatOrder)
                
                if flatOrder.topMeas is None and flatOrder.botMeas is None:
                    continue
                
                # find spatial trace from edge traces
                try:
                    self.findSpatialTrace(flatOrder)
                except DrpException as e:
                    logger.info('failed to find spatial trace: {}'.format(e.message))
                    flatOrder.valid = False
                    continue
   
                if flatOrder.spatialTraceFitResidual > config.params['max_spatial_trace_res']:
                    logger.info('spatial trace fit residual too large, limit = {}'.format(
                            config.params['max_spatial_trace_res']))
                    flatOrder.valid = False
                    continue    
                                 
                try:
                    self.cutOutOrder(flatOrder)
                except DrpException as e:
                    logger.warning('failed to extract flat order {}: {}'.format(
                            str(orderNum), e.message))
                    flatOrder.valid = False
                    continue
                
                try:
                    flatOrder.reduce()
                except DrpException as e:
                    logger.warning('failed to reduce flat order{}: {}'.format(
                        str(orderNum), e.message))
                    flatOrder.valid = False
                    continue
                
                flatOrder.valid = True
                logger.debug('flat order {} validated'.format(orderNum))
                self.flatOrders.append(flatOrder)
                        
        logger.info('flat reduction compete')
        logger.info('n orders expected = {}'.format(nOrdersExpected))
        logger.info('n orders found = {}'.format(len([p for p in self.flatOrders if p.valid == True])))
        return
        
    def getBaseName(self):
        return self.filename[self.filename.rfind('/') + 1:self.filename.rfind('.')]
        
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
        tall_peaks_i = np.where(peak_intensities > (np.amax(peak_intensities) * 0.10))
        
        return(peak_rows[tall_peaks_i[0]])
        
        
    def findOrder(self, flatOrder):
        
        flatOrder.topMeas = min((abs(flatOrder.topCalc - i), i) for i in self.topEdgePeaks)[1]
        flatOrder.botMeas = min((abs(flatOrder.botCalc - i), i) for i in self.botEdgePeaks)[1]
        
        max_delta = config.get_max_edge_location_error(self.filterName, self.slit)

        if flatOrder.topMeas is None or abs(flatOrder.topMeas - flatOrder.topCalc) > max_delta:
            logger.info('measured top edge location too far from expected location')
            logger.info('\tcalc={:.0f}, meas={:.0f}, delta={:.0f}, max delta={:.0f}'.format(
                    flatOrder.topCalc, flatOrder.topMeas, 
                    abs(flatOrder.topMeas - flatOrder.topCalc), max_delta))
            flatOrder.topMeas = None
            topStr = 'not found'
        else:
            topStr = str(flatOrder.topMeas)
            
        if flatOrder.botMeas is None or abs(flatOrder.botMeas - flatOrder.botCalc) > max_delta:
            logger.info('measured bottom edge location too far from expected location')
            logger.info('\tcalc={:.0f}, meas={:.0f}, delta={:.0f}, max delta={:.0f}'.format(
                    flatOrder.botCalc, flatOrder.botMeas, 
                    abs(flatOrder.botMeas - flatOrder.botCalc), max_delta))
            flatOrder.botMeas = None
            botStr = 'not found'
        else:
            botStr = str(flatOrder.botMeas)

        logger.info('measured y location:  top = ' + topStr + ', bottom = ' + botStr)
        
        return
        
        
    def findSpatialTrace(self, flatOrder):   
         
        if flatOrder.topMeas is not None:
            logger.debug('tracing top of order')
            flatOrder.topEdgeTrace = nirspec_lib.trace_order_edge(self.topEdgeImg, flatOrder.topMeas)
        
        if flatOrder.botMeas is not None:
            logger.debug('tracing bottom of order')
            flatOrder.botEdgeTrace = nirspec_lib.trace_order_edge(self.botEdgeImg, flatOrder.botMeas)
            
        if flatOrder.topEdgeTrace is None and flatOrder.botEdgeTrace is None:
            raise DrpException('could not trace top or bottom edge')
    
        if flatOrder.topEdgeTrace is not None and flatOrder.botEdgeTrace is not None:
            logger.info('using top and bottom trace')
            flatOrder.avgEdgeTrace = (flatOrder.topEdgeTrace + flatOrder.botEdgeTrace) / 2.0
    
        elif flatOrder.botEdgeTrace is None:
            logger.info('using top trace only')
            flatOrder.avgEdgeTrace = flatOrder.topEdgeTrace - \
                    ((flatOrder.topMeas - flatOrder.botCalc) / 2.0) + 1.0
            
        else:
            logger.info('using bottom trace only')
            flatOrder.avgEdgeTrace = flatOrder.botEdgeTrace + \
                    ((flatOrder.topCalc - flatOrder.botMeas) / 2.0) + 1.0
                    
        # apply long slit edge margin correction to raw traces
        if '24' in self.slit:
            logger.info('applying long slit edge margins of {} pixels'.format(
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
                
        logger.info('spatial trace smoothed, ' + \
                str(self.flatImg.shape[1] - np.count_nonzero(flatOrder.spatialTraceMask)) + 
                ' points ignored')
    
        flatOrder.spatialTraceFitResidual = np.sqrt(
                np.mean(np.square(flatOrder.avgEdgeTrace - flatOrder.smoothedSpatialTrace)))
        logger.info('spatial trace smoothing rms fit residual = {:.2f}'.format(
                flatOrder.spatialTraceFitResidual))
            
        return
        
        
    def cutOutOrder(self, flatOrder):
        
        # determine cutout padding
        flatOrder.cutoutPadding = config.get_cutout_padding(self.filterName, self.slit)
        
        # add extra padding for orders with large tilt
        tilt = abs(flatOrder.avgEdgeTrace[0] - flatOrder.avgEdgeTrace[-1])
        if  tilt > config.params['large_tilt_threshold']:
            logger.info('large order tilt detected, tilt = ' + str(round(tilt, 1)) + 
                ' threshold = ' + str(config.params['large_tilt_threshold']) + 
                ' extra padding = ' + str(config.params['large_tilt_extra_padding']))
            flatOrder.cutoutPadding += config.params['large_tilt_extra_padding']
        logger.debug('cutout padding = ' + str(round(flatOrder.cutoutPadding, 0)))
        
        # determine highest point of top trace (ignore edge)
        if flatOrder.topEdgeTrace is None:
            flatOrder.topEdgeTrace = flatOrder.botEdgeTrace + \
                (flatOrder.topCalc - flatOrder.botCalc) - 5
            
        flatOrder.highestPoint = np.amax(flatOrder.topEdgeTrace[0:-config.params['overscan_width']])
            
        if flatOrder.botEdgeTrace is None:
            flatOrder.botEdgeTrace = flatOrder.topEdgeTrace - \
                    (flatOrder.topCalc - flatOrder.botCalc) + 5
            
        flatOrder.lowestPoint = np.amin(flatOrder.botEdgeTrace[0:-config.params['overscan_width']])
             
        flatOrder.cutout = np.array(image_lib.cut_out(
                self.flatImg, flatOrder.highestPoint, flatOrder.lowestPoint, 
                flatOrder.cutoutPadding))
        
        flatOrder.shiftOffset = flatOrder.cutoutPadding + flatOrder.botMeas
        
        if float(flatOrder.lowestPoint) > float(flatOrder.cutoutPadding):
            flatOrder.onOrderMask, flatOrder.offOrderMask = get_masks(
                    flatOrder.cutout.shape, 
                    flatOrder.topEdgeTrace - flatOrder.lowestPoint + flatOrder.cutoutPadding, 
                    flatOrder.botEdgeTrace - flatOrder.lowestPoint + flatOrder.cutoutPadding)
        else:
            flatOrder.onOrderMask, flatOrder.offOrderMask = get_masks(
                    flatOrder.cutout.shape, flatOrder.topEdgeTrace, flatOrder.botEdgeTrace)
            
        flatOrder.cutout = np.ma.masked_array(flatOrder.cutout, mask=flatOrder.offOrderMask)
    
        return
    
    
def get_masks(shape, top_trace, bot_trace):
    
    y, x = np.indices(shape, dtype=np.float32)
    
    off_top = y > top_trace
    off_bot = y < bot_trace
    off_order = off_top | off_bot

    belowtop = y < top_trace
    abovebot = y > bot_trace
    on_order = belowtop & abovebot
        
    return on_order, off_order   
        
        
        