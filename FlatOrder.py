import logging
import numpy as np
import image_lib


# logger = logging.getLogger('obj')

class FlatOrder:
    """

    Top refers to high row numbers, longer wavelengths.
    Bottom refers to low row numbers, shorter wavelengths.
    LHS refers to left hand side of order, low column numbers, shorter wavelengths.
    """

    def __init__(self, baseName, orderNum, logger):

        self.flatBaseName = baseName
        self.orderNum = orderNum
        self.logger = logger

        self.valid = False

        self.topCalc = None              # LHS top row of order, according to grating eq
        self.botCalc = None              # LHS bottom row of order, according to grating eq
        self.gratingEqWaveScale = None  # wavelength scale, according to grating eq

        self.topMeas = None              # measured LHS top row of order
        self.botMeas = None              # measured LHS bottom row of order

        self.topEdgeTrace = None         # top edge trace
        self.botEdgeTrace = None         # bot edge trace
        self.avgEdgeTrace = None

        self.longSlitEdgeMargin = 0
        self.cutoutPadding = 0

        self.highestPoint = None
        self.lowestPoint = None
        self.topTrim = None
        self.botTrim = None

        self.onOrderMask = None
        self.offOrderMask = None

        self.mean = None

        self.cutout = None
#         self.flatImg = None
        self.normFlatImg = None
        self.rectFlatImg = None

        self.normalized = False
        self.spatialRectified = False
        self.spectralRectified = False

        self.smoothedSpatialTrace = None
        self.spatialTraceMask = None
        self.spatialTraceFitResidual = None

    def reduce(self):

        self.logger.info('reducing flat order {}'.format(self.orderNum))

        # normalize flat
        self.normFlatImg, self.mean = image_lib.normalize(
            self.cutout, self.onOrderMask, self.offOrderMask)
        self.normalized = True
        self.logger.info('flat normalized, flat mean = ' +
                         str(round(self.mean, 1)))

        # spatially rectify flat
        self.rectFlatImg = image_lib.rectify_spatial(
            self.normFlatImg, self.smoothedSpatialTrace)

        self.spatialRectified = True

        # compute top and bottom trim points
        self.calcTrimPoints()

        # trim rectified flat order images
        self.rectFlatImg = self.rectFlatImg[self.botTrim:self.topTrim, :]

        self.logger.debug(
            'reduction of flat order {} complete'.format(
                self.orderNum))

        return

    def calcTrimPoints(self):
        if self.lowestPoint > self.cutoutPadding:
            self.topTrim = self.highestPoint - self.lowestPoint + self.cutoutPadding - 3
        else:
            self.topTrim = self.highestPoint - 3
        h = np.amin(self.topEdgeTrace - self.botEdgeTrace)
        self.botTrim = self.topTrim - h + 3
        self.botTrim = int(max(0, self.botTrim))
        self.topTrim = int(min(self.topTrim, 1023))

        return
