# import numpy as np


class Order:

    def __init__(self, frames, baseNames, flatOrder):
        """
        Attributes:

            baseName:

            flatOrder:

            frames:

        """

        self.frames = frames
        self.baseNames = baseNames
        self.flatOrder = flatOrder

        self.integrationTime = 0.0

        """
        Attributes:

            topCalc: Pixel row of top left (x = 0) corner of order predicted by grating equation.

            topMean: Pixel row of top left corder as determined by edge detection.

            botCalc: Pixel row of bottom left corner of order predicted by grating equation.

            botMeas: Pixel row of bottom left corner of order determined by edge detection.

        """
        self.topCalc = 0.0
        self.topMeas = 0.0
        self.botCalc = 0.0

        self.spectralTrace = []

        self.objCutout = {}
        for frame in self.frames:
            self.objCutout[frame] = []

        self.onOrderMask = []
        self.offOrderMask = []

        #
        # fields pertaining to wavelength calibration
        #
        self.lines = []

        self.waveShift = None
        self.waveScale = []
        self.calMethod = 'unknown'

        self.orderCal = False        # true if per-order fit computed
        self.orderCalWaveScale = []  # wavelength scale based on per-order fit
        self.orderCalSlope = 0.0     # slope of linear per-order fit
        self.orderCalIncpt = 0.0     # y-intercept of linear per-order fit
        self.orderCalCorrCoeff = 0.0  # per-order linear fit r^2
        self.orderCalNLines = 0      # number of lines used in per-order fit
        self.orderCalRMSRes = 0.0    # RMS fit residual for per-order fit

        self.frameCalWaveScale = []  # wavelength scale based on per-frame fit

        self.mfCalScale = []

#         self.topEdgeProfiles = None
#         self.botEdgeProfiles = None

        self.flatNormalized = False
        self.flattened = False
        self.spatialRectified = False
        self.spectralRectified = False

        self.flatMean = 0.0

        """

        """
        self.flatImg = []
        self.normalizedFlatImg = []

        self.objImg = {}
        self.ffObjImg = {}

        for frame in self.frames:
            self.objImg[frame] = []
            self.ffObjImg[frame] = []

        self.noiseImg = []

        """
        These attributes are used for storing order images after spatial
        rectification but before spectral rectification for diagnostic purposes.

        Attribute:

            srNormFlatImg: Normalized flat order image before spectral rectification.

            srFlatObjAImg: Flat-fielded object A order image before spectral rectification.

            srFlatObjBImg: Flat-fielded object B order image before spectral rectification.

            srFlatObjABImg: Flat-fielded object A - B order image before spectral rectification.

        """
        self.srNormFlatImg = []

        self.srFfObjImg = {}
        for frame in self.frames:
            self.srFfObjImg[frame] = []

        """
        These attributes are used to store spatial profiles and spatial peak locations.
        In pair subtraction mode profiles and peaks of A and B images are used to
        extract sky spectra for sky line calibration. In single object frame mode,
        only the profile and peak for frame A are used.

        Attributes:

            spatialProfileA, B, AB:

            spatialPeakA, B, AB:

            centroid:

            gaussianParams:

        """

        self.spatialProfile = {}
        self.peakLocation = {}
        self.centroid = {}
        self.gaussianParams = {}

        for frame in frames:
            self.spatialProfile[frame] = []
            self.peakLocation[frame] = 0
            self.centroid[frame] = None
            self.gaussianParams[frame] = None

        """
        These attributes define the extraction windows as lists of pixel row numbers

        Attributes:
            objWindow: Window for extracting object flux spectrum.  In the AB pair
                subtraction case, this is centered on the positive peak in the
                subtracted frame.  In the single frame case, this is centered on the
                spatial peak in the single frame A.

            topSkyWindow: Window for extracting top (higher row numbers) sky window.
                This is always on frame A in pair or single frame cases.

            botSkyWindow: Window for extracting bottom (lower row numbers) sky window.
                Always refers to frame A.
        """
        self.objWindow = {}
        self.topSkyWindow = {}
        self.botSkyWindow = {}

        self.topBgMean = {}
        self.botBgMean = {}
        self.snr = {}

        self.objSpec = {}
        self.skySpec = {}

        for frame in frames:
            self.objWindow[frame] = []
            if frame == 'AB':
                self.topSkyWindow[frame] = None
                self.botSkyWindow[frame] = None
                self.skySpec[frame] = None
            else:
                self.topSkyWindow[frame] = []
                self.botSkyWindow[frame] = []
                self.skySpec[frame] = []
            self.topBgMean[frame] = None
            self.botBgMean[frame] = None
            self.snr[frame] = None
            self.objSpec[frame] = []
        """

        """

        """

        """
        self.flatSpec = []
        self.noiseSpec = {}
        for frame in frames:
            self.noiseSpec[frame] = []
        self.synthesizedSkySpec = []
