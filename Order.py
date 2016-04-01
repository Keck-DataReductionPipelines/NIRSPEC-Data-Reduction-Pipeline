#import numpy as np
from __builtin__ import False

class Order:
    
    def __init__(self, orderNum):
                
        self.orderNum = orderNum
        self.integrationTime = 0.0
        self.yOffset = -1
        
        # pixel location of left hand end of top and bottom of order 
        # calc are locations based on evaluation of grating equation
        # meas are locations from analyzing the image and finding the 
        # actual edge of the order on the detector
        self.topCalc = 0.0
        self.topMeas = 0.0
        self.botCalc = 0.0
        self.botMeas = 0.0
        
        self.wavelengthScaleCalc = []
        self.wavelengthShift = None
        self.wavelengthScaleMeas = None
        
        # per-order wavelength calibration  
        self.perOrderCal = False
        self.perOrderSlope = 0.0;
        self.perOrderIntercept = 0.0;
        self.perOrderCorrCoeff = 0.0;
        
        self.topEdgeProfiles = None
        self.botEdgeProfiles = None
        
        # wavelength calibration lines, array of class Line
        self.lines = []
        
        self.padding = 0
        self.objCutout = []
        self.flatCutout = []
        self.onOrderMask = []
        self.offOrderMask = []
        self.highestPoint = None
        self.lowestPoint = None
        
        self.topTrace = None
        self.botTrace = None
        self.avgTrace = []
        self.smoothedTrace = []
        self.traceMask = []
        self.spectralTrace = []
        
        self.flatNormalized = False
        self.flattened = False
        self.spatialRectified = False
        self.spectralRectified = False

        self.flatMean = 0.0
        
        self.flatImg = []
        self.normalizedFlatImg = []
        self.objImg = []
        self.flattenedObjImg = []
        self.noiseImg = []
        
        #
        self.srNormFlatImg = []
        self.srFlatObjImg = []
        #
        
        self.spatialProfile = []
        self.peakLocation = 0
        self.centroid = None
        self.gaussianParams = None
        
        self.objWindow = []
        self.topSkyWindow = []
        self.botSkyWindow = []
        
        self.topBgMean = None
        self.botBgMean = None
        self.snr = None
        
        self.objSpec = []
        self.flatSpec = []
        self.noiseSpec = []
        self.skySpec = []
        self.synthesizedSkySpec = []