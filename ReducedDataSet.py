import numpy as np
# import image_lib

# import Flat

class ReducedDataSet:
    """This class represents a reduced data set consisting of object, flat and dark file names
    and image data, a reduced flat, a list of reduced orders, and various other attributes
    describing the data set and reduction results.
    
    In the psuedo object oriented design of the NSDRP classes are used primarily as data 
    structures and include only trivial methods.  Objects of this class are operated on
    primarily by functions in the reduce_frame module.
    
    Attributes:
        baseName:
        header:
        pair:
        hasDark:
        darkKOAId:
        flatKOAIds:
        darkSubtracted:
        cosmicCleaned:
        objA:
        objB:
        objAB:
        nOrders:
        nOrdersReduced:
        snrMean:
        snrMin:
        wMean:
        wMax:
        orders:
        nLinesFound:
        nLinesUsed:
        frameCalAvailable:
        frameCalRmsRes:
        frameCalCoeffs:
        calFrame:
        Flat:
    
    """
    
    def __init__(self, raw):
                
        self.baseNames = raw.baseNames
        self.header = raw.objHeader
        self.isPair = raw.isPair
        
        if self.isPair:
            self.frames = ['A', 'B', 'AB']
        else:
            self.frames = ['A']
            
            
        self.objImg = {}
        for frame in self.frames:
            self.objImg[frame] = np.zeros(self.getShape())
            
        # save KOA IDs of first dark (if any) and flat(s), these are added
        # to FITS headers later.
        if len(raw.darkFns) > 0:
            self.darkKOAId = raw.darkFns[0]
        else:
            self.darkKOAId = 'none'
            
        self.flatKOAIds = []
        for flat_name in raw.flatFns:
            self.flatKOAIds.append(flat_name[flat_name.rfind('/') + 1:flat_name.rfind('.')])
        
        self.hasDark = False
        self.darkSubtracted = False
        self.cosmicCleaned = False
        

        self.flatImg = np.zeros(self.getShape())
        self.darkImg = np.zeros(self.getShape())
        
        self.nOrders = 0
        self.nOrdersReduced = 0
        
        self.snrMean = None
        self.snrMin = None
        self.wMean = None
        self.wMax = None
        
        self.orders = []
        
        self.nLinesFound = 0
        self.nLinesUsed = 0
        
        self.frameCalAvailable = False
        self.frameCalRmsRes = None  # rms per-frame fit residual
        self.frameCalCoeffs = None  # per-frame wavelength equation coefficients
        self.calFrame = None        # frame (name or KOAID) used for wavelength cal if not this one
        
        self.Flat = None
        
    def getBaseName(self):
        if self.isPair:
            return(self.baseNames['AB'])
        else:
            return(self.baseNames['A'])
        
    def getFileName(self):
        return self.fileName
    
    def getTargetName(self):
        return self.header['TARGNAME']
    
    def getShape(self):
        return self.header['NAXIS1'], self.header['NAXIS2']
        
    def getFilter(self):
        if self.header['filname'].startswith('NIRSPEC') and len(self.header['filname']) > 9:
            return self.header['filname'][:9]
        else:
            return self.header['filname']
        
    def getFullFilterName(self):
        return self.header['filname']
        
    def getEchPos(self):
        return self.header['echlpos']

    def getDispPos(self):
        return self.header['disppos']
    
    def getSlit(self):
        return self.header['slitname']
    
    def getITime(self):
        return self.header['itime']
    
    def getNCoadds(self):
        return self.header['coadds']
    
    def getDate(self):
        return self.header['DATE-OBS']
    
    def getTime(self):
        return self.header['UTC']
    
    def getIntegrationTime(self):
        try:
            return self.header['ELAPTIME']
        except KeyError:
            return(self.header['itime'])
            
    
    def getObjectName(self):
        return self.header['OBJECT']
    
    def subtractDark(self):   
        if self.hasDark:
            frames = ['A']
            if self.isPair:
                frames.append('B')
            for frame in frames:
                self.objImg[frame]  = np.subtract(self.objImg[frame], self.dark)
            self.flatImg = np.subtract(self.flatImg, self.dark)
            self.darkSubtracted = True
            
