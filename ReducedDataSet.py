import numpy as np
import image_lib

class ReducedDataSet:
    
    def __init__(self, fileName, header):
        
        self.fileName = fileName
        self.baseName = fileName[fileName.find('NS'):fileName.rfind('.')]
        self.header = header
       
        self.hasDark = False
        self.darkKOAId = None
        self.flatKOAIds = []
        self.darkSubtracted = False
        self.cosmicCleaned = False;
        
        self.obj = np.zeros(self.getShape())
        self.flat = np.zeros(self.getShape())
        self.dark = np.zeros(self.getShape())
        
        self.topEdgesImg = None
        self.botEdgesImg = None
        self.topEdgesProfile = None
        self.botEdgesProfile = None
        self.topEdgePeaks = None
        self.botEdgePeaks = None
        
        self.orders = []
        
        self.rmsFitRes = None
        self.coeffs = None
        
    def getFileName(self):
        return self.fileName
    
    def getTargetName(self):
        return self.header['TARGNAME']
    
    def getShape(self):
        return self.header['NAXIS1'], self.header['NAXIS2']
        
    def getFilter(self):
        return self.header['filname']
        
    def getEchPos(self):
        return self.header['echlpos']

    def getDispPos(self):
        return self.header['disppos']
    
    def getSlit(self):
        return self.header['slitname']
    
    def getDate(self):
        return self.header['DATE-OBS']
    
    def getIntegrationTime(self):
        return self.header['ELAPTIME']
    
    def getObjectName(self):
        return self.header['OBJECT']
    
    def subtractDark(self):   
        if self.hasDark:
            self.obj  = np.subtract(self.obj, self.dark)
            self.flat = np.subtract(self.flat, self.dark)
            self.darkSubtracted = True
            
    def cleanCosmicRayHits(self):
        self.obj = image_lib.cosmic_clean(self.obj)
#         self.flat = image_lib.cosmic_clean(self.flat)
        self.cosmicCleaned = True 
