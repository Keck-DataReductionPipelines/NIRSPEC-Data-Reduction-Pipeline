import numpy as np
import image_lib

import Flat

RDS_VERSION = 1;

class ReducedDataSet:
    
    def __init__(self, fileName, header):
        
        self.version = RDS_VERSION;
        
        self.fileName = fileName
        if fileName.find('NS.') == 0:
            self.baseName = fileName[fileName.find('NS'):].rstrip('.gz').rstrip('.fits')
        else:
            if fileName.find('/') >= 0:
                self.baseName = fileName[fileName.rfind('/')+1:fileName.lower().find('.fits')]
            else:
                self.baseName = fileName[:fileName.lower().find('.fits')]
        self.header = header
       
        self.hasDark = False
        self.darkKOAId = None
        self.flatKOAIds = []
        self.darkSubtracted = False
        self.cosmicCleaned = False;
        
        self.obj = np.zeros(self.getShape())
        self.flat = np.zeros(self.getShape())
        self.dark = np.zeros(self.getShape())
        
        self.orders = []
        
        self.rmsFitRes = None
        self.coeffs = None
        
        self.Flat = None
        
 
        
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
    
    def getIntegrationTime(self):
        try:
            return self.header['ELAPTIME']
        except KeyError:
            return(self.header['itime'])
            
    
    def getObjectName(self):
        return self.header['OBJECT']
    
    def subtractDark(self):   
        if self.hasDark:
            self.obj  = np.subtract(self.obj, self.dark)
            self.flat = np.subtract(self.flat, self.dark)
            self.darkSubtracted = True
            
