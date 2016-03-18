import numpy as np
from astropy.io import fits

class RawDataSet:
    
    def __init__(self, objFileName, objHeader):
        self.data = []
        self.objFileName = objFileName
        self.objHeader = objHeader
        self.flatFileNames = []
        self.darkFileNames = []
        
    def getObjFileName(self):
        return(self.objFileName)
    
    def getObjHeader(self):
        return(self.objHeader)
    
    def getShape(self):
        return((self.objHeader['NAXIS1'], self.objHeader['NAXIS2']))
    
#     def getFilter(self):
#         return self.objHeader['filname']
#         
#     def getEchPos(self):
#         return self.objHeader['echlpos']
# 
#     def getDispPos(self):
#         return self.objHeader['disppos']
#     
#     def getSlit(self):
#         return self.objHeader['slitname']
#     
#     def getDate(self):
#         return self.objHeader['DATE-OBS']
    
    
    def combineFlats(self):
        """
        """
        if len(self.flatFileNames) == 0:
            return None
        if len(self.flatFileNames) == 1:
            return(fits.getdata(self.flatFileNames[0], ignore_missing_end=True))
        flat_data_list = []
        for filename in self.flatFileNames:
            flat_data_list.append(fits.getdata(filename, ignore_missing_end=True))   
        return np.median(flat_data_list, axis=0)
         
    def combineDarks(self):       
        """
        """
        if len(self.darkFileNames) == 0:
            return None
        if len(self.darkFileNames) == 1:
            return(fits.getdata(self.darkFileNames[0], ignore_missing_end=True))
        if len(self.darkFileNames) > 0:
            dark_data_list = []
            for filename in self.darkFileNames:
                dark_data_list.append(fits.getdata(filename, ignore_missing_end=True))   
            return np.median(dark_data_list, axis=0)
