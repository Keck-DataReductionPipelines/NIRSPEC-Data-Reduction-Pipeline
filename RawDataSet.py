import numpy as np
from astropy.io import fits

class RawDataSet:
    
    def __init__(self, objFileName, objHeader):
        #self.data = []
        self.objFileName = objFileName
        
        if objFileName.find('NS.') == 0:
            self.baseName = objFileName[objFileName.find('NS'):].rstrip('.gz').rstrip('.fits')
        else:
            if objFileName.find('/') >= 0:
                self.baseName = objFileName[objFileName.rfind('/')+1:objFileName.lower().find('.fits')]
            else:
                self.baseName = objFileName[:objFileName.lower().find('.fits')]
                
        self.objHeader = objHeader
        self.flatFileNames = []
        self.darkFileNames = []
        
    def getObjFileName(self):
        return(self.objFileName)
    
    def getObjHeader(self):
        return(self.objHeader)
    
    def getShape(self):
        return((self.objHeader['NAXIS1'], self.objHeader['NAXIS2']))
    
    def toString(self):
        return 'raw data set: fn={}, ut={}, disppos={}, echlpos={}, filname={}, slitname={}'.format(
                self.baseName, self.objHeader['UTC'], self.objHeader['disppos'], 
                self.objHeader['echlpos'], self.objHeader['filname'], self.objHeader['slitname'])
    
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
