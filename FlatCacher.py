import os
from astropy.io import fits
import numpy as np
import image_lib
import Flat
import config
from DrpException import DrpException

class FlatCacher:
    """
    """
    
    def __init__(self, logger, flatDir):
        """
        """    
        logger.info('instantiating flat cacher')
        self.logger = logger
        self.flatDir = flatDir  # directory in which FITS and pickle jars are stored
        self.fnToRawFlatList = dict()
        self.fnToFlat = dict()
        
        if not os.path.exists(flatDir):
            try: 
                os.makedirs(flatDir)
            except: 
                msg = 'flat output dir {} does not exist and cannot be created'.format(flatDir)
                # logger.critical(msg) can't create log if no output directory
                raise IOError(msg)
        
    def getFlat(self, fns):
        """    
        """
        
        # If reduced flat exists in memory cache then return it
        
        # If reduced flat exists on disk, then restore it, cache it and return it
        
        # If reduced flat does not exist, create it, store in memory cache,
        # store on disk and return it
        
#         self.logger.debug('looking for flat: ' + 
#             ', '.join(str(x) for x in ([fn[fn.find("NS"):fn.rfind(".")] for fn in fns])))

        baseFns = self.fnsToBaseFns(fns)
        
        fn = self.getCachedFlat(baseFns)
        
        if fn is None:
            fn = self.constructFlatName(baseFns)
            hdu = self.combineFlats(fns)
            self.logger.debug('creating flat ' + fn + ' which combines ' + 
                    ', '.join(str(x) for x in baseFns))            
            flat = Flat.Flat(fn, baseFns, hdu.header, hdu.data, self.flatDir)
            self.fnToFlat[fn] = flat
            self.fnToRawFlatList[fn] = baseFns
            for k in ['GAIN.SPE', 'FREQ.SPE']:
                if hdu.header[k]:
                    hdu.header.pop(k)
            hdu.writeto(self.flatDir + '/' + fn, clobber=True)
            return(flat)
        else:
            self.logger.debug('reusing flat ' + fn)
            return(self.fnToFlat[fn])
        
    def getCachedFlat(self, baseFns):
        
        for k, v in self.fnToRawFlatList.items():
            if set(baseFns) == set(v):
                return k 
        
        return None
        
        
    def reducedFlatExists(self, fns):
        """
        Determine if a reduced flat constructed of the specified list of flats already exists.
        The flat may exist in memory cache or on disk.
        """
        return False
        
        
    def saveReducedFlat(self):
        """
        """ 
        
    def restoreReducedFlat(self):
        """
        """
        
    def constructFlatName(self, baseFns):
        """
        """
        for n in range(8):
            fn = sorted(baseFns)[0][:sorted(baseFns)[0].find('.fits')] + '_flat_' + str(n) + '.fits'
            if fn not in self.fnToRawFlatList:
                return fn
        raise DrpException('too many similar flats')
        
        
    def combineFlats(self, fns):
        """
        Takes a list of flat file names, combines them and returns the 
        resulting image in the form of a new FITS HDU.  
        Keywords are added the new FITS: CLEANED, FLAT0, FLAT1, ...
        """
        if len(fns) == 0:
            return None
    
        hdr = fits.getheader(fns[0])
        for n in range(len(fns)):
            hdr['FLAT' + str(n)] = fns[n][fns[n].rfind('/')+1:fns[n].rfind('.')]
        
        
        # If there is more than one flat in the list then median combine them.
        if len(fns) == 1:
            data = fits.getdata(fns[0], ignore_missing_end=True)
        else:
            data = np.median([fits.getdata(fn) for fn in fns], axis=0)

            
        # If there are fewer than 3 flats in the list, then filter image
        hdr['CLEANED'] = 'No'
        if len(fns) < 3:
            if config.params['no_cosmic']:
                self.logger.info('flat cosmic ray cleaning inhibited by command line flag')
            else:
                self.logger.info('cosmic cleaning flat because < 3 flats were median combined')
                data = image_lib.cosmic_clean(data)
                hdr['CLEANED'] = 'yes'
                self.logger.info('cosmic ray cleaning complete')
         
        return(fits.PrimaryHDU(data, hdr))
        
        
    def fnsToBaseFns(self, fns):
        return([f[f.rfind('/')+1:f.rfind('.')] for f in fns])

    
        
        
    