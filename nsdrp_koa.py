import logging
import os
import sys
from astropy.io import fits

import RawDataSet
import FlatCacher
import create_raw_data_sets
import config
import reduce_frame
import DrpException
import products
import dgn
import NirspecConfig


def process_dir(in_dir, base_out_dir):
    """
    NSDRP. Assembles raw data sets from FITS files in the input directory,
    then generates reduced data sets from raw data sets.  Level 1 data products
    are generated from the reduced data sets and saved in the output directory.
    """
        
    # use main logger for outer logging, reduce_frame makes the per frame logger
    logger = logging.getLogger('main')
    
    # open summary spreadsheet file
    ssFptr = open(base_out_dir + '/nsdrp_summary.csv', 'w')
    ssFptr.write(', '.join(['fn', 'date', 'target', 'filter', 'slit', 'disppos', 'echlpos',
            'itime', 'coadds', 'n orders', 'n reduced', 'snr mean', 'snr min', 'width mean', 
            'width max', 'n lines found', 'n lines used', 'r rms', 'cal frame']) + '\n')
    ssFptr.flush()
    
    # get list of raw data sets arranged in chronological order
    rawDataSets = create_raw_data_sets.create(in_dir)
    
    # keep track of how many data sets have been successfully reduced
    n_reduced = len(rawDataSets)
    logger.info(str(len(rawDataSets)) + " raw data set(s) assembled")
        
    # instantiate a flat cacher so each flat or set of flats is only reduced once 
    flatCacher = FlatCacher.FlatCacher(logger, base_out_dir + '/flats')
    
    reducedDataSets = []
    nirspecConfig = None
    
    # process each raw data set   
    for rawDataSet in rawDataSets:
        
        # determine the output directory to use for this frame
        out_dir = get_out_dir(rawDataSet.objFileName, base_out_dir);
        
        if nirspecConfig is None:
            nirspecConfig = NirspecConfig.NirspecConfig(rawDataSet.objHeader)
        else:
            if nirspecConfig.isTheSame(rawDataSet.objHeader) == False:
    
                if len(reducedDataSets) > 0:
    
                    if len(reducedDataSets) > 1:
                        mcal(reducedDataSets)
    
                    gen_data_products(reducedDataSets, nirspecConfig, out_dir, ssFptr)
    
                    del reducedDataSets[:] 
                    logger.info('starting new multi-frame set for nirpsec config: {}', 
                            nirspecConfig.toString())
                    nirspecConfig = NirspecConfig.NirspecConfig(rawDataSet.objHeader)
            
        try:
            reducedDataSets.append(reduce_frame.reduce(rawDataSet, out_dir, flatCacher))  
        except DrpException as e:
            n_reduced -= 1
            logger.error('failed to reduce {}: {}'.format(
                    rawDataSet.objFileName, e.message))
        except IOError as e:
            logger.critical('DRP failed due to I/O error: {}'.format(str(e)))
            sys.exit(1)
            
    if len(reducedDataSets) > 1:
        logger.info('doing multi-frame calibration on {} frames, nirspec config: {}'.format(
                len(reducedDataSets), nirspecConfig.toString()))
        mcal(reducedDataSets)

    gen_data_products(reducedDataSets, nirspecConfig, base_out_dir, ssFptr)
            
    if len(rawDataSets) > 0:
        logger.info('n object frames reduced = {}/{}'.format(
                n_reduced, len(rawDataSets)))   
        
    ssFptr.close()
        
    logger.info('end nsdrp')
    return    

def gen_data_products(reducedDataSets, nirspecConfig, base_out_dir, ssFptr):
    
    logger = logging.getLogger('main')

    if len(reducedDataSets) > 1:
        logger.info('generating data products for multi-frame set, nirspec config: {}'.format(
                nirspecConfig.toString()))
    else:
        logger.info('generating data products for single-frame set, nirspec config: {}'.format(
                 nirspecConfig.toString()))
    
    if config.params['no_products'] is True:
        logger.info('data product generation inhibited by command line switch')
    else:
        for reducedDataSet in reducedDataSets:
            products.gen(reducedDataSet, get_out_dir(reducedDataSet.fileName, base_out_dir))

    if config.params['dgn'] is True:
        logger.info('diagnostic mode enabled, generating diagnostic data products')
        for reducedDataSet in reducedDataSets:
            dgn.gen(reducedDataSet, get_out_dir(reducedDataSet.fileName, base_out_dir))    
            
    for reducedDataSet in reducedDataSets:
        append_to_summary_ss(reducedDataSet, ssFptr)    
                   
    return
    
def mcal(reducedDataSets):
     
    logger = logging.getLogger('main')

    logger.info('running mcal')
    
    allCalibrated = True
    for rds in reducedDataSets:
        if rds.frameCalAvailable is False:
            allCalibrated = False
            break
        
    if allCalibrated is True:
        logger.info('all frame have been calibrated, no adjacent frame calibration required')
        return
    
    coeffs = None
    minRes = 1000.0
    calFrame = None
    for rds in reducedDataSets:
        if rds.frameCalAvailable is True:
            if rds.frameCalRmsRes < minRes:
                minRes = rds.frameCalRmsRes
                coeffs = rds.frameCalCoeffs
                calFrame = rds.baseName
         
    if coeffs is None:
        logger.info('calibration frame not found to calibrate adjacent frames')
        return
    
    for rds in reducedDataSets:
        if rds.frameCalAvailable is False:
            rds.frameCalCoeffs = coeffs
            rds.frameCalRmsRes = minRes
            rds.calFrame = calFrame
            reduce_frame.apply_wavelength_soln(rds)
            logger.info('{} wavelength calibration applied to {}'.format(
                    calFrame, rds.baseName))
    
    return
        
        
def get_out_dir(fn, base_out_dir):
    
    if config.params['subdirs'] is True:
        fn = fn.rstrip('.gz').rstrip('.fits')
        if config.params['shortsubdir']:
            out_dir = base_out_dir + '/' + fn[fn[:fn.rfind('.')].rfind('.')+1:]
        else:
            out_dir = base_out_dir + '/' + fn[fn.rfind('/'):]
    else:
        out_dir = base_out_dir   
             
    if not os.path.exists(out_dir):
        try: 
            os.mkdir(out_dir)
        except: 
            msg = 'output directory {} does not exist and cannot be created'.format(out_dir)
            # logger.critical(msg) can't create log if no output directory
            raise IOError(msg)
        
    return(out_dir)

def append_to_summary_ss(reduced, ssFptr):
    v = []
    v.append(reduced.baseName)
    v.append(reduced.getDate())
    v.append(reduced.getTargetName())
    v.append(reduced.getFilter())
    v.append(reduced.getSlit())
    v.append('{:.2f}'.format(reduced.getDispPos()))
    v.append('{:.2f}'.format(reduced.getEchPos()))
    v.append('{:f}'.format(reduced.getITime()))
    v.append('{:d}'.format(reduced.getNCoadds()))
    v.append('{:d}'.format(reduced.Flat.nOrdersExpected))
    v.append('{:d}'.format(reduced.Flat.nOrdersFound))
    v.append('{:.1f}'.format(reduced.snrMean))
    v.append('{:.1f}'.format(reduced.snrMin))
    if reduced.wMean is None:
        v.append(' ')
    else:
        v.append('{:.1f}'.format(reduced.wMean))
    if reduced.wMax is None:
        v.append(' ')
    else:
        v.append('{:.1f}'.format(reduced.wMax))
    v.append('{:d}'.format(reduced.nLinesFound))
    v.append('{:d}'.format(reduced.nLinesUsed))
    if reduced.frameCalRmsRes is None:
        v.append(' ')
    else:
        v.append('{:.3f}'.format(reduced.frameCalRmsRes))
    v.append('{}'.format(reduced.calFrame))
    ssFptr.write(', '.join(v) + '\n')
    ssFptr.flush()
    return
