import logging
import os
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
        out_dir = get_out_dir(rawDataSet, base_out_dir);
        
        if nirspecConfig is None:
            nirspecConfig = NirspecConfig.NirspecConfig(rawDataSet.objHeader)
        else:
            if nirspecConfig.isTheSame(rawDataSet.objHeader) == False:
    
                if len(reducedDataSets) > 0:
    
                    if len(reducedDataSets) > 1:
                        logger.info('doing multi-frame calibration on {} frames, config = {}',
                                len(reducedDataSets), nirspecConfig.toString())
                        logger.info('MCAL GOES HERE')
    
                    gen_data_products(reducedDataSets, nirspecConfig, out_dir)
    
                    reducedDataSets.clear() 
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
        logger.info('MCAL GOES HERE')

    gen_data_products(reducedDataSets, nirspecConfig, out_dir)
            
    if len(rawDataSets) > 0:
        logger.info('n object frames reduced = {}/{}'.format(
                n_reduced, len(rawDataSets)))   
        
    logger.info('end nsdrp')
    return    

def gen_data_products(reducedDataSets, nirspecConfig, out_dir):
    
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
            products.gen(reducedDataSet, out_dir)

    if config.params['dgn'] is True:
        logger.info('diagnostic mode enabled, generating diagnostic data products')
        for reducedDataSet in reducedDataSets:
            dgn.gen(reducedDataSet, out_dir)                   
    return
    
# def reduce(rawDataSet, out_dir, flatCacher=None):    
#     
#     logger = logging.getLogger('main')
# 
#     # generate reduced data set by reducing raw data set
#     reducedDataSet = reduce_frame.reduce(rawDataSet, out_dir, flatCacher)
#         
#     # produce data products from reduced data set
# #     if config.params['no_products'] is True:
# #         logger.info('data product generation inhibited by command line switch')
# #     else:
# #         products.gen(reducedDataSet, out_dir)
# 
#     # if diagnostics mode is enabled, then produce diagnostic data products
# #     if config.params['dgn'] is True:
# #         logger.info('diagnostic mode enabled, generating diagnostic data products')
# #         dgn.gen(reducedDataSet, out_dir)
        
        
def get_out_dir(rawDataSet, base_out_dir):
    
    if config.params['subdirs'] is True:
        fn = rawDataSet.objFileName.rstrip('.gz').rstrip('.fits')
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