import os
import logging
import numpy as np
from astropy.io import fits

import config
import DrpException
import ReducedDataSet
import grating_eq
import extract_order
import reduce_order
import nirspec_constants as constants
import wavelength_utils

logger = logging.getLogger('obj')

def reduce_frame(raw, out_dir):
    """
    
    raw - RawDataSet object
    out_dir - 
    """
         
    # initialize per-object logger and check output directory
    init(raw.objFileName, out_dir)
    
    # if no flats in raw data set then fail
    if (len(raw.flatFileNames) < 1):
        logger.error("no flats")
        raise StandardError('no flats');
    
    # create reduced data set
    reduced = ReducedDataSet.ReducedDataSet(raw.getObjFileName(), raw.getObjHeader())
    
    # put object summary info into per-object log
    log_start_summary(reduced)
    
    # read raw object data into reduced data set object
    reduced.obj = fits.getdata(raw.objFileName)
    
    # combine flats and darks, if darks exist then subtract from obj and flat,
    # store results in processed data set
    process_darks_and_flats(raw, reduced)
        
    # clean cosmic ray hits
    if config.params['cosmic']:
        logger.info('starting cosmic ray cleaning')
        reduced.cleanCosmicRayHits()
        logger.info('cosmic ray cleaning complete')
    else:
        logger.info("not cleaning cosmic ray hits")
        
    # reduce orders
    reduce_orders(reduced)
    
    try:
        # find wavelength solution
        find_global_wavelength_soln(reduced)
    except DrpException as e:
        logger.info('not applying wavelength solution')
    else:
        # apply wavelength solution
        apply_wavelength_soln(reduced)
    
    return(reduced)
 
 
def reduce_orders(reduced):
    
    starting_order = constants.get_starting_order(reduced.getFilter())
    first_order_found = False
    
    for order_num in range(starting_order, 0, -1):
        
        logger.info('***** order ' + str(order_num) + ' *****')
        
        # get expected location of order on detector
        top_calc, bot_calc, wavelength_scale_calc = \
                grating_eq.evaluate(order_num, reduced.getFilter(), reduced.getSlit(), 
                reduced.getEchPos(), reduced.getDispPos(), reduced.getDate())
        
        logger.info('predicted y location: top = ' + 
                    '{:.0f}'.format(top_calc) + ', bottom = ' + '{:.0f}'.format(bot_calc))
        
        if not grating_eq.is_on_detector(top_calc, bot_calc):
            
            # order is not on the detector
            
            logger.info('order ' + str(order_num) + ' is not on the detector')
            if first_order_found:
                break;

        else:
           
            # order is on the detector
            
            first_order_found = True
            
            try:
                order = extract_order.extract_order(order_num, reduced.obj, reduced.flat, 
                        top_calc, bot_calc, reduced.getFilter(), reduced.getSlit())
            except DrpException as e:
                logger.warning('failed to extract order {}: {}'.format(
                            str(order_num), e.message))
                continue
                
            # put integration time and wavelength scale based on 
            # grating equation into order object
            order.integrationTime = reduced.getIntegrationTime() # used in noise calc
            order.wavelengthScaleCalc = wavelength_scale_calc
            order.wavelengthScaleMeas = wavelength_scale_calc
            
            # reduce order, i.e. rectify, extract spectra, identify sky lines
            reduce_order.reduce_order(order)
                
            # add reduced order to list of reduced orders in Reduced object
            try:
                reduced.orders.append(order)                      
            except Exception as e:
                logger.warning('failed to reduce order {}: {}'.format(
                        str(order_num), e.message))
        
        # end if order is on the detector
    # end for each order
    
    logger.info('{} orders reduced'.format(len(reduced.orders)))
    logger.info('end')

    return
            
def find_global_wavelength_soln(reduced):
    
    # create arrays of col, 1/order, accepted wavelength
    # in future will modify twodfit() to take list of lines rather than these constructed arrays
    col = []
    order_inv = []
    accepted = []
    
    for order in reduced.orders:
        for line in order.lines:
            col.append(line.col)
            order_inv.append(1.0 / order.orderNum)
            accepted.append(line.acceptedWavelength)
            
            
    logger.info('total number of identified lines = ' + str(len(col)))
    
#     print accepted
#     raw_input('waiting')
            
    reduced.coeffs, wave_fit, wave_exp = wavelength_utils.twodfit(
            np.asarray(col, dtype='float32'), 
            np.asarray(order_inv, dtype='float32'), 
            np.asarray(accepted, dtype='float32'))    
    
    if wave_fit is None:
        #raise DrpException.DrpException('cannot find wavelength solution')
        return
    
    logger.info('number of lines used in wavelength fit = ' + str(len(wave_fit)))
    
    # There must be a more pythonic way of doing this
    for i, exp in enumerate(wave_exp):
        found = False
        for order in reduced.orders:
            for line in order.lines:
                if abs(line.acceptedWavelength - exp) <= 0.1:
#                     print(str(line.acceptedWavelength) + " = " + str(exp))
                    line.usedInGlobalFit = True
                    line.globalFitWavelength = wave_fit[i]
                    found = True
                    break;
            if found:
                break
        if not found:
            logger.error('cannot find line for expected wavelength ' + str(exp))
    
#     raw_input('waiting')

    # for each line used in the fit, compute slope of wavelength solution at that point
    for order in reduced.orders:
        for line in order.lines:   
            if line.usedInGlobalFit:
                line.globalFitSlope = \
                        reduced.coeffs[1] + \
                        (2.0 * reduced.coeffs[2] * line.col) + \
                        (reduced.coeffs[4] / order.orderNum) + \
                        (2.0 * reduced.coeffs[5] * line.col / order.orderNum)  
    
#     print('coeffs: ' + str(reduced.coeffs))
#     print('wave_fit: ' + str(wave_fit))
#     print('wave_exp: ' + str(wave_exp))
#     raw_input('waiting')
    
    return

def apply_wavelength_soln(reduced):
    
    if reduced.coeffs is None:
        return
    
    for order in reduced.orders:
        x = np.arange(1024)
        y = 1.0 / order.orderNum
        order.wavelengthScaleMeas = np.ravel(
                reduced.coeffs[0] +
                reduced.coeffs[1] * x +
                reduced.coeffs[2] * x ** 2 +
                reduced.coeffs[3] * y +
                reduced.coeffs[4] * x * y +
                reduced.coeffs[5] * (x **2) * y)
        
    return
    
def init(objFileName, out_dir):
    """
    Sets up per-object logger
    
    Returns True if errors were encountered, othewise returns False
    """    
    
    # set up obj logger
    if config.params['debug']:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s ' +
                '%(levelname)s - %(filename)s:%(lineno)s - %(message)s')
    else:
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
    
    fn = out_dir + '/' + objFileName[objFileName.find("NS"):objFileName.rfind(".")]  + '.log'
        
    if os.path.exists(fn):
        os.rename(fn, fn + '.prev')
        
    fh = logging.FileHandler(filename=fn)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    if config.params['debug']:
        sformatter = logging.Formatter('%(asctime)s %(levelname)s - %(filename)s:%(lineno)s - %(message)s')
    else:
        sformatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
    sh = logging.StreamHandler()
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(sformatter)
    logger.addHandler(sh) 
    
    # if output directory does not exist try to create it
    if not os.path.exists(out_dir):
        try: 
            os.mkdir(out_dir)
        except: 
            logger.critical('output directory cannot be created')
        return    

def log_start_summary(reduced):
    logger.info('starting reduction of ' + reduced.getFileName())
    logger.info('  date of observation: ' + reduced.getDate() + ' UT')
    logger.info('          filter name: ' + reduced.getFilter())
    logger.info('            slit name: ' + reduced.getSlit())
    logger.info('        echelle angle: ' + str(reduced.getEchPos()) + ' deg')
    logger.info('cross disperser angle: ' + str(reduced.getDispPos()) + ' deg')

def process_darks_and_flats(raw, reduced):
    
    # combine flats and darks
    logger.info(str(len(raw.flatFileNames)) + ' flats: ' + 
            ', '.join(str(x) for x in  ([s[s.find("NS"):s.rfind(".")] for s in raw.flatFileNames])))
    reduced.flat = raw.combineFlats()
    logger.info(str(len(raw.flatFileNames)) + ' flats have been median combined')

    if len(raw.darkFileNames) > 0:
        reduced.hasDark = True
        logger.info(str(len(raw.darkFileNames)) + ' darks: ' + 
            ', '.join(str(x) for x in ([s[s.find("NS"):s.rfind(".")] for s in raw.darkFileNames])))
        reduced.dark = raw.combineDarks()
        logger.info(str(len(raw.darkFileNames)) + ' darks have been median combined')
    else:
        logger.info('no darks')

    # if dark(s) exist, subtract from obj and flat
    if reduced.hasDark:
        reduced.subtractDark()
        logger.info('dark subtracted from obj and flat')
        
    return
    