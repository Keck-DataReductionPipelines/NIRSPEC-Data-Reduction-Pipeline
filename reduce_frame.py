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
main_logger = logging.getLogger('main')

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
        
    # find order edge peak locations on flat, not fully moved here yet
    find_order_edge_peaks(reduced)
        
    # reduce orders
    try:
        reduce_orders(reduced)
    except IOError as e:
        # might want to do something else here
        raise
    
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
    """
    Successively reduces each order in the frame.  
    
    Starting order is determined from a lookup table indexed by filter name.
    
    The grating equation is evaluated for y-axis location of the short wavelength end
    of the order on the detector.
    
    If the order is on the detector then extract_order() is called to cut out from the full
    frame a rectangular array of pixels containing the entire order plus padding.  
    
    Then, reduce_order() is called to reduce the order.
    reduce_order() returns an order object which is and instance of the Order class 
    and contains all of the reduced data for this order.  The order object is then
    appended to the list of order objects in the ReducedDataSet object representing 
    the current frame.
    
    After the first order that should be on the detector is found, processing continues
    through each order in descending order number, working toward higher pixel row numbers
    on the detector, until the first off-detector order is found.
    """
    starting_order = constants.get_starting_order(reduced.getFilter())
    first_order_found = False
    n_orders_on_detector = 0
    
    for order_num in range(starting_order, 0, -1):
        
        logger.info('***** order ' + str(order_num) + ' *****')
        
        # get expected location of order on detector
        top_calc, bot_calc, wavelength_scale_calc = \
                grating_eq.solve(order_num, reduced.getFilter(), reduced.getSlit(), 
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
            
            n_orders_on_detector += 1
            
            order = extract_order.extract_order(order_num, reduced.obj, reduced.flat, 
                    top_calc, bot_calc, reduced.getFilter(), reduced.getSlit())
            
            if order is None:
                logger.warning('failed to extract order {}'.format(str(order_num)))
                continue

            # put integration time and wavelength scale based on 
            # grating equation into order object
            order.integrationTime = reduced.getIntegrationTime() # used in noise calc
            order.wavelengthScaleCalc = wavelength_scale_calc
            order.wavelengthScaleMeas = wavelength_scale_calc
            
            try:
                # reduce order, i.e. rectify, extract spectra, identify sky lines
                reduce_order.reduce_order(order)
                
                # add reduced order to list of reduced orders in Reduced object
                reduced.orders.append(order)                      
#             except DrpException as e:
            except Exception as e:
                logger.warning('failed to reduce order {}: {}'.format(
                        str(order_num), e.message))
 
                        
        # end if order is on the detector
    # end for each order
    
    msg = '{} orders on the detector'.format(n_orders_on_detector)
    logger.info(msg)
    main_logger.info(msg)
    msg = '{} orders reduced'.format(len(reduced.orders))
    logger.info(msg)
    main_logger.info(msg)

    logger.info('end of frame reduction')

    return
            
def find_order_edge_peaks(reduced):
    
    from scipy.signal import argrelextrema

    # make top and bottom edge images
    rolled = np.roll(reduced.flat, 5, axis=0)
    reduced.topEdgesImg = rolled - reduced.flat
    reduced.botEdgesImg = reduced.flat - rolled
    
    
    # take a vertical cut of edges
    reduced.topEdgesProfile = np.median(reduced.topEdgesImg[:, 40:50], axis=1)
    reduced.botEdgesProfile = np.median(reduced.botEdgesImg[:, 40:50], axis=1)

    # find the highest peaks in crosscut, search +/- 15 pixels to narrow down list
    top_extrema = argrelextrema(reduced.topEdgesProfile, np.greater, order=35)[0]
    bot_extrema = argrelextrema(reduced.botEdgesProfile, np.greater, order=35)[0]

    # find crosscut values at those extrema
    top_intensities = reduced.topEdgesProfile[top_extrema]
    bot_intensities = reduced.botEdgesProfile[bot_extrema]

    reduced.topEdgePeaks = zip(top_extrema, top_intensities)
    reduced.botEdgePeaks = zip(bot_extrema, bot_intensities)
    
    return

    
def find_global_wavelength_soln(reduced):
    
    # create arrays of col, 1/order, accepted wavelength
    # in future will modify twodfit() to take list of lines rather than these constructed arrays
    col = []
    centroid = []
    order_inv = []
    accepted = []
    
    for order in reduced.orders:
        for line in order.lines:
            col.append(line.col)
            centroid.append(line.centroid)
            order_inv.append(1.0 / order.orderNum)
            accepted.append(line.acceptedWavelength)
            
            
    logger.info('total number of identified lines = ' + str(len(col)))
    
    if config.params['int_c'] is True:
        logger.warning('using integer column numbers in wavelength fit')
        c = np.asarray(col, dtype='float32')
    else:
        c = np.asarray(centroid, dtype='float32')
            
    reduced.coeffs, wave_fit, wave_exp, reduced.rmsFitRes = wavelength_utils.twodfit(
            c, 
            np.asarray(order_inv, dtype='float32'), 
            np.asarray(accepted, dtype='float32'))    
    
    if wave_fit is None:
        #raise DrpException.DrpException('cannot find wavelength solution')
        return
    
    msg = '{} lines used in wavelength fit'.format(len(wave_fit))
    logger.info(msg)
    main_logger.info(msg)
    msg = 'rms wavelength fit residual = {:.3f}'.format(reduced.rmsFitRes)
    logger.info(msg)
    main_logger.info(msg)
    
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
    logger.handlers = []
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
    main_logger = logging.getLogger('main')
    s = 'starting reduction of ' + reduced.getFileName()[reduced.getFileName().rfind('/') + 1:]
    logger.info(s)
    main_logger.info(s)
    logger.info('  date of observation: ' + reduced.getDate() + ' UT')
    s = '          filter name: ' + reduced.getFilter()
    logger.info(s)
    main_logger.info(s)
    s = '            slit name: ' + reduced.getSlit()
    logger.info(s)
    main_logger.info(s)
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
    