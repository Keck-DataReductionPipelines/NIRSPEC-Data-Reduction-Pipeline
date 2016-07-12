import os
import logging
import numpy as np
from astropy.io import fits

import config
import DrpException
import ReducedDataSet
import reduce_order
import wavelength_utils
import nsdrp
import Flat
from logging import INFO
import Order
import image_lib

logger = logging.getLogger('obj')
# main_logger = logging.getLogger('main')
# main_logger = logging.getLogger('main')

def reduce(raw, out_dir, flatCacher=None):
    """
    raw - RawDataSet object
    out_dir - data product root directory
    """
         
    # initialize per-object logger and check output directory
    init(raw.objFileName, out_dir)
    
    # if no flats in raw data set then fail
    if (len(raw.flatFileNames) < 1):
        logger.error("no flats")
        raise DrpException.DrpException('no flats');
    
    # create reduced data set
    reduced = ReducedDataSet.ReducedDataSet(raw.getObjFileName(), raw.getObjHeader())
    
    # save KOA IDs of first dark (if any) and flat(s), these are added
    # to FITS headers later.
    if len(raw.darkFileNames) > 0:
        reduced.darkKOAId = raw.darkFileNames[0]
    else:
        reduced.darkKOAId = 'none'
        
    for flat_name in raw.flatFileNames:
        reduced.flatKOAIds.append(flat_name[flat_name.rfind('/') + 1:flat_name.rfind('.')])
        
    
    # put object summary info into per-object log
    log_start_summary(reduced)
    
    # read raw object data into reduced data set object
    reduced.obj = fits.getdata(raw.objFileName, ignore_missing_end=True)
    
    # combine flats and darks, if darks exist then subtract from obj and flat,
    # store results in processed data set
    #process_darks_and_flats(raw, reduced)
    
    # reduce flat
    if flatCacher is None:
        # in command line mode with only one flat
        if config.params['no_cosmic']:
            logger.info('cosmic ray rejection on flat inhibited by command line flat')
            flat_data = fits.PrimaryHDU.readfrom(raw.flatFileNames[0], ignore_missing_end=True).data
        else:
            logger.info('starting cosmic ray cleaning on flat')
            flat_data = image_lib.cosmic_clean(fits.PrimaryHDU.readfrom(
                    raw.flatFileNames[0], ignore_missing_end=True).data)
            logger.info('cosmic ray cleaning on flat complete')
        
        reduced.Flat = Flat.Flat(
                raw.flatFileNames[0],
                fits.PrimaryHDU.readfrom(raw.flatFileNames[0], ignore_missing_end=True).header, 
                flat_data)
    else:
        reduced.Flat = flatCacher.getFlat(raw.flatFileNames)
    
    # clean cosmic ray hits on object frame
    if config.params['no_cosmic']:
        logger.info("cosmic ray rejection on object frame inhibited by command line flat")

    else:
        logger.info('starting cosmic ray cleaning on object frame')
        reduced.obj = image_lib.cosmic_clean(reduced.obj)
        reduced.cosmicCleaned = True 
        logger.info('cosmic ray cleaning on object frame complete')
        
    # if darks are available, combine them if there are more than one
    # and subtract from object and flat
    process_darks(raw, reduced)
    
    # reduce orders
    try:
        reduce_orders(reduced)
    except IOError as e:
        # might want to do something else here
        raise
    
    # find wavelength solution
    reload(wavelength_utils)
    if find_global_wavelength_soln(reduced) is True:
        apply_wavelength_soln(reduced)
    else:
        logger.info('not applying wavelength solution')
        for order in reduced.orders:
            order.waveScale = order.gratingEqWaveScale
            order.calMethodUsed = 'grating equation'
    
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
    
    for flatOrder in reduced.Flat.flatOrders:
        if flatOrder.valid is not True:
            continue        
        
        logger.info('***** order ' + str(flatOrder.orderNum) + ' *****')
            
        order = Order.Order(reduced.baseName, flatOrder.orderNum)
            
        order.objCutout = np.array(image_lib.cut_out(reduced.obj, 
                flatOrder.highestPoint, flatOrder.lowestPoint, flatOrder.cutoutPadding))
        
        order.integrationTime = reduced.getIntegrationTime() # used in noise calc
        order.gratingEqWaveScale = flatOrder.gratingEqWaveScale
        #order.wavelengthScaleMeas = flatOrder.waveScaleCalc
        order.topTrace = flatOrder.topEdgeTrace
        order.botTrace = flatOrder.botEdgeTrace
        order.avgTrace = flatOrder.avgEdgeTrace
        order.smoothedTrace = flatOrder.smoothedSpatialTrace
        order.traceMask = flatOrder.spatialTraceMask
        order.flatCutout = flatOrder.cutout
        order.highestPoint = flatOrder.highestPoint
        order.lowestPoint = flatOrder.lowestPoint
        order.padding = flatOrder.cutoutPadding
        order.botMeas = flatOrder.botMeas
            
        try:
            
            # reduce order, i.e. rectify, extract spectra, identify sky lines
            reduce_order.reduce_order(order, flatOrder)
    
            # add reduced order to list of reduced orders in Reduced object
            reduced.orders.append(order)                      

        except DrpException.DrpException as e:
            logger.warning('failed to reduce order {}: {}'.format(
                     str(flatOrder.orderNum), e.message))
 
                        
        # end if order is on the detector
    # end for each order

    loggers = ['obj']
    if config.params['cmnd_line_mode'] is False:
        loggers.append('main')
#     for l in loggers:
#         logging.getLogger(l).log(INFO, 'n orders on the detector = {}'.format(n_orders_on_detector))
#         logging.getLogger(l).log(INFO, 'n orders reduced = {}'.format(len(reduced.orders)))
        
    if len(reduced.orders) == 0:
        return
    
    for l in loggers:
        reduced.snrMean = sum(reduced.orders[i].snr for i in range(len(reduced.orders))) / \
                len(reduced.orders)
        logging.getLogger(l).log(INFO, 'mean signal-to-noise ratio = {:.1f}'.format(
                reduced.snrMean))
    
    snr = []
    for i in range(len(reduced.orders)):
        snr.append(reduced.orders[i].snr)
    
    for l in loggers:
        reduced.snrMin = np.amin(snr)
        logging.getLogger(l).log(INFO, 'minimum signal-to-noise ratio = {:.1f}'.format(
                reduced.snrMin))
    
        try:
            reduced.wMean = sum(abs(reduced.orders[i].gaussianParams[2]) \
                    for i in range(len(reduced.orders))) /  len(reduced.orders)
            logging.getLogger(l).log(INFO, 'mean spatial peak width = {:.1f} pixels'.format(
                    reduced.wMean))
        except:
            logging.getLogger(l).log(logging.WARNING, 'mean spatial peak width = unknown') 
    
        w = []
        for i in range(len(reduced.orders)):
            if reduced.orders[i].gaussianParams is not None:
                w.append(reduced.orders[i].gaussianParams[2])
        try:
            reduced.wMax = np.amax(w)
            logging.getLogger(l).log(INFO, 'maximum spatial peak width = {:.1f} pixels'.format(
                    reduced.wMax))
        except:
            logging.getLogger(l).error('maximum spatial peak width cannot be determined')
            logging.getLogger(l).log(INFO, 'maximum spatial peak width = unknown')


    return
    

    
def find_global_wavelength_soln(reduced):
    
    loggers = ['obj']
    if config.params['cmnd_line_mode'] is False:
        loggers.append('main')
        
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
            accepted.append(line.waveAccepted)
            
    reduced.nLinesFound = len(col)
    for l in loggers:
        logging.getLogger(l).log(INFO, 'n sky lines identified = {:d}'.format(reduced.nLinesFound))
    
    if config.params['int_c'] is True:
        logger.warning('using integer column numbers in wavelength fit')
        c = np.asarray(col, dtype='float32')
    else:
        c = np.asarray(centroid, dtype='float32')
            
    reduced.frameCalCoeffs, wave_fit, wave_exp, reduced.frameCalRmsRes = wavelength_utils.twodfit(
            c, 
            np.asarray(order_inv, dtype='float32'), 
            np.asarray(accepted, dtype='float32'))    
    
    if wave_fit is None:
#         raise DrpException.DrpException('cannot find wavelength solution')
        reduced.frameCalAvailable = False
        return False
    
    reduced.frameCalAvailable = True
    reduced.calFrame = 'self'
    reduced.nLinesUsed = len(wave_fit)
    for l in ['main', 'obj']:
        logging.getLogger(l).log(INFO, 'n lines used in wavelength fit = {:d}'.format(
                reduced.nLinesUsed))
        logging.getLogger(l).log(INFO, 'rms wavelength fit residual = {:.3f}'.format(
                reduced.frameCalRmsRes))
    
    # There must be a more pythonic way of doing this
    for i, exp in enumerate(wave_exp):
        found = False
        for order in reduced.orders:
            for line in order.lines:
                if abs(line.waveAccepted - exp) <= 0.1:
#                     print(str(line.acceptedWavelength) + " = " + str(exp))
                    line.frameFitOutlier = False
                    line.frameFitWave = wave_fit[i]
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
            if (line.frameFitOutlier == False):
                line.frameFitSlope = \
                        reduced.frameCalCoeffs[1] + \
                        (2.0 * reduced.frameCalCoeffs[2] * line.col) + \
                        (reduced.frameCalCoeffs[4] / order.orderNum) + \
                        (2.0 * reduced.frameCalCoeffs[5] * line.col / order.orderNum)  
    
#     print('coeffs: ' + str(reduced.coeffs))
#     print('wave_fit: ' + str(wave_fit))
#     print('wave_exp: ' + str(wave_exp))
#     raw_input('waiting')
    
    return True

def apply_wavelength_soln(reduced):
    
    if reduced.frameCalCoeffs is None:
        return
    
    for order in reduced.orders:
        x = np.arange(1024)
        y = 1.0 / order.orderNum
        order.frameCalWaveScale = np.ravel(
                reduced.frameCalCoeffs[0] +
                reduced.frameCalCoeffs[1] * x +
                reduced.frameCalCoeffs[2] * x ** 2 +
                reduced.frameCalCoeffs[3] * y +
                reduced.frameCalCoeffs[4] * x * y +
                reduced.frameCalCoeffs[5] * (x **2) * y)
        
        dx = np.diff(order.frameCalWaveScale)
        if np.all(dx <= 0) or np.all(dx >= 0):
            # wavelength scale is monotonic
            order.waveScale = order.frameCalWaveScale
            order.calMethod = 'frame sky line cal'
        else:
            logger.warning('wavelength scale not monotonic for order {}'.format(order.orderNum))
            logger.warning('using theoretical wavelength scale ')
            order.waveScale = order.gratingEqWaveScale
            order.calMethod = 'grating equation'
                   
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
    
    if config.params['cmnd_line_mode'] is True:
        fn = out_dir + '/nsdrp.log'
    else:
        fn = out_dir + '/' + objFileName[objFileName.find("NS"):].rstrip('.gz').rstrip('.fits')  + '.log'

    if config.params['subdirs'] is False and config.params['cmnd_line_mode'] is False:
        parts = fn.split('/')
        parts.insert(len(parts)-1, 'log')
        fn = '/'.join(parts)
        
    if os.path.exists(fn):
        os.rename(fn, fn + '.prev')
        
    fh = logging.FileHandler(filename=fn)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    if config.params['verbose'] is True:
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
    logger.info('nsdrp version {}'.format(nsdrp.VERSION))
    
    loggers = ['obj']
    if config.params['cmnd_line_mode'] is False:
        loggers.append('main')

    for l in loggers:
        logging.getLogger(l).log(INFO, 'starting reduction of ' + reduced.baseName)
#           #reduced.getFileName()[reduced.getFileName().rfind('/') + 1:].rstrip('.gz').rstrip('.fits'))
        logging.getLogger(l).log(INFO, '   date of observation = ' + str(reduced.getDate() + ' UT'))
        logging.getLogger(l).log(INFO, '           target name = ' + str(reduced.getTargetName()))
        logging.getLogger(l).log(INFO, '           filter name = ' + str(reduced.getFullFilterName()))
        logging.getLogger(l).log(INFO, '             slit name = ' + str(reduced.getSlit()))
        logging.getLogger(l).log(INFO, '      integration time = ' + str(reduced.getITime()) + ' sec')
        logging.getLogger(l).log(INFO, '              n coadds = ' + str(reduced.getNCoadds()))
    logger.info('        echelle angle = ' + str(reduced.getEchPos()) + ' deg')
    logger.info('cross disperser angle = ' + str(reduced.getDispPos()) + ' deg')
    return

def process_darks_and_flats(raw, reduced):
    
    # combine flats and darks
    logger.info(str(len(raw.flatFileNames)) + ' flats: ' + 
            ', '.join(str(x) for x in ([s[s.find("NS"):s.rfind(".")] for s in raw.flatFileNames])))
    reduced.flat = raw.combineFlats()
    if len(raw.flatFileNames) > 1:
        logger.info(str(len(raw.flatFileNames)) + ' flats have been median combined')

    if len(raw.darkFileNames) > 0:
        reduced.hasDark = True
        logger.info(str(len(raw.darkFileNames)) + ' darks: ' + 
            ', '.join(str(x) for x in ([s[s.find("NS"):s.rfind(".")] for s in raw.darkFileNames])))
        reduced.dark = raw.combineDarks()
        if len(raw.darkFileNames) > 1:
            logger.info(str(len(raw.darkFileNames)) + ' darks have been median combined')
    else:
        logger.info('no darks')

    # if dark(s) exist, subtract from obj and flat
    if reduced.hasDark:
        reduced.subtractDark()
        logger.info('dark subtracted from obj and flat')
        
    return
    
def process_flats(raw, reduced):
    
    # combine flats and darks
    logger.info(str(len(raw.flatFileNames)) + ' flats: ' + 
            ', '.join(str(x) for x in ([s[s.find("NS"):s.rfind(".")] for s in raw.flatFileNames])))
    reduced.flat = raw.combineFlats()
    if len(raw.flatFileNames) > 1:
        logger.info(str(len(raw.flatFileNames)) + ' flats have been median combined')
        
    return

def process_darks(raw, reduced):

    if len(raw.darkFileNames) > 0:
        reduced.hasDark = True
        logger.info(str(len(raw.darkFileNames)) + ' darks: ' + 
            ', '.join(str(x) for x in ([s[s.find("NS"):s.rfind(".")] for s in raw.darkFileNames])))
        reduced.dark = raw.combineDarks()
        if len(raw.darkFileNames) > 1:
            logger.info(str(len(raw.darkFileNames)) + ' darks have been median combined')
    else:
        logger.info('no darks')

    # if dark(s) exist, subtract from obj and flat
    if reduced.hasDark:
        reduced.subtractDark()
        logger.info('dark subtracted from obj and flat')
        
    return
    
