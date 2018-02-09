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
import imp

logger = logging.getLogger('obj')
# main_logger = logging.getLogger('main')
# main_logger = logging.getLogger('main')

def reduce_frame(raw, out_dir, flatCacher=None):
    """
    
    Arguments:
        raw: RawDataSet object
        out_dir: Data product root directory
        
    """
         
    # initialize per-object logger and check output directory
    if raw.isPair:
        init(raw.baseNames['AB'], out_dir)
    else:
        init(raw.baseNames['A'], out_dir)
    
    # if no flats in raw data set then fail
    if (len(raw.flatFns) < 1):
        logger.error("no flats for {}".format(raw.baseName))
        raise DrpException.DrpException('no flats');
    
    # create reduced data set
    reduced = ReducedDataSet.ReducedDataSet(raw)
    
    # read raw object image data into reduced data set object
    reduced.objImg['A'] = fits.getdata(raw.objAFn, ignore_missing_end=True)
    
    if raw.isPair:
        reduced.objImg['B'] = fits.getdata(raw.objBFn, ignore_missing_end = True)
    
    # put object summary info into per-object log
    log_start_summary(reduced)
    
    # Get fully processed flat in the form of a Flat object
    reduced.Flat = getFlat(raw, flatCacher)
    logger.info('using flat {}'.format(reduced.Flat.getBaseName()))
        
    # clean cosmic ray hits on object frame(s)
    if config.params['no_cosmic']:
        logger.info("cosmic ray rejection on object frame inhibited by command line flag")

    else:
        logger.info('cosmic ray cleaning object frame A')
        reduced.objImg['A'], cosmicMethod = image_lib.cosmic_clean(reduced.objImg['A'])
        logger.debug('cosmic ray cleaning object frame A complete')
        if reduced.isPair:
            logger.info('cosmic ray cleaning object frame B')
            reduced.objImg['B'] = image_lib.cosmic_clean(reduced.objImg['B'])
            logger.debug('cosmic ray cleaning object frame B complete')
        reduced.cosmicCleaned = True 
        logger.info(cosmicMethod)
           
    # if darks are available, combine them if there are more than one
    # and subtract from object frame(s) and flat
    process_darks(raw, reduced)
    
    # if AB pair then subtract B from A
    if reduced.isPair:
        reduced.objImg['AB'] = np.subtract(reduced.objImg['A'], reduced.objImg['B'])

    # reduce orders
    try:
        reduce_orders(reduced)
    except IOError as e:
        # might want to do something else here
        raise
    
    # find and apply wavelength solution
    imp.reload(wavelength_utils)
    if find_global_wavelength_soln(reduced) is True:
        apply_wavelength_soln(reduced)
    else:
        logger.info('not applying wavelength solution')
        for order in reduced.orders:
            order.waveScale = order.flatOrder.gratingEqWaveScale
            order.calMethod = 'grating equation'
    
    return(reduced)
 
 
def getFlat(raw, flatCacher):
    """Given a raw data set and a flat cacher, creates and returns a Flat object
    
    If there is no flat cacher then we are in command line mode where there is a single
    flat and the flats are not reused.  In this case 1. the flat image data is read from the 
    flat file specified in the raw data set, 2. unless cosmic ray rejection is inhibited
    cosmic ray artifacts are removed from the flat, and 3. a Flat object is created from 
    the flat image data.
    
    If there is a flat cacher then a possibly cached Flat object is retrieved from it.
    
    Note that in either case, the Flat object represents a fully reduced flat, with orders on
    the flat identified, traced, cut out and rectified.
    
    Args:
        raw: A RawDataSet object containing, among other things, one or more flat file names
        flatCacher: A FlatCacher object or None 
        
    Returns:
        A Flat object
    
    """
    
    if flatCacher is None:
        # in command line mode with only one flat
        if config.params['no_cosmic']:
            logger.info('cosmic ray rejection on flat inhibited by command line flag')
            flat_data = fits.PrimaryHDU.readfrom(raw.flatFns[0], ignore_missing_end=True).data
        else:
            logger.info('starting cosmic ray cleaning flat')
            flat_data, cosmicMethod = image_lib.cosmic_clean(fits.PrimaryHDU.readfrom(
                    raw.flatFns[0], ignore_missing_end=True).data)
            logger.info(cosmicMethod)
            logger.info('cosmic ray cleaning on flat complete')
        
        return(Flat.Flat(
                raw.flatFns[0], [raw.flatFns[0]],
                fits.PrimaryHDU.readfrom(raw.flatFns[0], ignore_missing_end=True).header, 
                flat_data))
    else:
        return(flatCacher.getFlat(raw.flatFns))
     
 
def reduce_orders(reduced):
    """Reduces each order in the frame.  
    
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
            
        order = Order.Order(reduced.frames, reduced.baseNames, flatOrder)
        
        order.isPair = reduced.isPair
        
        for frame in order.frames:
            order.objCutout[frame] = np.array(image_lib.cut_out(reduced.objImg[frame], 
                    flatOrder.highestPoint, flatOrder.lowestPoint, flatOrder.cutoutPadding))  
        
        order.integrationTime = reduced.getIntegrationTime() # used in noise calc
            
        try:
            
            # reduce order, i.e. rectify, extract spectra, identify sky lines
            reduce_order.reduce_order(order)
    
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
        reduced.snrMean = sum(reduced.orders[i].snr['A'] for i in range(len(reduced.orders))) / \
                len(reduced.orders)
        logging.getLogger(l).log(INFO, 'mean signal-to-noise ratio = {:.1f}'.format(
                reduced.snrMean))
    
    snr = []
    for i in range(len(reduced.orders)):
        snr.append(reduced.orders[i].snr['A'])
    
    for l in loggers:
        reduced.snrMin = np.amin(snr)
        logging.getLogger(l).log(INFO, 'minimum signal-to-noise ratio = {:.1f}'.format(
                reduced.snrMin))
    
        try:
            reduced.wMean = sum(abs(reduced.orders[i].gaussianParams['A'][2]) \
                    for i in range(len(reduced.orders))) /  len(reduced.orders)
            logging.getLogger(l).log(INFO, 'mean spatial peak width = {:.1f} pixels'.format(
                    reduced.wMean))
        except:
            logging.getLogger(l).log(logging.WARNING, 'mean spatial peak width = unknown') 
    
        w = []
        for i in range(len(reduced.orders)):
            if reduced.orders[i].gaussianParams['A'] is not None:
                w.append(reduced.orders[i].gaussianParams['A'][2])
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
    # TODO: modify twodfit() to take list of lines rather than these constructed arrays
    col = []
    centroid = []
    order_inv = []
    accepted = []
    
    for order in reduced.orders:
        for line in order.lines:
            col.append(line.col)
            centroid.append(line.centroid)
            order_inv.append(1.0 / order.flatOrder.orderNum)
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
                    line.frameFitOutlier = False
                    line.frameWaveFit = wave_fit[i]
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
                        (reduced.frameCalCoeffs[4] / order.flatOrder.orderNum) + \
                        (2.0 * reduced.frameCalCoeffs[5] * line.col / order.flatOrder.orderNum)  
    
    return True

def apply_wavelength_soln(reduced):
    
    if reduced.frameCalCoeffs is None:
        return
    
    for order in reduced.orders:
        x = np.arange(1024)
        y = 1.0 / order.flatOrder.orderNum
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
            logger.warning('wavelength scale not monotonic for order {}'.format(
                    order.flatOrder.orderNum))
            logger.warning('using theoretical wavelength scale ')
            order.waveScale = order.flatOrder.gratingEqWaveScale
            order.calMethod = 'grating equation'
                   
    return
    
def init(baseName, out_dir):
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
        #fn = out_dir + '/' + objFileName[objFileName.find("NS"):].rstrip('.gz').rstrip('.fits')  + '.log'
        fn = out_dir + '/' + baseName  + '.log'

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
    """
    """
    #TODO: structure and format this like nsdrp_cmnd.write_summary()
    
    logger.info('nsdrp version {}'.format(nsdrp.VERSION))
    
    loggers = ['obj']
    if config.params['cmnd_line_mode'] is False:
        loggers.append('main')

    for l in loggers:
        if reduced.isPair:
            baseName = reduced.baseNames['AB']
        else:
            baseName = reduced.baseNames['A']
        logging.getLogger(l).log(INFO, 'starting reduction of ' + baseName)
        logging.getLogger(l).log(INFO, 'time of observation = ' + 
                str(reduced.getDate()) + ' ' + str(reduced.getTime()) + ' UT')
        logging.getLogger(l).log(INFO, 'target name = ' + str(reduced.getTargetName()))
        logging.getLogger(l).log(INFO, 'filter name = ' + str(reduced.getFullFilterName()))
        logging.getLogger(l).log(INFO, 'slit name = ' + str(reduced.getSlit()))
        logging.getLogger(l).log(INFO, 'integration time = ' + str(reduced.getITime()) + ' sec')
        logging.getLogger(l).log(INFO, 'n coadds = ' + str(reduced.getNCoadds()))
    logger.info('echelle angle = ' + str(reduced.getEchPos()) + ' deg')
    logger.info('cross disperser angle = ' + str(reduced.getDispPos()) + ' deg')
    return

def process_darks(raw, reduced):
    """
    
    """
    
    #TODO if there are < 3 darks should cosmic clean

    if len(raw.darkFns) > 0:
        reduced.hasDark = True
        logger.info(str(len(raw.darkFns)) + ' darks: ' + 
            ', '.join(str(x) for x in ([s[s.find("NS"):s.rfind(".")] for s in raw.darkFns])))
        reduced.dark = raw.combineDarks()
        if len(raw.darkFns) > 1:
            logger.info(str(len(raw.darkFns)) + ' darks have been median combined')
    else:
        logger.info('no darks')

    # if dark(s) exist, subtract from obj and flat
    if reduced.hasDark:
        reduced.subtractDark()
        logger.info('dark subtracted from obj and flat')
        
    return
    
