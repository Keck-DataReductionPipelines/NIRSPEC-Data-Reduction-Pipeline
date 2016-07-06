import logging
import numpy as np
import scipy.stats
import scipy.optimize
#import scipy.ndimage

import config
import image_lib
import nirspec_lib
import wavelength_utils
import Line
import DrpException

logger = logging.getLogger('obj')

def reduce_order(order, flat_order):
        
    # flatten obj but keep original for noise calc
    order.flattenedObjImg = np.array(order.objCutout / flat_order.normFlatImg)
    if np.amin(order.flattenedObjImg) < 0:
        order.flattenedObjImg -= np.amin(order.flattenedObjImg)
    order.flattened = True
    order.objImg = np.array(order.objCutout) # should probably use objImg instead of objCutout to begin with
    logger.info('order has been flat fielded')
 
    # rectify flat, normalized flat, obj and flattened obj in spatial dimension
    order.objImg = image_lib.rectify_spatial(order.objImg, flat_order.smoothedSpatialTrace)
    order.flattenedObjImg = image_lib.rectify_spatial(
            order.flattenedObjImg, flat_order.smoothedSpatialTrace)
 
    # trim rectified order
    order.objImg = order.objImg[flat_order.botTrim:flat_order.topTrim, :]
    order.flattenedObjImg = order.flattenedObjImg[flat_order.botTrim:flat_order.topTrim, :]
    
    order.srNormFlatImg = flat_order.rectFlatImg
    order.srFlatObjImg = order.flattenedObjImg
    
    order.spatialRectified = True
    
    # find spatial profile and peak
    find_spatial_profile_and_peak(order)
    
    # characterize spatial profile by fitting to Gaussian
    characterize_spatial_profile(order)


    # find and smooth spectral trace
    try:
        order.spectralTrace = nirspec_lib.smooth_spectral_trace(
                nirspec_lib.find_spectral_trace(
                        order.flattenedObjImg), order.flattenedObjImg.shape[0])
    except Exception as e:
        logger.warning('not rectifying order {} in spectral dimension'.format(order.orderNum))
 
    else:
        flat_order.rectFlatImg = image_lib.rectify_spectral(flat_order.rectFlatImg, order.spectralTrace)
        order.objImg = image_lib.rectify_spectral(order.objImg, order.spectralTrace)
        order.objImgFlattened = image_lib.rectify_spectral(order.flattenedObjImg, order.spectralTrace)
        order.spectralRectified = True
     
    # compute noise image
    order.noiseImg = nirspec_lib.calc_noise_img(
            order.objImg, flat_order.rectFlatImg, order.integrationTime)
    
    # extract spectra
    extract_spectra(order, flat_order)
            
    # calculate SNR
    bg = 0.0
    if order.topBgMean is not None:
        bg += order.topBgMean
    if order.botBgMean is not None:
        bg += order.botBgMean
    if order.topBgMean is not None and order.botBgMean is not None:
        bg /= 2
    order.snr = np.absolute(np.mean(order.flattenedObjImg[order.peakLocation:order.peakLocation + 1, :]) / bg)
    logger.info('signal-to-noise ratio = {:.1f}'.format(order.snr))
            
    # find and identify sky lines   
    line_pairs = None # line_pairs are (column number, accepted wavelength
    try:
        oh_wavelengths, oh_intensities = wavelength_utils.get_oh_lines()
    except IOError as e:
        logger.critical('cannot read OH line file: ' + str(e))
        raise
        
    try:
        # synthesize sky spectrum and store in order object
        order.synthesizedSkySpec = wavelength_utils.synthesize_sky(
                oh_wavelengths, oh_intensities, order.wavelengthScaleCalc)
         
        # identify lines and return list of (column number, accepted wavelength) tuples
        line_pairs = wavelength_utils.line_id(order, oh_wavelengths, oh_intensities)
        
    except (IOError, ValueError) as e:
        logger.warning('sky line matching failed: ' + str(e))
        
    if line_pairs is not None:
        
        logger.info(str(len(line_pairs)) + ' matched sky lines found in order')

        # add line pairs to Order object as Line objects
        for line_pair in line_pairs:
            col, waveAccepted = line_pair
            peak = order.skySpec[col]
            cent = image_lib.centroid(order.skySpec, 1024, 5, col)
            line = Line.Line(order.frame, order.orderNum, waveAccepted, col, cent, peak)
            order.lines.append(line)
            
        if len(order.lines) >= 3:
            # do local wavelength fit
            measured = []
            accepted = []
            for line in order.lines:
                measured.append(order.wavelengthScaleCalc[line.col])
                accepted.append(line.waveAccepted)
            (order.perOrderSlope, order.perOrderIntercept, order.perOrderCorrCoeff, p, e) = \
                    scipy.stats.linregress(np.array(measured), np.array(accepted))  
                    
            logger.info('per order wavelength fit: n = {}, a = {:.6f}, b = {:.6f}, r = {:.6f}'.format(
                    len(order.lines), order.perOrderIntercept, order.perOrderSlope, 
                    order.perOrderCorrCoeff))

            for line in order.lines:
                line.orderWaveFit = order.perOrderIntercept + \
                    (order.perOrderSlope * order.wavelengthScaleCalc[line.col])    
                line.orderFitRes = abs(line.orderWaveFit - line.waveAccepted)  
                line.orderFitSlope = (order.perOrderSlope * (order.wavelengthScaleCalc[1023] - order.wavelengthScaleCalc[0]))/1024.0
    else:
        logger.warning('no matched sky lines in order ' + str(order.orderNum))
                        
    return
         

        
    
    
def extract_spectra(order, flat_order):
    
    order.objWindow, order.topSkyWindow, order.botSkyWindow = \
        image_lib.get_extraction_ranges(order.objImg.shape[0], order.peakLocation,
        config.params['obj_window'], config.params['sky_window'], config.params['sky_separation'])
        
    logger.info('extraction window width = {}'.format(str(len(order.objWindow))))
    logger.info('top background window width = {}'.format(str(len(order.topSkyWindow))))
    if len(order.topSkyWindow) > 0:
        logger.info('top background window separation = {}'.format(
                str(order.topSkyWindow[0] - order.objWindow[-1])))
    logger.info('bottom background window width = {}'.format(str(len(order.botSkyWindow))))
    if len(order.botSkyWindow) > 0:
        logger.info('bottom background window separation = {}'.format(
                str(order.objWindow[0] - order.botSkyWindow[-1])))
    
    order.objSpec, order.flatSpec, order.skySpec, order.noiseSpec, order.topBgMean, \
            order.botBgMean = image_lib.extract_spectra(
                order.flattenedObjImg, flat_order.rectFlatImg, order.noiseImg, 
                order.objWindow, order.topSkyWindow, order.botSkyWindow)
            
    return

def characterize_spatial_profile(order):
    
    try:
        for w in range(10, 30, 10):
            logger.debug('gaussian window width = {}'.format(2 * w))
            x0 = max(0, order.peakLocation - w)
            x1 = min(len(order.spatialProfile) - 1, order.peakLocation + w)
            x = range(x1 - x0)
            order.gaussianParams, pcov = scipy.optimize.curve_fit(
                    image_lib.gaussian, x, order.spatialProfile[x0:x1] - np.amin(order.spatialProfile[x0:x1]))
            order.gaussianParams[1] += x0
            if order.gaussianParams[2] > 1.0:
                break
    except Exception as e:
        logger.warning('cannot fit spatial profile to Gaussian')
        order.gaussianParams = None
    else:
        logger.info('spatial peak width = {:.1f} pixels'.format(abs(order.gaussianParams[2])))
        
    return

def find_spatial_profile_and_peak(order):
    
    MARGIN = 5
    
    order.spatialProfile = order.flattenedObjImg.mean(axis=1)
    
    if len(order.spatialProfile) < (2 * MARGIN) + 2:
        raise DrpException.DrpException(
                'cannot find spatial profile for order {}'.format(order.orderNum))
        
    order.peakLocation = np.argmax(order.spatialProfile[MARGIN:-MARGIN]) + MARGIN
    logger.info('spatial profile peak intensity row {:d}'.format(order.peakLocation))
    p0 = order.peakLocation - (config.params['obj_window'] / 2)
    p1 = order.peakLocation + (config.params['obj_window'] / 2)
    order.centroid = (scipy.ndimage.measurements.center_of_mass(
            order.spatialProfile[p0:p1]))[0] + p0 
    logger.info('spatial profile peak centroid row {:.1f}'.format(float(order.centroid)))
    
    return