import logging
import numpy as np
import scipy.stats
import scipy.ndimage

import config
import image_lib
import nirspec_lib
import wavelength_utils
import Line
#import Order

logger = logging.getLogger('obj')

def reduce_order(order):
            
    # normalize flat
    order.normalizedFlatImg, order.flatMean =  image_lib.normalize(
            order.flatCutout, order.onOrderMask, order.offOrderMask)
    order.flatNormalized = True
    logger.info('flat normalized, flat mean = ' + str(round(order.flatMean, 1)))
        
    # flatten obj but keep original for noise calc
    order.flattenedObjImg = np.array(order.objCutout / order.normalizedFlatImg)
    order.flattened = True
    order.objImg = np.array(order.objCutout) # should probably use objImg instead of objCutout to begin with
    order.flatImg = np.array(order.flatCutout)
    logger.info('order has been flat fielded')
    
    # smooth spatial trace
    # this should probably be done where the trace is first found    
    order.smoothedTrace, order.traceMask = nirspec_lib.smooth_spatial_trace(order.avgTrace)
    logger.info('spatial trace smoothed, ' + str(order.objImg.shape[1] - np.count_nonzero(order.traceMask)) + 
            ' points ignored')

    
    # rectify flat, normalized flat, obj and flattened obj in spatial dimension
    order.flatImg = image_lib.rectify_spatial(order.flatImg, order.smoothedTrace)
    order.normalizedFlatImg = image_lib.rectify_spatial(order.normalizedFlatImg, order.smoothedTrace)
    order.objImg = image_lib.rectify_spatial(order.objImg, order.smoothedTrace)
    order.flattenedObjImg = image_lib.rectify_spatial(order.flattenedObjImg, order.smoothedTrace)
 
    if order.lowestPoint > order.padding:
        top = order.highestPoint - order.lowestPoint + order.padding - 3
    else:
        top = order.highestPoint - 3
    h = np.amin(order.topTrace - order.botTrace)
    bot = top - h + 3
 

    bot = max(0, bot)
    top = min(top, 1023)
    
    order.flatImg = order.flatImg[bot:top, :]
    order.normalizedFlatImg = order.normalizedFlatImg[bot:top, :]
    order.objImg = order.objImg[bot:top, :]
    order.flattenedObjImg = order.flattenedObjImg[bot:top, :]
    
    order.srNormFlatImg = order.normalizedFlatImg
    order.srFlatObjImg = order.flattenedObjImg
    
    order.spatialRectified = True
    
    # find spatial profile and peak
    order.spatialProfile = order.flattenedObjImg.mean(axis=1)
    order.peakLocation = np.argmax(order.spatialProfile[5:-5]) + 5
    logger.info('peak intensity row {:d}'.format(order.peakLocation))
    p0 = order.peakLocation - (config.params['obj_window_width'] / 2)
    p1 = order.peakLocation + (config.params['obj_window_width'] / 2)
    order.centroid = (scipy.ndimage.measurements.center_of_mass(order.spatialProfile[p0:p1]))[0] + p0 
    logger.info('centroid row {:.1f}'.format(float(order.centroid)))
    
    
    # find and smooth spectral trace
    try:
        order.spectral_trace = nirspec_lib.smooth_spectral_trace(
                nirspec_lib.find_spectral_trace(
                        order.flattenedObjImg, order.padding), order.flattenedObjImg.shape[0])
    except Exception as e:
        logger.warning('not rectifying order in spectral dimension')
 
    else:
        # rectify flat, normalized flat, obj and flattened obj in spectral dimension 
        order.flatImg = image_lib.rectify_spectral(order.flatImg, order.spectral_trace)
        order.normalizedFlatImg = image_lib.rectify_spectral(order.normalizedFlatImg, order.spectral_trace)
        order.objImg = image_lib.rectify_spectral(order.objImg, order.spectral_trace)
        order.objImgFlattened = image_lib.rectify_spectral(order.flattenedObjImg, order.spectral_trace)
        order.spectralRectified = True
     
    # compute noise image
    order.noiseImg = nirspec_lib.calc_noise_img(
            order.objImg, order.normalizedFlatImg, order.integrationTime)
    
    # extract spectra
    order.objWindow, order.topSkyWindow, order.botSkyWindow = \
        image_lib.get_extraction_ranges(order.objImg.shape[0], order.peakLocation,
        config.params['obj_window_width'], config.params['sky_window_width'], 
        config.params['sky_dist'])
        
    logger.info('extraction window width = {}'.format(str(len(order.objWindow))))
    logger.info('top background window width = {}'.format(str(len(order.topSkyWindow))))
    logger.info('top background window separation = {}'.format(str(order.topSkyWindow[0] - order.objWindow[-1])))
    logger.info('bottom background window width = {}'.format(str(len(order.botSkyWindow))))
    logger.info('bottom background window separation = {}'.format(str(order.objWindow[0] - order.botSkyWindow[-1])))
    
    order.objSpec, order.skySpec, order.noiseSpec, order.topBgMean, order.botBgMean = \
            image_lib.extract_spectra(
                order.flattenedObjImg, order.noiseImg, 
                order.objWindow, order.topSkyWindow, order.botSkyWindow)
    
    # find and identify sky lines   
    line_pairs = None # line_pairs are (column number, accepted wavelength
    try:
        oh_wavelengths, oh_intensities = wavelength_utils.get_oh_lines(config.params['oh_filename'])
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
            line = Line.Line()
            line.col, line.acceptedWavelength = line_pair
            line.peak = order.skySpec[line.col]
            
            # find centroid of peak
            w = 5
            p0 = max(0, line.col - (w / 2))
            p1 = min(1023, line.col + (w / 2)) + 1
            line.centroid = p0 + scipy.ndimage.center_of_mass(order.skySpec[p0:p1])[0]
            
            if abs(line.centroid - line.col) > 1:
                logger.error('centroid error {} {}'.format(line.col, line.centroid))
                line.centroid = line.col                
            
            order.lines.append(line)
            
        if len(order.lines) >= 3:
            # do local wavelength fit
            measured = []
            accepted = []
            for line in order.lines:
                measured.append(order.wavelengthScaleCalc[line.col])
                accepted.append(line.acceptedWavelength)
            (order.perOrderSlope, order.perOrderIntercept, order.perOrderCorrCoeff, p, e) = \
                    scipy.stats.linregress(np.array(measured), np.array(accepted))  
                    
            logger.info('per order wavelength fit: n={}, a={:.6f}, b={:.6f}, r={:.6f}'.format(
                    len(order.lines), order.perOrderIntercept, order.perOrderSlope, 
                    order.perOrderCorrCoeff))

            for line in order.lines:
                line.localFitWavelength = order.perOrderIntercept + \
                    (order.perOrderSlope * order.wavelengthScaleCalc[line.col])    
                line.localFitResidual = abs(line.localFitWavelength - line.acceptedWavelength)  
                line.localFitSlope = (order.perOrderSlope * (order.wavelengthScaleCalc[1023] - order.wavelengthScaleCalc[0]))/1024.0
    else:
        logger.warning('no matched sky lines in order ' + str(order.orderNum))
                        
    return
         
    