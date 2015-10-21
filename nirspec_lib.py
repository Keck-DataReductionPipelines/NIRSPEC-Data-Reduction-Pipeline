import numpy as np
from scipy.signal._peak_finding import argrelextrema

import logging

import nirspec_constants as constants
import image_lib
import tracer
from __builtin__ import False

logger = logging.getLogger('obj')
      
        
def calc_noise_img(obj, flat, integration_time):
    """
    flat is expected to be normalized and both obj and flat are expected to be rectified
    """
    
    G = 5.8     
    RN = 23.0
    DC = 0.8

#     # normalize flat
#     flat_n = image_lib.normalize(flat, on_order, off_order)
#     
#     # rectify obj and flat
#     obj_r = image_lib.rectify_spatial(obj, spatial_curve)
#     obj_r = image_lib.rectify_spectral(obj_r, spectral_curve)
#     flat_r = image_lib.rectify_spatial(flat_n, spatial_curve)
#     flat_r = image_lib.rectify_spectral(flat_r, spectral_curve)
    
    # calculate photon noise
    noise = obj / G
    
    # add read noise
    noise += np.square(RN / G)
    
    # add dark current noise
    noise += (DC / G) * integration_time
    
    # divide by normalized flat squared
    noise /= np.square(flat)
    
    return noise
    

#     pl.figure("")
#     pl.plot(curve, "r", label="curve")
#     pl.plot(curve_p, "b", label="curve_p")
#     pl.legend(loc='best', prop={'size': 8})
#     pl.show()


ORDER_EDGE_SEARCH_WIDTH = 10
ORDER_EDGE_BG_WIDTH = 30
ORDER_EDGE_JUMP_THRESH = 1.9
ORDER_EDGE_JUMP_LIMIT = 200

def trace_order_edge(data, start):
    trace, nJumps =  tracer.trace_edge(
            data, start, ORDER_EDGE_SEARCH_WIDTH, ORDER_EDGE_BG_WIDTH, ORDER_EDGE_JUMP_THRESH)
    if trace is None:
        logger.warning('trace failed')
        return None
    if nJumps > ORDER_EDGE_JUMP_LIMIT:
        logger.warning('order edge trace jump limit exceeded: n jumps=' + 
                str(nJumps) + ' limit=' + str(ORDER_EDGE_JUMP_LIMIT))
        return None
    return trace
    
SKY_LINE_SEARCH_WIDTH = 3
SKY_LINE_BG_WIDTH = 0
SKY_LINE_JUMP_THRESH = 0.8
SKY_LINE_JUMP_LIMIT = 10
        
def trace_sky_line(data, start):
    trace, nJumps =  tracer.trace_edge(
            data, start, SKY_LINE_SEARCH_WIDTH, SKY_LINE_BG_WIDTH, SKY_LINE_JUMP_THRESH)
    if trace is None:
        logger.warning('sky line trace failed')
        return None
    if nJumps > SKY_LINE_JUMP_LIMIT:
        logger.debug('sky line trace jump limit exceeded: n jumps=' + 
                str(nJumps) + ' limit=' + str(SKY_LINE_JUMP_LIMIT))        
        return None
    return trace


def smooth_spatial_trace(y_raw):
    """
    """
    
    deg = 3
    n_end_ignore = 20
    threshold = 3
    
    mask = np.ones(y_raw.shape[0] - n_end_ignore, dtype=bool)
    mask = np.append(mask, np.zeros(n_end_ignore, dtype=bool))
    
    x = np.arange(y_raw.shape[0])
    
    coeffs = np.polyfit(x[mask], y_raw[mask], deg)
    y_fit = np.polyval(coeffs, x)
    res1 = y_raw - y_fit
    stdev1 = np.std(res1)
    
    greater = np.greater(np.absolute(res1), threshold * stdev1)
    mask = np.logical_and(mask, np.logical_not(greater))
    
    coeffs = np.polyfit(x[mask], y_raw[mask], deg)
    y_fit = np.polyval(coeffs, x)
    res2 = y_raw - y_fit
    stdev2 = np.std(res2)
    
    #pl.figure()
    #pl.plot(np.absolute(res1), "r-")
    #pl.plot(np.absolute(res2), "b-")
 
    return y_fit, mask

SKY_SIGMA = 2.25
EXTRA_PADDING = 5
MIN_LINE_SEPARATION = 5

def find_spectral_trace(data, padding):
    
    # transpose the array because spectroid can only read horizontal peaks for now
    npsts = data.transpose()

    # The order cutout has padding on each side. In order to find the sky lines we should 
    # only look at the central section of the cut out array
    npsts = npsts[:, padding + 5:npsts.shape[1] - 5 - padding]
    cc = np.sum(npsts[:, 0:5], axis=1)
    locpeaks = argrelextrema(cc, np.greater)     
    locmaxes = np.where(cc[locpeaks[0]] > SKY_SIGMA * cc.mean())
    maxes = np.array(locpeaks[0][locmaxes[0]])

    deletelist = []
   
    # remove adjacent sky lines that are closer than MIN_LINE_SEPARATION pixels
    for i in range(1, len(maxes)):
        if abs(maxes[i] - maxes[i - 1]) < MIN_LINE_SEPARATION:
            deletelist.append(i)
    maxes = np.delete(maxes, deletelist, None)

    peaks = cc[maxes]        

    sortorder = np.argsort(peaks)
            
    maxes = maxes[sortorder]
    maxes = maxes[::-1]

    skydict = {}
    centroid_sky_sum = np.array([])
    fitnumber = 0

    for maxskyloc in maxes:
        if 10 < maxskyloc < 1010:
            
            centroid_sky = trace_sky_line(npsts, maxskyloc)
           
            if centroid_sky is None:
                continue  # skip this skyline

            # average up the good ones
            # if badfit < 10:
            if True:
                skydict[fitnumber] = centroid_sky
                fitnumber += 1
                if centroid_sky_sum.any():
                    centroid_sky_sum = centroid_sky_sum + centroid_sky - centroid_sky[0]
                else:
                    centroid_sky_sum = centroid_sky - centroid_sky[0]
            if fitnumber > 2:
                break

    if centroid_sky_sum.any():
        logger.info(str(fitnumber) + ' sky lines used for spectral rectification')
        return centroid_sky_sum / fitnumber
    
    logger.warning('failed to find sky line trace')
    raise StandardError('failed to find sky line trace')
    
    
def smooth_spectral_trace(data, l):
    p0 = np.polyfit(np.arange(len(data) - 10), data[:-10], deg=1)  # end always drops off
    logger.info('spectral tilt is {:.3f} pixels/pixel'.format(p0[0]))
    fit = np.polyval(p0, np.arange(l))
    return fit