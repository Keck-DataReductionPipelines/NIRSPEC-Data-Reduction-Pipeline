import numpy as np
import scipy.ndimage as ndimage
#from skimage.feature.peak import peak_local_max
import cosmics
#from __builtin__ import None
import os
import logging


logger = logging.getLogger('main')

def rectify_spatial(data, curve):
    """
    Shift data, column by column, along y-axis according to curve.
    
    Returns shifted image.
    
    Throws IndexError exception if length of curve 
    is not equal to number of columns in data. 
    """
    
    # shift curve to be centered at middle of order 
    # and change sign so shift is corrective
#     curve_p = -1.0 * (curve - (data.shape[0] / 2))
    curve_p = -1.0 * curve
    curve_p = curve_p - np.amin(curve_p)
    
#     import pylab as pl
#     pl.figure()
#     pl.cla()
#     pl.plot(curve, 'r-')
#     pl.plot(curve_p, 'g-')
#     pl.show()
#     
#     pl.figure()
#     pl.cla()
#     pl.imshow(data, vmin=0, vmax=256)
#     pl.show()
    
    rectified = []
    for i in range(0, len(curve_p)):
        s = data[:, i]
        rectified.append(ndimage.interpolation.shift(
                s, curve_p[i], order=3, mode='nearest', prefilter=True))
        
#     pl.figure()
#     pl.cla()
#     pl.imshow((np.array(rectified)).transpose(), vmin=0, vmax=256)
#     pl.show()
    
    return((np.array(rectified)).transpose())


def rectify_spectral(data, curve):
    """
    Shift data, row by row, along x-axis according to curve.
    
    Returns shifted image.
    
    Throws IndexError exception if length of curve
    is not equal to number of rows in data.
    """
    
    # pivot curve around peak 
    # and change sign so shift is corrective
    profile = data.sum(axis=1)
    peak = np.argmax(profile)
    curve_p = -1.0 * (curve - curve[peak])
        
    rectified = []
    for i in range(0, len(curve_p)):
        s = data[i, :]
        rectified.append(ndimage.interpolation.shift(
                s, curve_p[i], order=3, mode='nearest', prefilter=True))  
    
    return(np.array(rectified))


def normalize(data, on_order, off_order):
    """
    data is the image cut-out plus padding
    
    on_order is array of same size as data with 
    on-order pixels set to 1.0 and off order (padding) pixels set to 0.0.
    
    off_order is array of same size as data with 
    off-order (padding) pixels set to 1.0 and on order pixels sto to 0.0.
    
    returns normalized data array and mean of the on-order pixels
    """
    
    m = np.mean(data)
    non = np.count_nonzero(on_order)
    noff = np.count_nonzero(off_order)
    
    data_copy = data
    
    if np.count_nonzero(on_order) == 0:
        return

    # ignore pixels beyond column 1000 by setting value to 1.0
    data_copy[:, 1000:] = 1.0

    # take mean of only the on-order pixels
    mean = np.ma.masked_array(data_copy, mask=off_order).mean()
    
    # create normalized data array
    normalized = (data_copy * on_order) / mean

    # around the edges of the order can blow up when div by mean, set those to one
    normalized[np.where(normalized > 10.0)] = 1.0

    # avoid zeroes (not to sure about these)
    normalized[np.where(normalized == 0.0)] = 1.0
    normalized[np.where(normalized < 0.2)] = 1.0

    return normalized, mean


def cosmic_clean(data):
    """
    """
    max_iter = 3
    sig_clip = 5.0
    sig_frac = 0.3
    obj_lim = 5.0
    
    try:
        import astroscrappy
        msg = 'using Cython implementation of LACOSMIC'
        cosmicRayMask, cleanedArray = astroscrappy.detect_cosmics(data, pssl=0.0, gain=2.2, sigclip=sig_clip, sigfrac=sig_frac, objlim=obj_lim, readnoise=10.0, satlevel=np.inf, inmask=None, sepmed=False, cleantype='medmask', fsmode='median')
    except ImportError:
        msg = 'using original LACOSMIC'
        c = cosmics.cosmicsImage(data, sigclip=sig_clip, sigfrac=sig_frac, objlim=obj_lim, verbose=False)
        c.run(max_iter)
        cleanedArray = c.cleanarray

    logger.info(msg)

    return(cleanedArray)

def get_extraction_ranges(image_width, peak_location, obj_w, sky_w, sky_dist):
    """
    This function was modified so it can be used to define object window only or object and sky 
    windows.  If sky_w and sky_dist are None then only image window pixel list is computed and
    returned.
    
    Truncate windows that extend beyond top or bottom of order.
    
    Args:
        image_width:
        peak_location:
        obj_w:
        sky_w:
        sky_dist:
        
    Returns:
        three element tuple consisting of:
            0: Extraction range list.
            1: Top sky range list or None.
            2: Bottom sky range list or None.
    """
    
    if obj_w % 2:
        ext_range = np.array(range(int((1 - obj_w) / 2.0), int((obj_w + 1) / 2.0))) + peak_location
    else:  
        ext_range = np.array(range((-obj_w) / 2, obj_w / 2)) + peak_location
    ext_range = np.ma.masked_less(ext_range, 0).compressed()
    ext_range = np.ma.masked_greater_equal(ext_range, image_width).compressed()
    ext_range_list = ext_range.tolist()


    if sky_w is not None and sky_dist is not None:
        sky_range_top = np.array(range(ext_range[-1] + sky_dist, ext_range[-1] + sky_dist + sky_w))
        sky_range_top = np.ma.masked_less(sky_range_top, 0).compressed()
        sky_range_top = np.ma.masked_greater_equal(sky_range_top, image_width).compressed()
        sky_range_top_list = sky_range_top.tolist()
    
        sky_range_bot = np.array(range(ext_range[0] - sky_dist - sky_w + 1,
                ext_range[0] - sky_dist + 1))
        sky_range_bot = np.ma.masked_less(sky_range_bot, 0).compressed()
        sky_range_bot = np.ma.masked_greater_equal(sky_range_bot, image_width).compressed()
        sky_range_bot_list = sky_range_bot.tolist()
    else:
        sky_range_top_list = None
        sky_range_bot_list = None
    
    return ext_range_list, sky_range_top_list, sky_range_bot_list


   
def extract_spectra(obj, flat, noise, obj_range, sky_range_top, sky_range_bot):
    
    """
    """
    
    obj_sum = np.sum(obj[i, :] for i in obj_range)
    flat_sum = np.sum(flat[i, :] for i in obj_range)
    
    flat_sp = flat_sum / len(obj_range)

    sky_top_sum = np.sum(obj[i, :] for i in sky_range_top)
    sky_bot_sum = np.sum(obj[i, :] for i in sky_range_bot)
    
    if len(sky_range_top) > 0:
        top_bg_mean = (sky_top_sum / len(sky_range_top)).mean()
    else:
        top_bg_mean = None
    if len(sky_range_bot) > 0:
        bot_bg_mean = (sky_bot_sum / len(sky_range_bot)).mean()
    else:
        bot_bg_mean = None
    
    sky_mean = (sky_top_sum + sky_bot_sum) / (len(sky_range_top) + len(sky_range_bot))

#     sky_mean -= np.median(sky_mean) 

    obj_sp = obj_sum - (len(obj_range) * sky_mean)

    sky_sp = sky_mean - sky_mean.mean() # why this?
    
    obj_noise_sum = np.sum(noise[i, :] for i in obj_range)
    sky_noise_top_sum = np.sum(noise[i, :] for i in sky_range_top)
    sky_noise_bot_sum = np.sum(noise[i, :] for i in sky_range_bot)
    
    k = np.square(len(obj_range)) / np.square((len(sky_range_top) + len(sky_range_bot)))
    noise_sp = np.sqrt(obj_noise_sum + (k * (sky_noise_top_sum + sky_noise_bot_sum)))
    
    return obj_sp, flat_sp, sky_sp, noise_sp, top_bg_mean, bot_bg_mean

def gaussian(x, a, b, c):
    return(a * np.exp(-(x - b)**2 / c**2))

def cut_out(data, top, bot, padding):
    top = int(top)
    bot = int(bot)
    padding = int(padding)
    
    try:    
        if bot > padding:
            return data[bot - padding:top + padding, :]
        else:
            return data[0:top + padding, :]
    except TypeError:
        if bot > padding:
            return data[bot - padding:top + padding, :]
        else:
            return data[0:top + padding, :]

def centroid(spec, width, window, approx):
    p0 = max(0, approx - (window / 2))
    p1 = min(width - 1, approx + (window / 2)) + 1
    c = p0 + ndimage.center_of_mass(spec[p0:p1])[0]
    
    if abs(c - approx) > 1:
        #logger.debug('centroid error, approx = {}, centroid = {:.3f}'.format(approx, c))
        return(approx)    
    
    return(c)            
