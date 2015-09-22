import numpy as np
import pylab as pl
import scipy
from skimage.feature.peak import peak_local_max
import cosmics


def rectify_spatial(data, curve):
    """
    Shift data, column by column, along y-axis according to curve.
    
    Returns shifted image.
    
    Throws IndexError exception if length of curve 
    is not equal to number of columns in data. 
    """
    
    # shift curve to be centered at middle of order 
    # and change sign so shift is corrective
    curve_p = -1.0 * (curve - (data.shape[0] / 2))
    
    rectified = []
    for i in range(0, len(curve_p)):
        s = data[:, i]
        rectified.append(scipy.ndimage.interpolation.shift(
                s, curve_p[i], order=3, mode='nearest', prefilter=True))
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
        rectified.append(scipy.ndimage.interpolation.shift(
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
    
    c = cosmics.cosmicsImage(data, sigclip=sig_clip, sigfrac=sig_frac, objlim=obj_lim, 
            verbose=False)
    c.run(max_iter)
    return(c.cleanarray)


EXT_WINDOW = 6
SKY_WINDOW = 8
SKY_DISTANCE = 10

def get_extraction_ranges(image_width, peak_location):
    """

    """
    if EXT_WINDOW % 2:
        ext_range = range(int((1 - EXT_WINDOW) / 2.0), int((EXT_WINDOW + 1) / 2.0))
    else:  
        ext_range = range((-EXT_WINDOW) / 2, EXT_WINDOW / 2)  

    sky_range_top = range(ext_range[-1] + SKY_DISTANCE, 
                          ext_range[-1] + SKY_DISTANCE + SKY_WINDOW)
    sky_range_bot = range(ext_range[0] - SKY_DISTANCE - SKY_WINDOW + 1,
                          ext_range[0] - SKY_DISTANCE + 1) 

    if ((peak_location + sky_range_bot[-1] + 1) > image_width) or \
            ((peak_location + sky_range_bot[0]) < 0):

        sky_distance = min(peak_location - sky_range_bot, image_width - peak_location)
        sky_height = 2

        if ext_range[0] - sky_distance + 1 < image_width:
            sky_range_bot = range(ext_range[0] - sky_distance - sky_height + 1, 
                    ext_range[0] - sky_distance + 1)
        else:
            sky_range_bot = None

    if (peak_location + sky_range_top[-1] + 1) > image_width:
        sky_distance = 1
        sky_height = 4
        if peak_location + 1 + 4 < image_width:
            sky_range_top = range(ext_range[-1] + 1, ext_range[-1] + 1 + 4)
        else:
            sky_range_top = None

    if image_width - peak_location < 2 or peak_location < 3:
        return None, None, None
    
    return ext_range, sky_range_top, sky_range_bot


def extract_spectra(obj, noise, peak, obj_range, sky_range_top, sky_range_bot):
    
    """
    """
    
    obj_sum = np.sum(obj[peak - i, :] for i in obj_range)
    sky_top_sum = np.sum(obj[peak + i, :] for i in sky_range_top)
    sky_bot_sum = np.sum(obj[peak + i, :] for i in sky_range_bot)
    
    obj_mean = obj_sum / len(obj_range)
    sky_mean = (sky_top_sum + sky_bot_sum) / (len(sky_range_top) + len(sky_range_bot))
    sky_mean -= np.median(sky_mean) # why this?
    obj_sp = obj_mean - sky_mean
    
    sky_sp = sky_mean - sky_mean.mean() # why this?
    
    obj_noise_sum = np.sum(noise[peak - i, :] for i in obj_range)
    sky_noise_top_sum = np.sum(noise[peak + i, :] for i in sky_range_top)
    sky_noise_bot_sum = np.sum(noise[peak + i, :] for i in sky_range_bot)
    
    k = np.square(len(obj_range)) / np.square((len(sky_range_top) + len(sky_range_bot)))
    noise_sp = np.sqrt(obj_noise_sum + (k * (sky_noise_top_sum + sky_noise_bot_sum)))
    
    return obj_sp, sky_sp, noise_sp