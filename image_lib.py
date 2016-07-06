import numpy as np
import scipy.ndimage as ndimage
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
    
    c = cosmics.cosmicsImage(data, sigclip=sig_clip, sigfrac=sig_frac, objlim=obj_lim, 
            verbose=False)
    c.run(max_iter)
    return(c.cleanarray)


# EXT_WINDOW = 8
# SKY_WINDOW = 8
# SKY_DISTANCE = 2

def get_extraction_ranges(image_width, peak_location, obj_w, sky_w, sky_dist):
    """

    """
    if obj_w % 2:
        ext_range = np.array(range(int((1 - obj_w) / 2.0), int((obj_w + 1) / 2.0))) + peak_location
    else:  
        ext_range = np.array(range((-obj_w) / 2, obj_w / 2)) + peak_location

    sky_range_top = np.array(range(ext_range[-1] + sky_dist, 
                          ext_range[-1] + sky_dist + sky_w))
    sky_range_bot = np.array(range(ext_range[0] - sky_dist - sky_w + 1,
                          ext_range[0] - sky_dist + 1))
        
    ext_range = np.ma.masked_less(ext_range, 0).compressed()
    sky_range_top = np.ma.masked_less(sky_range_top, 0).compressed()
    sky_range_bot = np.ma.masked_less(sky_range_bot, 0).compressed()
    
    ext_range = np.ma.masked_greater_equal(ext_range, image_width).compressed()
    sky_range_top = np.ma.masked_greater_equal(sky_range_top, image_width).compressed()
    sky_range_bot = np.ma.masked_greater_equal(sky_range_bot, image_width).compressed()

#     if ((peak_location + sky_range_bot[-1] + 1) > image_width) or \
#             ((peak_location + sky_range_bot[0]) < 0):
# 
#         sky_distance = min(peak_location - sky_range_bot, image_width - peak_location)
#         sky_height = 2
# 
#         if ext_range[0] - sky_distance + 1 < image_width:
#             sky_range_bot = range(ext_range[0] - sky_distance - sky_height + 1, 
#                     ext_range[0] - sky_distance + 1)
#         else:
#             sky_range_bot = None
# 
#     if (peak_location + sky_range_top[-1] + 1) > image_width:
#         sky_distance = 1
#         sky_height = 4
#         if peak_location + 1 + 4 < image_width:
#             sky_range_top = range(ext_range[-1] + 1, ext_range[-1] + 1 + 4)
#         else:
#             sky_range_top = None
# 
#     if image_width - peak_location < 2 or peak_location < 3:
#         return None, None, None
    
    return ext_range.tolist(), sky_range_top.tolist(), sky_range_bot.tolist()


# def extract_spectra(obj, noise, peak, obj_range, sky_range_top, sky_range_bot):
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