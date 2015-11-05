import logging
import numpy as np

import nirspec_lib
import Order
import nirspec_constants as constants

logger = logging.getLogger('obj')

LARGE_TILT_THRESHOLD = 20
LARGE_TILT_EXTRA_PADDING = 10
OVERSCAN_WIDTH = 10

def extract_order(order_num, obj, flat, top_calc, bot_calc, filter_name, slit_name):
    """
    
    filter and slit arguments are required to select filter and slit-specific 
    extraction tuning parameters
    """
    params = get_extraction_params(filter_name, slit_name)

    order = Order.Order(order_num)
    
    order.topCalc = top_calc
    order.botCalc = bot_calc
    
    order.padding = params['padding']
    
    # create top and bottom edge images for edge location and tracing
    tops, bots = make_top_and_bots(flat)
    
    # find actual top and bottom of order using edge detection
    determine_edge_locations(tops, bots, order, params['sigma'], params['thresh'])
    
    if order.topMeas is None and order.botMeas is None:
        msg = 'could not find top or bottom of order'
        logger.debug(msg)
#         raise DrpException.DrpException(msg)
        return None
    
    # find order edge traces   
    find_edge_traces(tops, bots, order, slit_name)
    
    if order.topTrace is None and order.botTrace is None:
        logger.info('could not trace top or bottom of order edge on flat')
        return None
    
    if order.botTrace is not None:
        order.botMeas = order.botTrace[1] 
        
    if slit_name.endswith('24'):
        logger.info('applying long slit edge margins of {} pixels'.format(
                constants.LONG_SLIT_EDGE_MARGIN))
        if order.topTrace is not None:
            order.topTrace -= constants.LONG_SLIT_EDGE_MARGIN
        if order.botTrace is not None:
            order.botTrace += constants.LONG_SLIT_EDGE_MARGIN
    
    # cut out order from object frame and flat and compute on and off order masks
    cut_out_order(obj, flat, order)
             
    return order


def cut_out_order(obj, flat, order):
    """
    
    obj - the full frame object image
    flat - the full flat, possibly combined, flat image
    order - the Order object in which the cutouts will be stored
    padding - the amount of padding to use to account for curvature
    """
    
    # add extra padding for orders with large tilt

    tilt = abs(order.avgTrace[0] - order.avgTrace[-1])
    if  tilt > LARGE_TILT_THRESHOLD:
        logger.info('large order tilt detected, tilt = ' + str(round(tilt, 1)) + 
            ' threshold = ' + str(LARGE_TILT_THRESHOLD) + 
            ' extra padding = ' + str(LARGE_TILT_EXTRA_PADDING))
        order.padding += LARGE_TILT_EXTRA_PADDING
    logger.debug('cutout padding = ' + str(round(order.padding, 0)))
    
    # determine highest point of top trace (ignore edge)
    if order.topTrace is None:
        order.topTrace = order.botTrace + (order.topCalc - order.botCalc) - 5
        
    order.highestPoint = np.amax(order.topTrace[0:-OVERSCAN_WIDTH])
        
    if order.botTrace is None:
        order.botTrace = order.topTrace - (order.topCalc - order.botCalc) + 5
        
    order.lowestPoint = np.amin(order.botTrace[0:-OVERSCAN_WIDTH])
         
    order.objCutout = np.array(cut_out(obj, order.highestPoint, order.lowestPoint, order.padding))
    order.flatCutout = np.array(cut_out(flat, order.highestPoint, order.lowestPoint, order.padding))
    order.shiftOffset = order.padding + order.botMeas
    
    if float(order.lowestPoint) > float(order.padding):
        order.onOrderMask, order.offOrderMask = get_masks(
                order.objCutout.shape, 
                order.topTrace - order.lowestPoint + order.padding, 
                order.botTrace - order.lowestPoint + order.padding)
    else:
        order.onOrderMask, order.offOrderMask = get_masks(
                order.objCutout.shape, order.topTrace, order.botTrace)
        
    order.objCutout = np.ma.masked_array(order.objCutout, mask=order.offOrderMask)
    order.flatCutout = np.ma.masked_array(order.flatCutout, mask=order.offOrderMask)

    return
    
    
def get_masks(shape, top_trace, bot_trace):
    
    y, x = np.indices(shape, dtype=np.float32)
    
    off_top = y > top_trace
    off_bot = y < bot_trace
    off_order = off_top | off_bot

    belowtop = y < top_trace
    abovebot = y > bot_trace
    on_order = belowtop & abovebot
        
    return on_order, off_order
    
    
def cut_out(data, top, bot, padding):
        
    if bot > padding:
        return data[bot - padding:top + padding, :]
    else:
        return data[0:top + padding, :]
    
    
def find_edge_traces(tops, bots, order, slit_name):  
    
    if order.topMeas is not None:
        logger.debug('tracing top of order')
        order.topTrace = nirspec_lib.trace_order_edge(tops, order.topMeas, slit_name)
        
    if order.botMeas is not None:
        logger.debug('tracing bottom of order')
        order.botTrace = nirspec_lib.trace_order_edge(bots, order.botMeas, slit_name)
        
    if order.topTrace is None and order.botTrace is None:
        return

    if order.topTrace is not None and order.botTrace is not None:
        logger.info('using top and bottom trace')
        order.avgTrace = (order.topTrace + order.botTrace) / 2.0

    elif order.botTrace is None:
        logger.info('using top trace only')
        order.avgTrace = order.topTrace - ((order.topMeas - order.botCalc) / 2.0) + 1.0
        
    else:
        logger.info('using bottom trace only')
        order.avgTrace = order.botTrace + ((order.topCalc - order.botMeas) / 2.0) + 1.0
        
    return
    
    
def determine_edge_locations(tops, bots, order, sigma, thresh):
    
    # find top edge
    
    order.topMeas = find_peak(tops, order.topCalc, sigma)
    
    if order.topMeas is None or abs(order.topMeas - order.topCalc) > thresh:
        logger.info('reducing edge detection threshold')
        order.topMeas = find_peak(tops, order.topCalc, sigma / 2)
        
    if order.topMeas is not None:
        if (order.topMeas < 1) or (abs(order.topMeas - order.topCalc) > (2 * thresh)):
            s = 'top edge too far off: meas={:.0f}, diff={:.0f}, thresh={:.0f}'.format(
                    order.topMeas, abs(order.topMeas - order.topCalc), thresh)
            logger.warning(s)
            order.topMeas = None
            
    # find bottom edge
    
    order.botMeas = find_peak(bots, order.botCalc, sigma)
    
    if order.botMeas is None or abs(order.botMeas - order.botCalc) > thresh:
        order.botMeas = find_peak(bots, order.botCalc, sigma / 2) 
        
    if order.botMeas is not None:
        if (order.botMeas < 1) or (abs(order.botMeas - order.botCalc) > (2 * thresh)):
            s = 'bottom edge too far off: meas={:.0f}, diff={:.0f}, thresh={:.0f}'.format(
                    order.botMeas, abs(order.botMeas - order.botCalc), thresh)
            logger.warning(s)
            order.botMeas = None
        
    # log results
    
    if order.topMeas is None:
        top_str = 'not found'
    else:
        top_str = str(order.topMeas)

    if order.botMeas is None:
        bot_str = 'not found'
    else:
        bot_str = str(order.botMeas)

    logger.info('measured y location:  top = ' + top_str + ', bottom = ' + bot_str)
    
    return

    
def find_peak(edges, row_approximate , threshold):
    """
    Takes an order edge image, formed by shifting and subtracting a flat field image, 
    and finds the edge closest row_approximate.  Order edges appear as intensity peaks in
    the shifted and subtracted image.  row_approximate is the approximate 
    location of the order edge found by evaluation of the grating equation.  
    
    edges is a 2-d image of order edges generated by shifting and subtracting
    a flat field image
    
    row_approximate is the approximate location of the desired edge
    
    threshold - peaks of intensity less than threshold are ignored 
    """
    from scipy.signal import argrelextrema

    # find the edge profile along the y axis near the short wavelength 
    # (x -> 0) edge of the detector.  Indices in the profile array 
    # correspond to detector row numbers.
    profile = np.median(edges[:, 40:50], axis=1)

    # find the indices in profile of maxima.  Order is the number of points
    # on each side to use for the comparison.
    peak_rows = argrelextrema(profile, np.greater, order=35)[0]

    # find peak intensities at extrema
    peak_intensities = profile[peak_rows]
    
    # find the indices in peak_intensities of peaks with 
    # intensities greater than threshold
    tall_peaks_i = np.where(peak_intensities > threshold)

    # narrow peak_rows to those corresponding to tall peaks
    peak_rows = peak_rows[tall_peaks_i[0]]
     
    if peak_rows.any():
        # find the tall peak closest to row_approximate
        return min((abs(row_approximate - i), i) for i in peak_rows)[1]
    else:
        return None

# def find_peak(edges, start, sigma):
#     
#     from scipy.signal import argrelextrema
# 
#     # take a vertical cut of edges
#     magcrosscut = np.median(edges[:, 40:50], axis=1)
# 
#     # find the highest peaks in crosscut, search +/- 15 pixels to narrow down list
#     extrema = argrelextrema(magcrosscut, np.greater, order=35)[0]
# 
#     # find crosscut values at those extrema
#     magcrosscutatextrema = magcrosscut[extrema]
# 
#     # narrow down extrema list to only ones over sigma
#     peaks = np.where(magcrosscutatextrema > sigma)
# 
#     actualpeaks = extrema[peaks[0]]
#      
#     if actualpeaks.any():
#         return min((abs(start - i), i) for i in actualpeaks)[1]
#     else:
#         return None
       
def make_top_and_bots(data):
    rolled = np.roll(data, 5, axis=0)
    return rolled - data, data - rolled


def get_extraction_params(filterName, slitName):
    
    if '0.288x24' in slitName:
        pad_mod = 0.0
    else:
        pad_mod = 1.0
            
    if 'NIRSPEC-1' in filterName.upper():
        params = {'padding': 0 * pad_mod, 
                  'sigma': 300.0, 
                  'thresh': 50.0, 
                  'spw': 5.0, 
                  'trace_width': 1.5}
    
    elif 'NIRSPEC-2' in filterName.upper():
        params = {'padding': 0 * pad_mod, 
                  'sigma': 300.0, 
                  'thresh': 50.0, 
                  'spw': 5.0, 
                  'trace_width': 1.5}

    elif 'NIRSPEC-3' in filterName.upper():
        params = {'padding': 10 * pad_mod, 
                  'sigma': 300.0, 
                  'thresh': 50.0, 
                  'spw': 3.0, 
                  'trace_width': 1.1}
            
    elif 'NIRSPEC-4' in filterName.upper():
        params = {'padding': 10 + 5 * pad_mod, 
                  'sigma': 600.0, 
                  'thresh': 20.0, 
                  'spw': 3.0, 
                  'trace_width': 1.1}
            
    elif 'NIRSPEC-5' in filterName.upper():
        params = {'padding': 10 + 5 * pad_mod, 
                  'sigma': 600.0, 
                  'thresh': 50.0, 
                  'spw': 3.0, 
                  'trace_width': 1.1}
            
    elif 'NIRSPEC-6' in filterName.upper():
        params = {'padding': 15.0, 
                  'sigma': 500.0, 
                  'thresh': 20.0, 
                  'spw': 3.0, 
                  'trace_width': 1.1}
            
    elif 'NIRSPEC-7' in filterName.upper():
        params = {'padding': 30, 
                  'sigma': 100.0,
                  'thresh': 20.0, 
                  'spw': 3.0, 
                  'trace_width': 1.1}
        
    if '24' in slitName:
        params['thresh'] = 30.0
        
    return params
    
    