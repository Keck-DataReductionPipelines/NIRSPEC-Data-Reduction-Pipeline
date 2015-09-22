import logging
import numpy as np

import nirspec_lib
import Order

logger = logging.getLogger('obj')

LARGE_TILT_THRESHOLD = 20
LARGE_TILT_EXTRA_PADDING = 10
OVERSCAN_WIDTH = 10

def extract_order(order_num, obj, flat, top_calc, bot_calc, filterName, slitName):
    """
    
    filter and slit arguments are required to select filter and slit-specific 
    extraction tuning parameters
    """
    
    order = Order.Order(order_num)
    
    order.topCalc = top_calc
    order.botCalc = bot_calc
    
    params = get_extraction_params(filterName, slitName)
    
    # create top and bottom edge images for edge location and tracing
    tops, bots = make_top_and_bots(flat)
    
    # find actual top and bottom of order using edge detection
    determine_edge_locations(tops, bots, order, params['sigma'], params['thresh'])
    
    if order.topMeas is None and order.botMeas is None:
        logger.error('could not find top or bottom of order')
        return None
    
    # find order edge traces   
    find_edge_traces(tops, bots, order)
    if order.botTrace is not None:
        order.botMeas = order.botTrace[1] 

    crosscut = np.median(flat, axis=1)

#     import pylab as pl
#     pl.figure('crosscut', facecolor='white')
#     pl.cla()
#     pl.plot(crosscut, 'k-', linewidth=2)
#     pl.xlim(0, 1023)
#     pl.show()
#     
#     pl.figure('crosscut', facecolor='white')
#     pl.cla()
#     pl.plot(flat[:, 511:512], 'k-', linewidth=2)
#     pl.xlim(0, 1023)
# 
#     pl.show()
    
#     import pylab as pl
#     pl.figure('cutout', facecolor='white')
#     pl.cla()
#     pl.imshow(obj, vmin=0, vmax=512)
#     pl.plot(order.topTrace, 'w-', linewidth=2)
#     pl.plot(order.botTrace, 'w-', linewidth=2)
#     pl.xlim(0, 1023)
#     pl.ylim(1023, 0)
# 
#     pl.show()
#     
#     import pylab as pl
#     pl.figure('cutout', facecolor='white')
#     pl.cla()
#     pl.imshow(tops, vmin=0, vmax=1024)
#     pl.set_cmap('gray')
# #     pl.plot(order.topTrace, 'w-', linewidth=2)
# #     pl.plot(order.botTrace, 'w-', linewidth=2)
#     pl.colorbar()
#     pl.xlim(0, 1023)
#     pl.ylim(1023, 0)
#     pl.show()
    
    # cut out order from object frame and flat and compute on and off order masks
    cut_out_order(obj, flat, order, params['padding'])
             
    return order


def cut_out_order(obj, flat, order, padding):
    
    # add extra padding for orders with large tilt

    tilt = abs(order.avgTrace[0] - order.avgTrace[-1])
    if  tilt > LARGE_TILT_THRESHOLD:
        logger.info('large order tilt detected, tilt = ' + str(round(tilt, 1)) + 
            ' threshold = ' + str(LARGE_TILT_THRESHOLD) + 
            ' extra padding = ' + str(LARGE_TILT_EXTRA_PADDING))
        padding += LARGE_TILT_EXTRA_PADDING
    order.padding = padding
    logger.info('cutout padding = ' + str(round(padding, 0)))
    
    # determine highest point of top trace (ignore edge)
    if order.topTrace is None:
        top_trace_approx = order.botTrace + order.topMeas + 1
        highest_point = max(top_trace_approx[0], top_trace_approx[-OVERSCAN_WIDTH])
    else:
        highest_point = max(order.topTrace[0], order.topTrace[-OVERSCAN_WIDTH])
        
    # get cut outs and masks
    if order.botMeas is None:
        bot = order.botCalc
    else:
        bot = order.botMeas
        
    order.objCutout = np.array(cut_out(obj, highest_point, bot, padding))
    
    order.flatCutout = np.array(cut_out(flat, highest_point, bot, padding))
    
    order.shiftOffset = padding + order.botMeas

    if float(order.botMeas) > float(padding):
        order.onOrderMask, order.offOrderMask = get_masks(order.objCutout.shape, 
                order.topTrace - order.botTrace[0] + padding, 
                order.botTrace - order.botTrace[0] + padding)
#                 order.topTrace - order.botMeas + padding, 
#                 order.botTrace - order.botMeas + padding)
    else:
        order.onOrderMask, order.offOrderMask = get_masks(order.objCutout.shape, 
            order.topTrace, 
            order.botTrace)
    
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
    
    
def find_edge_traces(tops, bots, order):  
      
    if order.topMeas is not None:
        logger.debug('tracing top of order')
        order.topTrace = nirspec_lib.trace_order_edge(tops, order.topMeas)
        
    if order.botMeas is not None:
        logger.debug('tracing bottom of order')
        order.botTrace = nirspec_lib.trace_order_edge(bots, order.botMeas)

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

    
def find_peak(edges, start, sigma):
    
    from scipy.signal import argrelextrema

    # take a vertical cut of edges
    magcrosscut = np.median(edges[:, 40:50], axis=1)

#     import pylab as pl
#     pl.figure('crosscut', facecolor='white')
#     pl.cla()
#     pl.plot(magcrosscut, 'k-', linewidth=2)
#     pl.xlim(0, 1023)
# 
#     pl.show()

    # find the highest peaks in crosscut, search +/- 15 pixels to narrow down list
    extrema = argrelextrema(magcrosscut, np.greater, order=35)[0]

    # find crosscut values at those extrema
    magcrosscutatextrema = magcrosscut[extrema]

    # narrow down extrema list to only ones over sigma
    peaks = np.where(magcrosscutatextrema > sigma)

    actualpeaks = extrema[peaks[0]]
     
    if actualpeaks.any():
        return min((abs(start - i), i) for i in actualpeaks)[1]
    else:
        return None

       
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
    
    