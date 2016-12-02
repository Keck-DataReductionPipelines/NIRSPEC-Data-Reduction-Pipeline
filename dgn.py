import matplotlib

#matplotlib.use('Agg')
    
import pylab as pl
import logging
import os
import errno
import numpy as np
from skimage import exposure
import config

import scipy.misc
from scipy.misc.pilutil import imresize

# import Order
# import ReducedDataSet

UNDERLINE = False

logger = logging.getLogger('obj')

subdirs = dict([
                ('traces.png',          'traces'),
                ('trace.npy',           'traces'),
                ('cutouts.png',         'cutouts'),
                ('obj_cutout.npy',      'cutouts'),
                ('flat_cutout.npy',     'cutouts'),
                ('spatrect.png',        'spatrect'),
                ('specrect.png',        'specrect'),
                ('edges.png',           'edges'),
                ('top_bot_edges.png',   'edges'),
                ('localcalids.txt',     'localcal'),
                ('skylines.png',        'skylines'),
                ('skylines.txt',        'skylines'),
                ('wavelength_scale.png','wavelength')
                ])

def gen(reduced, out_dir):
    

    logger.info('generating diagnostic data products for {}...'.format(reduced.getBaseName()))
    
    # make sub directories
    for k, v in subdirs.iteritems():
        try:
            os.makedirs(out_dir + '/diagnostics/' + v)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                logger.critical('failed to create output directory ' + v)
                logger.critical('failed to create output directory ' + v)
                raise IOError('failed to create output directory ' + v)
            
    
    # construct arrays for per-order wavelength table
    order_num = []
    col = []
    centroid = []
    source = []
    wave_exp = []
    wave_fit = []
    res = []
    peak = []
    slope = []
    
    for order in reduced.orders:
        if order.orderCalSlope > 0.95 and order.orderCalSlope < 1.05:
            for line in order.lines:
                order_num.append(order.flatOrder.orderNum)
                col.append(line.col)
                centroid.append(line.centroid)
                source.append('sky')
                wave_exp.append(line.waveAccepted)
                wave_fit.append(line.orderWaveFit)
                res.append(line.orderFitRes)
                peak.append(line.peak)
                slope.append(line.orderFitSlope)
                
    # per-order wavelength table
    perOrderWavelengthCalAsciiTable(
            out_dir, reduced.baseNames['A'], order_num, col, centroid, source, wave_exp, wave_fit, 
            res, peak, slope)
    
    # per-frame order edge profiles
    edges_plot(out_dir, reduced.Flat.baseName, reduced.Flat.topEdgeProfile, 
            reduced.Flat.botEdgeProfile, reduced.Flat.topEdgePeaks, reduced.Flat.botEdgePeaks)
    
    # per-frame order edge ridge images
    tops_bots_plot(out_dir, reduced.Flat.baseName, reduced.Flat.topEdgeImg, reduced.Flat.botEdgeImg)
    
    # per-frame order edge traces and order ID
    order_location_plot(out_dir, reduced.baseNames['A'], reduced.Flat.baseName, 
            reduced.Flat.flatImg, reduced.objImg['A'], reduced.orders)
    
    for order in reduced.orders:
        
        traces_plot(out_dir, reduced.baseNames['A'], order.flatOrder.flatBaseName, 
                order.flatOrder.orderNum, reduced.objImg['A'], reduced.Flat.flatImg, 
                order.flatOrder.topEdgeTrace, order.flatOrder.botEdgeTrace)
        
        # save smoothed trace cutout to numpy text file
        if config.params['npy'] is True:
            fn = constructFileName(out_dir, reduced.getBaseName(), order.flatOrder.orderNum, 'trace.npy')
            np.savetxt(fn, order.smoothedTrace)
        
        #
        # object and flat order cutout plot
        #
        if order.isPair:
            cutout = order.objCutout['AB']
        else:
            cutout = order.objCutout['A']
        if order.flatOrder.lowestPoint > order.flatOrder.cutoutPadding:
            cutouts_plot(out_dir, reduced.getBaseName(), reduced.Flat.baseName, 
                    order.flatOrder.orderNum, cutout, order.flatOrder.cutout, 
                    order.flatOrder.topEdgeTrace - order.flatOrder.lowestPoint + 
                            order.flatOrder.cutoutPadding, 
                    order.flatOrder.botEdgeTrace - order.flatOrder.lowestPoint + 
                            order.flatOrder.cutoutPadding,
                    order.flatOrder.smoothedSpatialTrace - order.flatOrder.lowestPoint + 
                            order.flatOrder.cutoutPadding)
        else:
            cutouts_plot(out_dir, reduced.getBaseName(), reduced.Flat.baseName, 
                    order.flatOrder.orderNum, cutout, order.flatOrder.cutout, 
                    order.flatOrder.topEdgeTrace, order.flatOrder.botEdgeTrace, 
                    order.flatOrder.smoothedSpatialTrace)
            
        #
        # object and flat spatially rectified order plot
        #
        if order.isPair:
            img = order.srFfObjImg['AB']
        else:
            img = order.srFfObjImg['A']
        spatrect_plot(out_dir, reduced.getBaseName(), order.flatOrder.orderNum, img, order.srNormFlatImg)
        
#         specrect_plot(out_dir, reduced.baseName, order.orderNum, 
#                 order.srFlatObjAImg, order.flattenedObjAImg)
       
        #
        # sky lines plot and table
        #
        skyLinesPlot(out_dir, order)
        skyLinesAsciiTable(out_dir, reduced.baseNames['A'], order)
        
        wavelengthScalePlot(out_dir, reduced.baseNames['A'], order)
        
    logger.info('done generating diagnostic data products')
    
    
def constructFileName(outpath, base_name, order, fn_suffix):
    if order is None:
        return outpath + '/diagnostics/' + subdirs[fn_suffix] + '/' + base_name + '_' + fn_suffix
    else:
        #return fn[:fn.rfind('_') + 1] + str(order) + fn[fn.rfind('_'):]    
        return outpath + '/diagnostics/' + subdirs[fn_suffix] + '/' + base_name + \
            '_' + str(order) + '_' + fn_suffix


def edges_plot(outpath, base_name, top_profile, bot_profile, top_peaks, bot_peaks):
    
    pl.figure('order edge profiles', facecolor='white', figsize=(6, 8))
    pl.cla()
    pl.suptitle('order edge profiles, {}'.format(base_name), fontsize=14)
    
    tops_plot = pl.subplot(2, 1, 1)
    tops_plot.set_title('tops')
    tops_plot.plot(top_profile, 'k-', linewidth=1.0)
    for peak in top_peaks:
        pl.annotate(str(peak), (peak, top_profile[peak]), size=8)
    tops_plot.set_xlim([0, 1023])

    bots_plot = pl.subplot(2, 1, 2)
    bots_plot.set_title('bottoms')
    bots_plot.plot(bot_profile, 'k-', linewidth=1.0)
    for peak in bot_peaks:
        pl.annotate(str(peak), (peak, bot_profile[peak]), size=8)
    bots_plot.set_xlim([0, 1023])

    pl.savefig(constructFileName(outpath, base_name, None, 'edges.png'))
    pl.close()    
    
    
def tops_bots_plot(outpath, base_name, tops, bots):
    
    pl.figure('edges', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('top and bottom order edges, {}'.format(base_name), fontsize=14)
#     pl.set_cmap('Blues_r')
    pl.rcParams['ytick.labelsize'] = 8

    obj_plot = pl.subplot(1, 2, 1)
    try:
        obj_plot.imshow(exposure.equalize_hist(tops))
    except:
        obj_plot.imshow(tops)

    obj_plot.set_title('top edges')
    obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])

    
    flat_plot = pl.subplot(1, 2, 2)
    try:
        flat_plot.imshow(exposure.equalize_hist(bots))
    except:
        flat_plot.imshow(bots)
    flat_plot.set_title('bottom edges')
    flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, None, 'top_bot_edges.png'))
    pl.close()
    
def traces_plot(outpath, obj_base_name, flat_base_name, order_num, obj_img, flat_img, top_trace, 
            bot_trace):
    
    pl.figure('traces', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('order edge traces, {}, order {}'.format(obj_base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')
    pl.rcParams['ytick.labelsize'] = 8

    obj_plot = pl.subplot(1, 2, 1)
    try:
        obj_plot.imshow(exposure.equalize_hist(obj_img))
    except:
        obj_plot.imshow(obj_img)
    obj_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)

    obj_plot.set_title('object ' + obj_base_name)
    obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])

    
    flat_plot = pl.subplot(1, 2, 2)
    try:
        flat_plot.imshow(exposure.equalize_hist(flat_img))
    except:
        flat_plot.imshow(flat_img)
    flat_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    flat_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)    
    flat_plot.set_title('flat ' + flat_base_name)
    flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, obj_base_name, order_num, 'traces.png'))
    pl.close()
    
def order_location_plot(outpath, obj_base_name, flat_base_name, flat_img, obj_img, orders):
    
    pl.figure('orders', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('order location and identification', fontsize=14)
    pl.set_cmap('Blues_r')
    pl.rcParams['ytick.labelsize'] = 8

    obj_plot = pl.subplot(1, 2, 1)
    try:
        obj_plot.imshow(exposure.equalize_hist(obj_img))
    except:
        obj_plot.imshow(obj_img)
    obj_plot.set_title('object ' + obj_base_name)
    obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])

    
    flat_plot = pl.subplot(1, 2, 2)
    try:
        flat_plot.imshow(exposure.equalize_hist(flat_img))
    except:
        flat_plot.imshow(flat_img)
    flat_plot.set_title('flat ' + flat_base_name)
    flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
    
    for order in orders:
        obj_plot.plot(np.arange(1024), order.flatOrder.topEdgeTrace, 'k-', linewidth=1.0)
        obj_plot.plot(np.arange(1024), order.flatOrder.botEdgeTrace, 'k-', linewidth=1.0)
        obj_plot.plot(np.arange(1024), order.flatOrder.smoothedSpatialTrace, 'y-', linewidth=1.0)
        obj_plot.text(10, order.flatOrder.topEdgeTrace[0] - 10, str(order.flatOrder.orderNum), 
                fontsize=10)
        
        flat_plot.plot(np.arange(1024), order.flatOrder.topEdgeTrace, 'k-', linewidth=1.0)
        flat_plot.plot(np.arange(1024), order.flatOrder.botEdgeTrace, 'k-', linewidth=1.0)  
        flat_plot.plot(np.arange(1024), order.flatOrder.smoothedSpatialTrace, 'y-', linewidth=1.0)  
        flat_plot.text(10, order.flatOrder.topEdgeTrace[0] - 10, str(order.flatOrder.orderNum), 
                fontsize=10)

    pl.tight_layout()
    pl.savefig(constructFileName(outpath, obj_base_name, None, 'traces.png'))
    pl.close()
    
def cutouts_plot(outpath, obj_base_name, flat_base_name, order_num, obj_img, flat_img, 
            top_trace, bot_trace, trace):
    
    pl.figure('traces', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('order cutouts, {}, order {}'.format(obj_base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')

    obj_plot = pl.subplot(2, 1, 1)
    try:
        obj_plot.imshow(exposure.equalize_hist(obj_img))
    except:
        obj_plot.imshow(obj_img)
    obj_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), trace, 'y-', linewidth=1.5)
    obj_plot.set_title('object ' + obj_base_name)
    obj_plot.set_xlim([0, 1023])
    
    flat_plot = pl.subplot(2, 1, 2)
    try:
        flat_plot.imshow(exposure.equalize_hist(flat_img))
    except:
        flat_plot.imshow(flat_img)
    flat_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    flat_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)    
    flat_plot.plot(np.arange(1024), trace, 'y-', linewidth=1.5)    
    flat_plot.set_title('flat ' + flat_base_name)
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, obj_base_name, order_num, 'cutouts.png'))
    pl.close()
    
def spatrect_plot(outpath, base_name, order_num, obj, flat):

    pl.figure('spatially rectified', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('spatially rectified, {}, order {}'.format(base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')

    obj_plot = pl.subplot(2, 1, 1)
    try:
        obj_plot.imshow(exposure.equalize_hist(obj))
    except:
        obj_plot.imshow(obj)
    obj_plot.set_title('object')
#     obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])
    
    flat_plot = pl.subplot(2, 1, 2)
    try:
        flat_plot.imshow(exposure.equalize_hist(flat))
    except:
        flat_plot.imshow(flat)
    flat_plot.set_title('flat')
#     flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, order_num, 'spatrect.png'))
    pl.close()
    
# def specrect_plot(outpath, base_name, order_num, before, after):
# 
#     pl.figure('before, after spectral rectify', facecolor='white', figsize=(8, 10))
#     pl.cla()
#     pl.suptitle('before, after spectral rectify, {}, order {}'.format(
#             base_name, order_num), fontsize=14)
#     pl.set_cmap('Blues_r')
# 
#     before_plot = pl.subplot(2, 1, 1)
#     
#     before = imresize(before, (500, 1024), interp='bilinear')
#     
#     try:
#         before_plot.imshow(exposure.equalize_hist(before))
#     except:
#         before_plot.imshow(before)
#     before_plot.set_title('before')
# #     obj_plot.set_ylim([1023, 0])
#     before_plot.set_xlim([0, 1023])
#     
#     after_plot = pl.subplot(2, 1, 2)
#     
#     after = imresize(after, (500, 1024), interp='bilinear')
#     
#     try:
#         after_plot.imshow(exposure.equalize_hist(after))
#     except:
#         after_plot.imshow(after)
#     after_plot.set_title('after')
# #     flat_plot.set_ylim([1023, 0])
#     after_plot.set_xlim([0, 1023])
#  
#     pl.tight_layout()
#     pl.savefig(constructFileName(outpath, base_name, order_num, 'specrect.png'))
#     pl.close()
    
def specrect_plot(outpath, base_name, order_num, before, after):
    
    pl.figure('spectral rectify', facecolor='white')
    pl.cla()
    pl.title('spectral rectify, ' + base_name + ', order ' + str(order_num), fontsize=14)

    pl.xlabel('column (pixels)')
    pl.ylabel('intensity (counts)')
    
    pl.minorticks_on()
    pl.grid(True)
    
    pl.xlim(0, 1023)
    
    pl.plot(before[10, :], "k-", mfc='none', ms=3.0, linewidth=1, 
            label='before')

    pl.plot(after[10, :], "b-", mfc='none', ms=3.0, linewidth=1, 
            label='after')
        
    pl.legend(loc='best', prop={'size': 8})
    
    fn = constructFileName(outpath, base_name, order_num, 'specrect.png')
    pl.savefig(fn)
    pl.close()

    return
    
def perOrderWavelengthCalAsciiTable(outpath, base_name, order, col, centroid, source, wave_exp, wave_fit, res, peak, slope):
    
    names = ['order', 'source', 'col', 'centroid', 'wave_exp', 'wave_fit', 'res', 'peak', 'disp'] 
    units = ['', '', 'pixels', 'pixels', 'Angstroms', 'Angstroms',  'Angstroms', 'counts', 'Angstroms/pixel']
    formats = [ 'd', '', 'd', '.3f', '.6e', '.6e', '.3f', 'd', '.3e']
    nominal_width = 10
    widths = []
    
    for name in names:
        widths.append(max(len(name), nominal_width))
            
    for i in range(len(names)):
        widths[i] = max(len(units[i]), widths[i])
       
    widths[0] = 6 
    widths[1] = 6
    widths[2] = 6
    widths[3] = 8
    widths[4] = 13
    widths[5] = 13
    widths[6] = 9
    widths[7] = 6

    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
      
    if UNDERLINE:  
        line = []
        for i, unit in enumerate(units):
            line.append('{:->{w}}'.format('', w=widths[i]))
        buff.append('| {} |'.format('--'.join(line)))
    
    for i in range(len(order)):
        data = [order[i], source[i], col[i], centroid[i], wave_exp[i], wave_fit[i], res[i], (int)(peak[i]), slope[i]]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
        
    #print('\n'.join(buff))
        
    fn = constructFileName(outpath, base_name, None, 'localcalids.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
    return

    
def skyLinesPlot(outpath, order):
    """
    Always uses frame A.
    """
    
    pl.figure('sky lines', facecolor='white', figsize=(8, 6))
    pl.cla()
    pl.suptitle("sky lines" + ', ' + order.baseNames['A'] + ", order " + 
            str(order.flatOrder.orderNum), fontsize=14)
#     pl.rcParams['ytick.labelsize'] = 8

    syn_plot = pl.subplot(2, 1, 1)
    syn_plot.set_title('synthesized sky')
    syn_plot.set_xlim([0, 1024])
    ymin = np.amin(order.synthesizedSkySpec) - ((np.amax(order.synthesizedSkySpec) - np.amin(order.synthesizedSkySpec)) * 0.1)
    ymax = np.amax(order.synthesizedSkySpec) + ((np.amax(order.synthesizedSkySpec) - np.amin(order.synthesizedSkySpec)) * 0.1)
    syn_plot.set_ylim(ymin, ymax)
    syn_plot.plot(order.synthesizedSkySpec, 'g-', linewidth=1)
    syn_plot.annotate('shift = ' + "{:.3f}".format(order.waveShift), 
                 (0.3, 0.8), xycoords="figure fraction")
    
    sky_plot = pl.subplot(2, 1, 2)
    sky_plot.set_title('sky')
    sky_plot.set_xlim([0, 1024])
    ymin = np.amin(order.skySpec['A']) - ((np.amax(order.skySpec['A']) - np.amin(order.skySpec['A'])) * 0.1)
    ymax = np.amax(order.skySpec['A']) + ((np.amax(order.skySpec['A']) - np.amin(order.skySpec['A'])) * 0.1)
    sky_plot.set_ylim([ymin, ymax])
    sky_plot.plot(order.skySpec['A'], 'b-', linewidth=1)
    
    ymin, ymax = sky_plot.get_ylim()
    dy = (ymax - ymin) / 4
    y = ymin + dy/8
    for line in order.lines:
        if line.frameFitOutlier == False:
            c = 'k--'
        else:
            c = 'r--'
        sky_plot.plot([line.col, line.col], [ymin, ymax], c, linewidth=0.5)
        pl.annotate(str(line.waveAccepted), (line.col, y), size=8)
        pl.annotate(str(line.col) + ', ' + '{:.3f}'.format(
                order.flatOrder.gratingEqWaveScale[line.col]), (line.col, y + (dy / 4)), size=8)
        y += dy
        if y > (ymax - dy):
            y = ymin + dy/8

    fn = constructFileName(outpath, order.baseNames['A'], order.flatOrder.orderNum, 'skylines.png')
        
    pl.savefig(fn)
    pl.close()
    return

def skyLinesAsciiTable(outpath, base_name, order):
    
    names = ['order', 'col', 'calc', 'accepted'] 
    units = [' ', 'pixels', '$\AA$', '$\AA$']
    formats = ['.0f', '.0f', '.1f', '.1f']
    widths = [6, 6, 6, 9] 
    
    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
        
    if UNDERLINE:
        line = []
        for i, unit in enumerate(units):
            line.append('{:->{w}}'.format('', w=widths[i]))
        buff.append('--'.join(line))
    
    for l in order.lines:
        data = [order.flatOrder.orderNum, l.col, order.flatOrder.gratingEqWaveScale[l.col], 
                l.waveAccepted]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
                
    fn = constructFileName(outpath, base_name, order.flatOrder.orderNum, 'skylines.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
 
    return

def wavelengthScalePlot(out_dir, base_name, order):

    pl.figure('wavelength scale', facecolor='white')
    pl.cla()
    pl.title('wavelength scale, ' + base_name + ', order ' + 
             str(order.flatOrder.orderNum), fontsize=14)

    pl.xlabel('column (pixels)')
    pl.ylabel('wavelength ($\AA$)')
    
    pl.minorticks_on()
    pl.grid(True)
    
    pl.xlim(0, 1023)
    
    pl.plot(order.flatOrder.gratingEqWaveScale, "k-", mfc='none', ms=3.0, linewidth=1, 
            label='grating equation')
    if order.waveScale is not None:
        pl.plot(order.waveScale, "b-", mfc='none', ms=3.0, linewidth=1, 
            label='sky lines')
        
    pl.legend(loc='best', prop={'size': 8})
    
    fn = constructFileName(out_dir, base_name, order.flatOrder.orderNum, 'wavelength_scale.png')
    pl.savefig(fn)
    pl.close()

    return