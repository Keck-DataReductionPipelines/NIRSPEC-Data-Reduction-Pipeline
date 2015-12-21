import pylab as pl
import logging
import os
import errno
import numpy as np
from skimage import exposure

import Order
import ReducedDataSet
import config

UNDERLINE = False

logger = logging.getLogger('obj')

subdirs = dict([
                ('traces.png',          'traces'),
                ('trace.npy',           'traces'),
                ('cutouts.png',         'cutouts'),
                ('obj_cutout.npy',      'cutouts'),
                ('flat_cutout.npy',     'cutouts'),
                ('sparect.png',         'sparect'),
                ('edges.png',           'edges'),
                ('top_bot_edges.png',   'edges'),
                ('localcalids.txt',     'localcal'),
                ('skylines.png',        'skylines'),
                ('skylines.txt',        'skylines')
                ])

def gen(reduced, out_dir):
    
    logger.info('generating diagnostic data products...')
    
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
        if order.perOrderSlope > 0.95 and order.perOrderSlope < 1.05:
            for line in order.lines:
                order_num.append(order.orderNum)
                col.append(line.col)
                centroid.append(line.centroid)
                source.append('sky')
                wave_exp.append(line.acceptedWavelength)
                wave_fit.append(line.localFitWavelength)
                res.append(line.localFitResidual)
                peak.append(line.peak)
                slope.append(line.localFitSlope)
                
    # per-order wavelength table
    perOrderWavelengthCalAsciiTable(
            out_dir, reduced.baseName, order_num, col, centroid, source, wave_exp, wave_fit, res, peak, slope)
    
    edges_plot(out_dir, reduced.baseName, reduced.topEdgesProfile, reduced.botEdgesProfile,
            reduced.topEdgePeaks, reduced.botEdgePeaks)
    
    tops_bots_plot(out_dir, reduced.baseName, reduced.topEdgesImg, reduced.botEdgesImg)
    
    for order in reduced.orders:
        
        traces_plot(out_dir, reduced.baseName, order.orderNum, reduced.obj, reduced.flat, 
                order.topTrace, order.botTrace)
        
        # save smooted trace cutout to numpy text file
        if config.params['npy'] is True:
            fn = constructFileName(out_dir, reduced.baseName, order.orderNum, 'trace.npy')
            np.savetxt(fn, order.smoothedTrace)
        
        # save obj and flat cutouts to numpy text files
        if config.params['npy'] is True:
            fn = constructFileName(out_dir, reduced.baseName, order.orderNum, 'obj_cutout.npy')
            np.savetxt(fn, order.objCutout)
            fn = constructFileName(out_dir, reduced.baseName, order.orderNum, 'flat_cutout.npy')
            np.savetxt(fn, order.flatCutout)
        
        if order.lowestPoint > order.padding:
            cutouts_plot(out_dir, reduced.baseName, order.orderNum, order.objCutout, order.flatCutout, 
                    order.topTrace - order.lowestPoint + order.padding, 
                    order.botTrace - order.lowestPoint + order.padding,
                    order.smoothedTrace - order.lowestPoint + order.padding)
        else:
            cutouts_plot(out_dir, reduced.baseName, order.orderNum, order.objCutout, order.flatCutout, 
                    order.topTrace, order.botTrace, order.smoothedTrace)
            
        sparect_plot(out_dir, reduced.baseName, order.orderNum, 
                order.srFlatObjImg, order.srNormFlatImg)
        
        skyLinesPlot(out_dir, reduced.baseName, order)
        skyLinesAsciiTable(out_dir, reduced.baseName, order)
        
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
#         pl.plot([peak[0], peak[0]], [(tops_plot.get_ylim())[0], (tops_plot.get_ylim())[1]], "g-", linewidth=1.0, label='peak')
        pl.annotate(str(peak[0]), (peak[0], peak[1]), size=8)
    tops_plot.set_xlim([0, 1023])

    bots_plot = pl.subplot(2, 1, 2)
    bots_plot.set_title('bottoms')
    bots_plot.plot(bot_profile, 'k-', linewidth=1.0)
    for peak in bot_peaks:
#         pl.plot([peak[0], peak[0]], [(bots_plot.get_ylim())[0], (bots_plot.get_ylim())[1]], "g-", linewidth=1.0, label='peak')    
        pl.annotate(str(peak[0]), (peak[0], peak[1]), size=8)
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
    obj_plot.imshow(exposure.equalize_hist(tops))
    obj_plot.set_title('top edges')
    obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])

    
    flat_plot = pl.subplot(1, 2, 2)
    flat_plot.imshow(exposure.equalize_hist(bots))
    flat_plot.set_title('bottom edges')
    flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, None, 'top_bot_edges.png'))
    pl.close()
    
def traces_plot(outpath, base_name, order_num, obj, flat, top_trace, bot_trace):
    
    pl.figure('traces', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('order edge traces, {}, order {}'.format(base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')
    pl.rcParams['ytick.labelsize'] = 8

    obj_plot = pl.subplot(1, 2, 1)
    obj_plot.imshow(exposure.equalize_hist(obj))
    obj_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)

    obj_plot.set_title('object')
    obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])

    
    flat_plot = pl.subplot(1, 2, 2)
    flat_plot.imshow(exposure.equalize_hist(flat))
    flat_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    flat_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)    
    flat_plot.set_title('flat')
    flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, order_num, 'traces.png'))
    pl.close()
    
def cutouts_plot(outpath, base_name, order_num, obj, flat, top_trace, bot_trace, trace):
    
    pl.figure('traces', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('order cutouts, {}, order {}'.format(base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')

    obj_plot = pl.subplot(2, 1, 1)
    obj_plot.imshow(exposure.equalize_hist(obj))
    obj_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), trace, 'y-', linewidth=1.5)
    obj_plot.set_title('object')
#     obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])
    
    flat_plot = pl.subplot(2, 1, 2)
    flat_plot.imshow(exposure.equalize_hist(flat))
    flat_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    flat_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)    
    flat_plot.plot(np.arange(1024), trace, 'y-', linewidth=1.5)    
    flat_plot.set_title('flat')
#     flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, order_num, 'cutouts.png'))
    pl.close()
    
def sparect_plot(outpath, base_name, order_num, obj, flat):

    pl.figure('spatially rectified', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('spatially rectified, {}, order {}'.format(base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')

    obj_plot = pl.subplot(2, 1, 1)
    obj_plot.imshow(exposure.equalize_hist(obj))
    obj_plot.set_title('object')
#     obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])
    
    flat_plot = pl.subplot(2, 1, 2)
    flat_plot.imshow(exposure.equalize_hist(flat))
    flat_plot.set_title('flat')
#     flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, order_num, 'sparect.png'))
    pl.close()
    
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

def skyLinesPlot(outpath, base_name, order):
    pl.figure('sky lines', facecolor='white', figsize=(8,5))
    pl.cla()
    pl.title("sky lines" + ', ' + base_name + ", order " + str(order.orderNum), fontsize=14)
    pl.xlabel('x (pixels)', fontsize=12)
    #pl.ylabel('row (pixel)', fontsize=12)
    
    synmax = np.amax(order.synthesizedSkySpec);
    skymax = np.amax(order.skySpec);
    

    if synmax > skymax:
        syn = order.synthesizedSkySpec
#         sky = (order.skySpec - np.amin(order.skySpec)) * (synmax / skymax)
        sky = order.skySpec * (synmax / skymax)
    else:
        syn = order.synthesizedSkySpec * (skymax / synmax)
#         sky = (order.skySpec - np.amin(order.skySpec))
        sky = order.skySpec

    pl.plot((syn * 0.4) + (max(synmax, skymax) / 2), 'g-', linewidth=1, label='synthesized sky')
    pl.plot(sky * 0.4, 'b-', linewidth=1, label='sky')

    ymin, ymax = pl.ylim()
    ymax = max(synmax, skymax) * 1.1
    ymin = -200
    pl.ylim(ymin, ymax)
    pl.xlim(0, 1024)

    pl.legend(loc='best', prop={'size': 8})
    
    pl.annotate('shift = ' + "{:.3f}".format(order.wavelengthShift), 
                (0.3, 0.8), xycoords="figure fraction")
    
    dy = 0
    for line in order.lines:
        pl.plot([line.col, line.col], [ymin, ymax], "k--", linewidth=0.5)
        pl.annotate(str(line.acceptedWavelength), (line.col, ((ymax - ymin) / 2) + dy), size=8)
        pl.annotate(str(line.col) + ', ' + '{:.3f}'.format(order.wavelengthScaleCalc[line.col]), (line.col, ((ymax - ymin) / 3) - dy), size=8)
        dy += 300
        if dy > 1500:
            dy = 0

    
    #pl.minorticks_on()
    #pl.grid(True)

    fn = constructFileName(outpath, base_name, order.orderNum, 'skylines.png')
        
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
        data = [order.orderNum, l.col, order.wavelengthScaleCalc[l.col], l.acceptedWavelength]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
                
    fn = constructFileName(outpath, base_name, order.orderNum, 'skylines.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
 
    return