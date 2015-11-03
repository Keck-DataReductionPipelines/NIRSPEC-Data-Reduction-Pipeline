import pylab as pl
import logging
import os
import errno
import numpy as np
from skimage import exposure

import Order
import ReducedDataSet

logger = logging.getLogger('obj')

subdirs = dict([
                ('traces.png',          'traces'),
                ('trace.npy',           'traces'),
                ('cutouts.png',         'cutouts'),
                ('obj_cutout.npy',      'cutouts'),
                ('flat_cutout.npy',     'cutouts'),
                ('sprect.png',          'sprect'),
                ('edges.png',           'edges'),
                ('top_bot_edges.png',   'edges')
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

    edges_plot(out_dir, reduced.baseName, reduced.topEdgesProfile, reduced.botEdgesProfile,
            reduced.topEdgePeaks, reduced.botEdgePeaks)
    
    tops_bots_plot(out_dir, reduced.baseName, reduced.topEdgesImg, reduced.botEdgesImg)
    
    for order in reduced.orders:
        
        traces_plot(out_dir, reduced.baseName, order.orderNum, reduced.obj, reduced.flat, 
                order.topTrace, order.botTrace)
        
        # save smooted trace cutout to numpy text file
        fn = constructFileName(out_dir, reduced.baseName, order.orderNum, 'trace.npy')
        np.savetxt(fn, order.smoothedTrace)
        
        # save obj and flat cutouts to numpy text files
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
            
        sprect_plot(out_dir, reduced.baseName, order.orderNum, 
                order.srFlatObjImg, order.srNormFlatImg)
        
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
    
def sprect_plot(outpath, base_name, order_num, obj, flat):

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
    pl.savefig(constructFileName(outpath, base_name, order_num, 'sprect.png'))
    pl.close()