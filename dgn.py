import pylab as pl
import logging
import os
import errno
import numpy as np


import Order
import ReducedDataSet

logger = logging.getLogger('obj')

subdirs = dict([
                ('traces.png',      'traces'),
                ('cutouts.png',     'cutouts')
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
            
    for order in reduced.orders:
        traces_plot(out_dir, reduced.baseName, order.orderNum, reduced.obj, reduced.flat, 
                order.topTrace, order.botTrace)
        cutouts_plot(out_dir, reduced.baseName, order.orderNum, order.objCutout, order.flatCutout, 
                order.topTrace - order.lowestPoint + order.padding, 
                order.botTrace - order.lowestPoint + order.padding)
        
    logger.info('done generating diagnostic data products')
    
def constructFileName(outpath, base_name, order, fn_suffix):
    fn = outpath + '/diagnostics/' + subdirs[fn_suffix] + '/' + base_name + '_' + fn_suffix
    if order is None:
        return fn
    else:
        return fn[:fn.rfind('_') + 1] + str(order) + fn[fn.rfind('_'):]    

def traces_plot(outpath, base_name, order_num, obj, flat, top_trace, bot_trace):
    
    pl.figure('traces', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('order edge traces, {}, order {}'.format(base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')
    pl.rcParams['ytick.labelsize'] = 8

    obj_plot = pl.subplot(1, 2, 1)
    obj_plot.imshow(obj, vmin=0, vmax=256)
    obj_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)
    obj_plot.set_title('object')
    obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])

    
    flat_plot = pl.subplot(1, 2, 2)
    flat_plot.imshow(flat, vmin=0, vmax=256)
    flat_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    flat_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)    
    flat_plot.set_title('flat')
    flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, order_num, 'traces.png'))
    pl.close()
    
def cutouts_plot(outpath, base_name, order_num, obj, flat, top_trace, bot_trace):
    
    pl.figure('traces', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.suptitle('order edge traces, {}, order {}'.format(base_name, order_num), fontsize=14)
    pl.set_cmap('Blues_r')

    obj_plot = pl.subplot(2, 1, 1)
    obj_plot.imshow(obj, vmin=0, vmax=256)
    obj_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    obj_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)
    obj_plot.set_title('object')
#     obj_plot.set_ylim([1023, 0])
    obj_plot.set_xlim([0, 1023])
    
    flat_plot = pl.subplot(2, 1, 2)
    flat_plot.imshow(flat, vmin=0, vmax=256)
    flat_plot.plot(np.arange(1024), top_trace, 'y-', linewidth=1.5)
    flat_plot.plot(np.arange(1024), bot_trace, 'y-', linewidth=1.5)    
    flat_plot.set_title('flat')
#     flat_plot.set_ylim([1023, 0])
    flat_plot.set_xlim([0, 1023])
 
    pl.tight_layout()
    pl.savefig(constructFileName(outpath, base_name, order_num, 'cutouts.png'))
    pl.close()
    


