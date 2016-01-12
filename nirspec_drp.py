import argparse
import os
import sys
import logging

import config
import dgn
import create_raw_data_sets
#import ReducedDataSet
import reduce_frame
import products
import DrpException

cosmic_clean = True
config_params = {}

def nirspec_drp(in_dir, base_out_dir):
    """
    NIRSPEC DRP. Assembles raw data sets from FITS files in the input directory,
    then generates reduced data sets from raw data sets.  Level 1 data products
    are generated from the reduced data sets and saved in the output directory.
    """
        
    logger = logging.getLogger('main')
    
    rawDataSets = create_raw_data_sets.create(in_dir)
    n_reduced = len(rawDataSets)
    
    logger.info(str(len(rawDataSets)) + " raw data set(s) assembled")
        
    # process each raw data set
    
    for rawDataSet in rawDataSets:

        if config.params['subdirs'] is True:
            fn = rawDataSet.objFileName
            out_dir = base_out_dir + '/' + fn[fn[:fn.rfind('.')].rfind('.')+1:fn.rfind('.')]
            if not os.path.exists(out_dir):
                try: 
                    os.mkdir(out_dir)
                except: 
                    msg = 'output directory {} does not exist and cannot be created'.format(out_dir)
                    # logger.critical(msg) can't create log if no output directory
                    raise IOError(msg)
        else:
            out_dir = base_out_dir
            
        try:
            # generate reduced data set by reducing raw data set
            reducedDataSet = reduce_frame.reduce_frame(rawDataSet, out_dir)
            
#             logger.info('mean snr = {:.3f}'.format(
#                     sum(reducedDataSet.orders[i].snr for i in range(len(reducedDataSet.orders))) / 
#                     len(reducedDataSet.orders)))
            
            # produce data products from reduced data set
            if config.params['products'] is True:
                products.gen(reducedDataSet, out_dir)
            else:
                logger.info('data product generation inhibited by command line switch')
            
            # if diagnostics mode is enabled, then produce diagnostic data products
            if config.params['dgn'] is True:
                logger.info('diagnostic mode enabled, generating diagnostic data products')
                dgn.gen(reducedDataSet, out_dir)
                
        except DrpException.DrpException as e:
            n_reduced -= 1
            logger.error('failed to reduce {}: {}'.format(
                    rawDataSet.objFileName, e.message))
        except IOError as e:
            logger.critical('DRP failed due to I/O error: {}'.format(str(e)))
            sys.exit(1)
            
    if len(rawDataSets) > 0:
        logger.info('n object frames reduced = {}'.format(
                n_reduced, len(rawDataSets)))   
        
    logger.info('end nirspec drp')
    return    
        
        
def init(in_dir, out_dir):
    """
    Sets up main logger, checks for existence of input directory, and checks that
    output directory either exists or can be created.
    
    """
    
    # create output directory and log subdirectory if -subdirs not set
    if not os.path.exists(out_dir):
        try: 
            os.makedirs(out_dir)
        except: 
            msg = 'output directory {} does not exist and cannot be created'.format(out_dir)
            # logger.critical(msg) can't create log if no output directory
            raise IOError(msg)
    if config.params['subdirs'] is False:
        log_dir = out_dir + '/log'
        if not os.path.exists(log_dir):
            try:
                os.makedirs(log_dir)
            except:
                raise IOError('log directory {} does not exist and cannot be created'.format(log_dir))
        
    # set up main logger
    logger = logging.getLogger('main')
    if config.params['debug']:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    if config.params['debug']:
        formatter = logging.Formatter('%(asctime)s %(levelname)s - %(filename)s:%(lineno)s - %(message)s')
    else:
        formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
    
#     log_fn = None
#     fns = os.listdir(in_dir)
#     for fn in fns:
#         if fn.startswith('NS.'):
#             log_fn = out_dir + '/' + fn[:fn.find('.', fn.find('.') + 1)] + '.log'
#             break
#     if log_fn is None:
#         log_fn = out_dir + '/nirspec_drp.log'
#     if config.params['subdirs'] is False:
#         parts = log_fn.split('/')
#         parts.insert(len(parts)-1, 'log')
#         log_fn = '/'.join(parts)

    log_fn = get_log_fn(in_dir, out_dir)
            
    if os.path.exists(log_fn):
        os.rename(log_fn, log_fn + '.prev')
        
    fh = logging.FileHandler(filename=log_fn)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    if config.params['debug']:
        sformatter = logging.Formatter('%(asctime)s %(levelname)s - %(filename)s:%(lineno)s -  %(message)s')
    else:   
        sformatter = logging.Formatter('%(asctime)s %(levelname)s -  %(message)s')

    if config.params['verbose'] is True:
        sh = logging.StreamHandler()
        sh.setLevel(logging.DEBUG)
        sh.setFormatter(sformatter)
        logger.addHandler(sh)
    
    logger.info('start nirspec drp')
    logger.info('cwd: {}'.format(os.getcwd()))
    logger.info('input dir: {}'.format(in_dir.rstrip('/')))
    logger.info('output dir: {}'.format(out_dir.rstrip('/')))
    
    # confirm that input directory exists
    if not os.path.exists(in_dir):
        msg = 'input directory {} does not exist'.format(in_dir)
        logger.critical(msg)
        raise IOError(msg)

    return

def get_log_fn(in_dir, out_dir):
    log_fn = None
    fns = os.listdir(in_dir)
    for fn in fns:
        if fn.startswith('NS.'):
            log_fn = out_dir + '/' + fn[:fn.find('.', fn.find('.') + 1)] + '.log'
            break
    if log_fn is None:
        log_fn = out_dir + '/nirspec_drp.log'
    if config.params['subdirs'] is False:
        parts = log_fn.split('/')
        parts.insert(len(parts)-1, 'log')
        log_fn = '/'.join(parts)
    return(log_fn)
     
def main():
    """
    Main entry point for DRP.  Parses command line arguments, 
    calls init() to initialize environment and then calls drp()
    to process data.
    
    Expects name of input directory containing raw FITS files and
    name of root output directory to be supplied on the command line.
    
    Run with -h to see command line arguments
    """
     
#     global cosmic_clean
#     global config_params
    
    # parse command line arguments
    parser = argparse.ArgumentParser(description="NIRSPEC DRP")
    parser.add_argument('in_dir', help='input directory')
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('-debug', 
            help='enables additional logging for debugging', 
            action='store_true')
    parser.add_argument('-verbose', 
            help='enables output of all log messages to stdout',
            action='store_true')
    parser.add_argument('-subdirs',
            help='enables creation of per object frame subdirectories for data products',
            action='store_true')
    parser.add_argument('-dgn', 
            help='enables storage of diagnostic data products',
            action='store_true')
    parser.add_argument('-npy', 
            help='enables generation of numpy text files for certain diagnostic data products',
            action='store_true')
    parser.add_argument('-cosmic', help='inhibits cosmic ray artifact rejection', 
            action='store_false')
    parser.add_argument('-products', help='inhibits data product generation', 
            action='store_false')
#     , default=config.DEFAULT_COSMIC)
    parser.add_argument('-obj_window', help='object extraction window width in pixels')
    #default=config.DEFAULT_OBJ_WINDOW)
    parser.add_argument('-sky_window', help='background extraction window width in pixels')
    #default=config.DEFAULT_SKY_WINDOW)
    parser.add_argument('-sky_separation', help='separation between object and sky windows in pixels')
    #default=config.DEFAULT_SKY_DIST)
    parser.add_argument('-oh_filename', help='path and filename of OH emission line catalog file')
    parser.add_argument('-int_c', help='revert to using integer column values in wavelength fit',
            action='store_true')
    parser.add_argument('-lla', type=int, default=2, 
            help='calibration line location algorithm, 1 or [2]')
    args = parser.parse_args()
    config.params['debug'] = args.debug
    config.params['verbose'] = args.verbose
    config.params['subdirs'] = args.subdirs
    config.params['dgn'] = args.dgn
    config.params['npy'] = args.npy
    config.params['cosmic'] = args.cosmic
    config.params['products'] = args.products
    if args.obj_window is not None:
        config.params['obj_window'] = int(args.obj_window_width)
    if args.sky_window is not None:
        config.params['sky_window'] = int(args.sky_window_width)
    if args.sky_separation is not None:
        config.params['sky_separation'] = int(args.sky_separation)
    if args.oh_filename is not None:
        config.params['oh_filename'] = args.oh_filename
    config.params['int_c'] = args.int_c
    config.params['lla'] = args.lla

    # initialize environment, setup main logger, check directories
    try:
        init(args.in_dir, args.out_dir)
    except Exception as e:
        print(e)
        sys.exit(2)    
    
    # process data
#     try: 
#         nirspec_drp(args.in_dir, args.out_dir)
#     except Exception as e:
#         print(e)
#         sys.exit(1)
        
    nirspec_drp(args.in_dir, args.out_dir)

    sys.exit(0)
    
    
if __name__ == "__main__":
    """
    NIRSPEC DRP v2.0
    """
    main()   
    
    
    
    