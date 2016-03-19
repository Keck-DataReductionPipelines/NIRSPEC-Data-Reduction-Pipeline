import argparse
import os
import sys
import traceback
import logging
from astropy.io import fits
import warnings

import config
import dgn
import create_raw_data_sets
import reduce_frame
import products
import nirspec_constants as constants
import RawDataSet
from DrpException import DrpException

VERSION = '0.9.9'

warnings.filterwarnings('ignore', category=UserWarning, append=True)

def nsdrp_koa(in_dir, base_out_dir):
    """
    NSDRP. Assembles raw data sets from FITS files in the input directory,
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
            fn = rawDataSet.objFileName.rstrip('.gz').rstrip('.fits')
            if config.params['shortsubdir']:
                out_dir = base_out_dir + '/' + fn[fn[:fn.rfind('.')].rfind('.')+1:]
            else:
                out_dir = base_out_dir + '/' + fn[fn.rfind('/'):]
        else:
            out_dir = base_out_dir        
        if not os.path.exists(out_dir):
                try: 
                    os.mkdir(out_dir)
                except: 
                    msg = 'output directory {} does not exist and cannot be created'.format(out_dir)
                    # logger.critical(msg) can't create log if no output directory
                    raise IOError(msg)
                
        try:
            reduce_data_set(rawDataSet, out_dir)    
        except DrpException as e:
            n_reduced -= 1
            logger.error('failed to reduce {}: {}'.format(
                    rawDataSet.objFileName, e.message))
        except IOError as e:
            logger.critical('DRP failed due to I/O error: {}'.format(str(e)))
            sys.exit(1)
            
    if len(rawDataSets) > 0:
        logger.info('n object frames reduced = {}'.format(
                n_reduced, len(rawDataSets)))   
        
    logger.info('end nsdrp')
    return    
        
def nsdrp_cmnd(fn1, fn2, out_dir):
    
    flat_fn = None
    obj_fn = None
    
    fn1_header = fits.PrimaryHDU.readfrom(fn1, ignore_missing_end=True).header
    fn2_header = fits.PrimaryHDU.readfrom(fn2, ignore_missing_end=True).header

    
    if fn1_header['flat'] == 1 and fn1_header['calmpos'] == 1:
        flat_fn = fn1
        obj_fn = fn2
        obj_header = fn2_header
        flat_header = fn1_header
    if fn2_header['flat'] == 1 and fn2_header['calmpos'] == 1:
        if flat_fn is not None:
            raise DrpException('two flats, no object frame')
        else:
            flat_fn = fn2
            obj_fn = fn1
            obj_header = fn1_header
            flat_header = fn2_header
    if flat_fn is None:
        raise DrpException('no flat')
    
    try:
        if obj_header['DISPERS'].lower() != 'high':
            raise DrpException('DISPERS != high')
    except KeyError:
        print('WARNING: DISPERS header card missing, continuing')
        
    if obj_header['NAXIS1'] != constants.N_COLS:
        raise DrpException('NAXIS1 != {}'.format(constants.N_COLS))
    if obj_header['NAXIS2'] != constants.N_ROWS:
        raise DrpException('NAXIS2 != {}'.format(constants.N_COLS))
    if obj_header['FILNAME'].lower().find('nirspec') < 0:
        raise DrpException('unsupported filter: {}'.format(obj_header['FILNAME']))
    
    if create_raw_data_sets.flat_criteria_met(obj_header, flat_header, ignore_dispers=True) is False:
        raise DrpException('flat is not compatible with object frame')
    
    rawDataSet = RawDataSet.RawDataSet(obj_fn, obj_header)
    rawDataSet.flatFileNames.append(flat_fn)
     
    if not os.path.exists(out_dir):
        try: 
            os.mkdir(out_dir)
        except: 
            msg = 'output directory {} does not exist and cannot be created'.format(out_dir)
            raise IOError(msg)
                
    reduce_data_set(rawDataSet, out_dir)
        
    return        
    
def reduce_data_set(rawDataSet, out_dir):    
    
    logger = logging.getLogger('main')

    # generate reduced data set by reducing raw data set
    reducedDataSet = reduce_frame.reduce_frame(rawDataSet, out_dir)
    
    # produce data products from reduced data set
    if config.params['no_products'] is True:
        logger.info('data product generation inhibited by command line switch')
    else:
        products.gen(reducedDataSet, out_dir)

    # if diagnostics mode is enabled, then produce diagnostic data products
    if config.params['dgn'] is True:
        logger.info('diagnostic mode enabled, generating diagnostic data products')
        dgn.gen(reducedDataSet, out_dir)
                         
        
def init(out_dir, in_dir = None):
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
    if config.params['subdirs'] is False and config.params['cmnd_line_mode'] is False:
        log_dir = out_dir + '/log'
        if not os.path.exists(log_dir):
            try:
                os.makedirs(log_dir)
            except:
                raise IOError('log directory {} does not exist and cannot be created'.format(log_dir))
        
    # set up main logger
    logger = logging.getLogger('main')
    if (config.params['cmnd_line_mode'] is False):
        setup_main_logger(logger, in_dir, out_dir)
    
        logger.info('start nsdrp version {}'.format(VERSION))
        logger.info('cwd: {}'.format(os.getcwd()))
        logger.info('input dir: {}'.format(in_dir.rstrip('/')))
        logger.info('output dir: {}'.format(out_dir.rstrip('/')))
        
        # confirm that input directory exists
        if not os.path.exists(in_dir):
            msg = 'input directory {} does not exist'.format(in_dir)
            logger.critical(msg)
            raise IOError(msg)

    return

def setup_main_logger(logger, in_dir, out_dir):
    if config.params['debug']:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    if config.params['debug']:
        formatter = logging.Formatter('%(asctime)s %(levelname)s - %(filename)s:%(lineno)s - %(message)s')
    else:
        formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
     
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
        
    return

    
def get_log_fn(in_dir, out_dir):
    log_fn = None
    
    if config.params['ut'] is not None:
        log_fn = '{}/NS.{}.log'.format(out_dir, config.params['ut'])
    else:
        fns = os.listdir(in_dir)
        for fn in fns:
            if fn.startswith('NS.'):
                log_fn = out_dir + '/' + fn[:fn.find('.', fn.find('.') + 1)] + '.log'
                break
        if log_fn is None:
            log_fn = out_dir + '/nsdrp.log'
            
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
     
    # parse command line arguments
    parser = argparse.ArgumentParser(description="NSDRP")
    parser.add_argument('arg1', help='input directory (KOA mode) | flat file name (cmnd line mode)')
    parser.add_argument('arg2', help='output directory (KOA mode) | object file name (cmnd line mode)')
    parser.add_argument('-debug', 
            help='enables additional logging for debugging', 
            action='store_true')
    parser.add_argument('-verbose', 
            help='enables output of all log messages to stdout, always true in command line mode',
            action='store_true')
    parser.add_argument('-subdirs',
            help='enables creation of per object frame subdirectories for data products,' +
            'ignored in command line mode',
            action='store_true')
    parser.add_argument('-dgn', 
            help='enables storage of diagnostic data products',
            action='store_true')
    parser.add_argument('-npy', 
            help='enables generation of numpy text files for certain diagnostic data products',
            action='store_true')
    parser.add_argument('-no_cosmic', help='inhibits cosmic ray artifact rejection', 
            action='store_true')
    parser.add_argument('-no_products', help='inhibits data product generation', 
            action='store_true')
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
    parser.add_argument('-pipes', 
            help='enables pipe character seperators in ASCII table headers',
            action='store_true')
    parser.add_argument('-shortsubdir',
            help='use file ID only, rather than full KOA ID, for subdirectory names, ' +
            'ignored in command line mode',
            action='store_true')
    parser.add_argument('-ut',
            help='specify UT to be used for summary log file, overrides auto based on UT in first frame')
    parser.add_argument('-gunzip',
            help='gunzip .gz files (not necessary for processing)', 
            action='store_true')
    parser.add_argument('-out_dir', 
            help='output directory in command line mode [.], ignored in KOA mode')
    args = parser.parse_args()
    config.params['debug'] = args.debug
    config.params['verbose'] = args.verbose
    config.params['subdirs'] = args.subdirs
    config.params['dgn'] = args.dgn
    config.params['npy'] = args.npy
    config.params['no_cosmic'] = args.no_cosmic
    config.params['no_products'] = args.no_products
    if args.obj_window is not None:
        config.params['obj_window'] = int(args.obj_window)
    if args.sky_window is not None:
        config.params['sky_window'] = int(args.sky_window)
    if args.sky_separation is not None:
        config.params['sky_separation'] = int(args.sky_separation)
    if args.oh_filename is not None:
        config.params['oh_filename'] = args.oh_filename
        config.params['oh_envar_override'] = True
    config.params['int_c'] = args.int_c
    config.params['lla'] = args.lla
    config.params['pipes'] = args.pipes
    config.params['shortsubdir'] = args.shortsubdir
    if args.ut is not None:
        config.params['ut'] = args.ut
    config.params['gunzip'] = args.gunzip
    if args.out_dir is not None:
        config.params['out_dir'] = args.out_dir

    
    # determine if we are in command line mode or KOA mode
    try:
#         fits.getheader(args.arg1)
#         fits.getheader(args.arg2)
        fits.PrimaryHDU.readfrom(args.arg1, ignore_missing_end=True)
        fits.PrimaryHDU.readfrom(args.arg2, ignore_missing_end=True)
    except IOError:
        # these aren't FITS files so must be in KOA mode
        print('KOA mode')
    else:
        # command line mode
        config.params['cmnd_line_mode'] = True
        config.params['verbose'] = True

    

    # initialize environment, setup main logger, check directories
    try:
        if config.params['cmnd_line_mode'] is True:
            init(config.params['out_dir'])
            nsdrp_cmnd(args.arg1, args.arg2, config.params['out_dir'])
        else:
            init(args.arg2, args.arg1)
            nsdrp_koa(args.arg1, args.arg2)
    except Exception as e:
        print('ERROR: ' + e.message)
        if config.params['debug'] is True:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
            traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)
        sys.exit(2)    

    sys.exit(0)
    
    
if __name__ == "__main__":
    """
    NIRSPEC DRP v2.0
    """
    main()   
    
    
    
    
