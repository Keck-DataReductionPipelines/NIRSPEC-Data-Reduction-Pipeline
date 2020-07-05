import os
import check_modules

requiredModules = check_modules.is_missing()
if requiredModules:
    print(('Missing modules: ' + ', '.join(requiredModules)))
    os.sys.exit()

import config
import nsdrp_cmnd
import nsdrp_koa
import argparse
import logging
from astropy.io import fits
import warnings

# from DrpException import DrpException
# import FlatCacher

VERSION = '0.9.17'

warnings.filterwarnings('ignore', category=UserWarning, append=True)


def main():
    """
    Main entry point for DRP.

    Run with -h to see command line arguments
    """

    args = parse_cmnd_line_args()

    # determine if we are in command line mode or KOA mode
    try:
        fits.PrimaryHDU.readfrom(args.arg1, ignore_missing_end=True)
        fits.PrimaryHDU.readfrom(args.arg2, ignore_missing_end=True)
    except IOError:
        # these aren't FITS files so must be in KOA mode
        print('KOA mode')
    else:
        # command line mode
        config.params['cmnd_line_mode'] = True
        config.params['verbose'] = True

    # setup configuration parameters based on command line args
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
    config.params['spatial_jump_override'] = args.spatial_jump_override
    if args.out_dir is not None:
        config.params['out_dir'] = args.out_dir
    config.params['jpg'] = args.jpg
    config.params['sowc'] = args.sowc

    # initialize environment, setup main logger, check directories
#     try:
    if config.params['cmnd_line_mode'] is True:
        init(config.params['out_dir'])
        nsdrp_cmnd.process_frame(
            args.arg1,
            args.arg2,
            args.b,
            config.params['out_dir'])
    else:
        init(args.arg2, args.arg1)
        nsdrp_koa.process_dir(args.arg1, args.arg2)
#     except Exception as e:
#         print('ERROR: ' + e.message)
#         if config.params['debug'] is True:
#             exc_type, exc_value, exc_traceback = sys.exc_info()
#             traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
#             traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)
#         sys.exit(2)

    os.sys.exit(0)


def init(out_dir, in_dir=None):
    """
    Sets up main logger, checks for existence of input directory, and checks that
    output directory either exists or can be created.

    """

    # create output directory and log subdirectory if -subdirs not set
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except BaseException:
            msg = 'output directory {} does not exist and cannot be created'.format(
                out_dir)
            # logger.critical(msg) can't create log if no output directory
            raise IOError(msg)
    if config.params['subdirs'] is False and config.params['cmnd_line_mode'] is False:
        log_dir = out_dir + '/log'
        if not os.path.exists(log_dir):
            try:
                os.makedirs(log_dir)
            except BaseException:
                raise IOError(
                    'log directory {} does not exist and cannot be created'.format(log_dir))

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
        formatter = logging.Formatter(
            '%(asctime)s %(levelname)s - %(filename)s:%(lineno)s - %(message)s')
    else:
        formatter = logging.Formatter(
            '%(asctime)s %(levelname)s - %(message)s')

    log_fn = get_log_fn(in_dir, out_dir)

    if os.path.exists(log_fn):
        os.rename(log_fn, log_fn + '.prev')

    fh = logging.FileHandler(filename=log_fn)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    if config.params['debug']:
        sformatter = logging.Formatter(
            '%(asctime)s %(levelname)s - %(filename)s:%(lineno)s - %(message)s')
    else:
        sformatter = logging.Formatter(
            '%(asctime)s %(levelname)s - %(message)s')

    if config.params['verbose'] is True:
        sh = logging.StreamHandler()
        sh.setLevel(logging.DEBUG)
        sh.setFormatter(sformatter)
        logger.addHandler(sh)

    return


def get_log_fn(in_dir, out_dir):
    """
    """
    log_fn = None

    if config.params['ut'] is not None:
        # use UT specified as command line argument
        log_fn = '{}/NS.{}.log'.format(out_dir, config.params['ut'])
    else:
        # try to get UT from filenames in input directory
        fns = os.listdir(in_dir)
        for fn in fns:
            if fn.startswith('NS.'):
                log_fn = out_dir + '/' + \
                    fn[:fn.find('.', fn.find('.') + 1)] + '.log'
                break
        if log_fn is None:
            # if all else fails, use canned log file name
            log_fn = out_dir + '/nsdrp.log'

    if config.params['subdirs'] is False:
        # if not in "subdirs" mode than put log file in log subdirectory
        parts = log_fn.split('/')
        parts.insert(len(parts) - 1, 'log')
        log_fn = '/'.join(parts)

    return(log_fn)


def parse_cmnd_line_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description="NSDRP")
    parser.add_argument(
        'arg1',
        help='input directory (KOA mode) | flat file name (cmnd line mode)')
    parser.add_argument(
        'arg2',
        help='output directory (KOA mode) | object file name (cmnd line mode)')
    parser.add_argument('-debug',
                        help='enables additional logging for debugging',
                        action='store_true')
    parser.add_argument('-verbose',
                        help='enables output of all log messages to stdout, always true in command line mode',
                        action='store_true')
    parser.add_argument('-subdirs',
                        help='enables creation of per-object-frame subdirectories for data products,' +
                        'ignored in command line mode',
                        action='store_true')
    parser.add_argument('-dgn',
                        help='enables saving diagnostic data products in ...out/diagnostics',
                        action='store_true')
    parser.add_argument('-npy',
                        help='enables generation of numpy text files for certain diagnostic data products',
                        action='store_true')
    parser.add_argument('-no_cosmic', help='inhibits cosmic ray artifact rejection',
                        action='store_true')
    parser.add_argument('-no_products', help='inhibits data product generation',
                        action='store_true')
#     , default=config.DEFAULT_COSMIC)
    parser.add_argument(
        '-obj_window',
        help='object extraction window width in pixels')
    # default=config.DEFAULT_OBJ_WINDOW)
    parser.add_argument(
        '-sky_window',
        help='background extraction window width in pixels')
    # default=config.DEFAULT_SKY_WINDOW)
    parser.add_argument(
        '-sky_separation',
        help='separation between object and sky windows in pixels')
    # default=config.DEFAULT_SKY_DIST)
    parser.add_argument(
        '-oh_filename',
        help='path and filename of OH emission line catalog file')
    parser.add_argument('-int_c', help='user integer column values rather than fractional values \
            determined by centroiding in wavelength fit',
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
                        help='specify UT to be used for summary log file, overrides automatic UT \
            determination based on UT of first frame')
    parser.add_argument('-gunzip',
                        help='forces decompression of compressed FITS files, leaves both the .gz and .fits \
                files in the source directory.  Note that the compressed files can be read \
                directly so it is not necessary to decompress them.', action='store_true')
    parser.add_argument('-spatial_jump_override',
                        help='inhibit rejection of order edge traces based on \'jump\' limit)',
                        action='store_true')
    parser.add_argument('-out_dir',
                        help='output directory path used in command line mode, default is current working \
            directory, ignored in KOA mode')
    parser.add_argument('-b',
                        help='filename of frame B in AB pair')
    parser.add_argument('-jpg', help='store preview plots in JPG format instead of PNG',
                        action='store_true')
    parser.add_argument(
        '-sowc',
        help='enable simple order width calculation',
        action='store_true')

    return(parser.parse_args())


if __name__ == "__main__":
    """
    NIRSPEC DRP
    """
    main()
