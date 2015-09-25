import argparse
import os
import sys
import logging

#import RawDataSet
import create_raw_data_sets
#import ReducedDataSet
import reduce_frame
import products


def nirspec_drp(in_dir, out_dir):
    """
    NIRSPEC DRP. Assembles raw data sets from FITS files in the input directory,
    then generates reduced data sets from raw data sets.  Level 1 data products
    are generated from the reduced data sets and saved in the output directory.
    """
    
    logger = logging.getLogger('main')
    
    rawDataSets = create_raw_data_sets.create(in_dir)
    n_reduced = len(rawDataSets)
    
    logger.info(str(len(rawDataSets)) + " raw data set(s) assembled")
        
    for rawDataSet in rawDataSets:
        try:
            reducedDataSet = reduce_frame.reduce_frame(rawDataSet, out_dir)
            products.gen(reducedDataSet, out_dir)
        except Exception as e:
            n_reduced -= 1
            logger.error('failed to reduce {}: {}'.format(
                    rawDataSet.objFileName, e.message))
            
    logger.info('{} out of {} object frames successfully reduced'.format(
            n_reduced, len(rawDataSets)))   
    return    
        
        
def init(in_dir, out_dir):
    """
    Sets up main logger, checks for existance of input directory, and checks that
    output directory either exists or can be created.
    
    """
    
    
    if not os.path.exists(out_dir):
        try: 
            os.mkdir(out_dir)
        except: 
            msg = 'output directory {} does not exist and cannot be created'.format(out_dir)
            # logger.critical(msg) can't create log if no output directory
            raise IOError(msg)
        
        
    # set up main logger
    logger = logging.getLogger('main')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - ' +
            '%(levelname)s - %(filename)s:%(lineno)s - %(message)s')
#     formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    fn = out_dir + '/nirspec_drp.log'
    if os.path.exists(fn):
        os.rename(fn, fn + '.prev')
        
    fh = logging.FileHandler(filename=fn)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    sformatter = logging.Formatter('%(levelname)s - %(filename)s:%(lineno)s -  %(message)s')
    sh = logging.StreamHandler()
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(sformatter)
    logger.addHandler(sh)
    
    logger.info('start')
    logger.info('cwd: {}'.format(os.getcwd()))
    logger.info('input dir: {}'.format(in_dir.rstrip('/')))
    logger.info('output dir: {}'.format(out_dir.rstrip('/')))
    
    # confirm that input directory exists
    if not os.path.exists(in_dir):
        msg = 'input directory {} does not exist'.format(in_dir)
        logger.critical(msg)
        raise IOError(msg)
    


    return
     
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
    parser = argparse.ArgumentParser(description="NIRSPEC DRP")
    parser.add_argument('in_dir', help='input directory')
    parser.add_argument('out_dir', help='output directory')
    args = parser.parse_args()
    
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
    
    
    
    