import os
#import subprocess
import fnmatch
import logging
from astropy.io import fits

import RawDataSet
import nirspec_constants as constants
import config

from __builtin__ import False

failed2reduce = {'itype':0, 'dispmode':0, 'n1':0, 'n2':0, 'fil':0} #, 'dmode':0}

def create(in_dir):
    """
    Given an input directory path, creates and returns a list of RawDataSet objects.
    A raw data set consists of an object file name, it's header, and associated flats 
    and darks.
    
    param 
        in_dir Input directory path.
        
    return 
        list of RawDataSet objects.
    
    """
    
    logger = logging.getLogger('main')
    logger.debug('creating raw data sets from files in ' + in_dir)
    
    # get list of fits files
    headers = get_headers(in_dir)
            
    if (len(headers) == 0):
        logger.critical('no fits files found')
        print "CRITICAL - no fits files found in: " + in_dir
        os.sys.exit(2)
        return
    
    logger.info('n fits files found = {}'.format(str(len(headers))))
    
    rawDataSets = []

    for filename, header in headers.items():
            if obj_criteria_met(header, failed2reduce):
                rawDataSets.append(RawDataSet.RawDataSet(filename, None, header))

                
    if obj_criteria_met(header, failed2reduce) is False:
#        if failed2reduce.get('itype') > 0:
#            logger.info('Ignored {} files because they are not object frames'.format(
#                    failed2reduce.get('itype')))
        if failed2reduce.get('dismode') > 0:
            logger.info('Failed to reduced {} files because of low dispersion mode'.format(
                    failed2reduce.get('dispmode')))
        if failed2reduce.get('n1') > 0:
            logger.info('Failed to reduced {} files because NAXIS1 != 1024'.format(
                    failed2reduce.get('n1')))
        if failed2reduce.get('n2') > 0:
            logger.info('Failed to reduced {} files because NAXIS2 != 1024'.format(
                    failed2reduce.get('n2')))
        if failed2reduce.get('fil') > 0:
            logger.info('Failed to reduce {} files because of unsupported filter'.format(
                    failed2reduce.get('fil')))
#        if failed2reduce.get('dmode') > 0:
#            logger.info('Failed to reduce {} files because SCAM mode'.format(
#                    failed2reduce.get('dmode')))
    
    logger.info('n object frames found = {}'.format(str(len(rawDataSets))))
    
    # associate darks and flats with each object frame
    for rawDataSet in rawDataSets:
        for filename, header in headers.items():
            if (header['IMAGETYP'] == 'flatlamp'):
                if len(rawDataSet.flatFns) < config.params['max_n_flats']:
                    if flat_criteria_met(rawDataSet.objHeader, header):
                        rawDataSet.flatFns.append(filename)
            elif (header['IMAGETYP'] == 'dark'):
                if len(rawDataSet.darkFns) < config.params['max_n_darks']:
                    if dark_criteria_met(rawDataSet.objHeader, header):
                        rawDataSet.darkFns.append(filename)
        rawDataSet.flatFns.sort()
        rawDataSet.darkFns.sort()
        
    rawDataSets.sort(key=lambda x: x.baseNames['A'])
                 
    rawDataSetsWithFlats = []
    
    # remove data sets for which no flats are available
    for rawDataSet in rawDataSets:
        if len(rawDataSet.flatFns) == 0:
            logger.debug('no flats for {}'.format(rawDataSet.baseNames['A']))
        else:
            logger.debug('{} has {} flats'.format(rawDataSet.baseNames['A'], len(rawDataSet.flatFns)))
            rawDataSetsWithFlats.append(rawDataSet)
            
    if len(rawDataSets) < 1:
        logger.warning('no raw data sets could be assembled')
        return None
            
    # sort raw data set list by ascending UTC
    rawDataSetsWithFlats.sort(key=lambda x:x.objHeader['UTC'])
    
    logger.info('n raw data sets with flats = {}'.format(len(rawDataSetsWithFlats)))
 
    return(rawDataSetsWithFlats)


def get_headers(in_dir):
    """
    Makes a list of FITS files found in in_dir, unzips them as needed and.
    returns a dictionary of headers indexed by file name.
    
    param
        in_dir Input directory path.
        
    return
        Dictionary of headers indexed by file name.
        
    """
    
    headers = dict()
        
    for filename in os.listdir(in_dir):
        full_filename = in_dir + '/' + filename
        if fnmatch.fnmatch(filename, '*.fits*'):
            if config.params['gunzip'] is True and filename.endswith('gz'):
                os.system('gunzip ' + full_filename)
                full_filename = full_filename.rstrip('.gz')
            headers[full_filename] = fits.getheader(full_filename)
        
    return headers

def obj_criteria_met(header, failed2reduce):
    """
    Takes an object frame header and determines if it is a frame that can be 
    reduced by the DRP.
    
    param
        header Object file header.
        
    return
        True if the object file can be reduced by the DRP, False otherwise
        
    """

#    if header['DETMODE'].lower != 'spec':
#        failed2reduce['dmode'] += 1
#        return False
    if config.params['cmnd_line_mode'] == False:
        if header['IMAGETYP'].lower() != 'object':
#           failed2reduce['itype'] += 1
            return False
    if config.params['cmnd_line_mode'] == False:
        if header['DISPERS'].lower() != 'high':
            failed2reduce['dispmode'] += 1
            return False
    if header['NAXIS1'] != constants.N_COLS:
        failed2reduce['n1'] += 1
        return False
    if header['NAXIS2'] != constants.N_ROWS:
        failed2reduce['n2'] += 1
        return False
    if header['FILNAME'].upper() not in constants.filter_names:
        failed2reduce['fil'] += 1
        return False
    return True
    
def flat_criteria_met(obj_header, flat_header, ignore_dispers=False):
    """Takes an object frame header and a flat frame header and determines if 
    the flat satisfies the criteria for association with the object frame
    
    params
        obj_header Object frame header.
        flat_header Flat frame header.
        
    return
        True if the flat corresponds to the object frame, False otherwise.
        
    """
    eq_kwds = ['disppos', 'echlpos', 'filname', 'slitname']
    if ignore_dispers is not True:
        eq_kwds.append('dispers')
    for kwd in eq_kwds:
        if obj_header[kwd] != flat_header[kwd]:
            return False
    return True

def is_valid_pair(obj_A_header, obj_B_header):
    """Determines if the frames associated with the two headers passed as arguments
    comprise a valid AB pair.
    
    args:
        obj_A_header: FITS header of A
        obj_B_header: FITS header of B
        
    Returns:
        True if valid pair.
        False if not valid pair.
    """
    kwds = ['disppos', 'echlpos', 'filname', 'slitname', 'dispers', 'itime']
    for kwd in kwds:
        if obj_A_header[kwd] != obj_B_header[kwd]:
            return False
    return True

def dark_criteria_met(obj_header, dark_header):
    """
    Takes an object frame header and a dark field frame header and determines if 
    the dark satisfies the criteria for association with the object frame
    
    params
        obj_header Object frame header.
        dark_header Dark frame header.
        
    return
        True if the dark corresponds to the object frame, False otherwise.
        
    """
    eq_kwds = ['elaptime']
    for kwd in eq_kwds:
        if obj_header[kwd] != dark_header[kwd]:
            return False
    return True
