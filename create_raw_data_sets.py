import os
import subprocess
import logging
from astropy.io import fits

import RawDataSet
import nirspec_constants as constants

from __builtin__ import False


def create(in_dir):
    """
    Given an input directory path, creates and returns a list of raw data sets.
    A raw data set is an instance of the class RawDataSet and consists of an
    object file name, it's header, and associated flatas and darks.
    """
    
    logger = logging.getLogger('main')
    logger.debug("creating raw data sets from files in " + in_dir)
    
    # get list of fits files
    headers = get_headers(in_dir)
            
    if (len(headers) == 0):
        logger.critical('no fits files found')
        return
    
    logger.info(str(len(headers)) + ' fits files found')
    
    rawDataSets = []

    for filename, header in headers.items():
            if obj_criteria_met(header):
                rawDataSets.append(RawDataSet.RawDataSet(filename, header))
            else:
                logger.info('{} is in low dispersion mode, not reduced'.format(
                        filename[filename.rfind('/') + 1:]))
    
    logger.info(str(len(rawDataSets)) + " object frame(s) found")
    
    # associate darks and flats with each object frame
    for rawDataSet in rawDataSets:
        for filename, header in headers.items():
            if (header['IMAGETYP'] == 'flatlamp'):
                if flat_criteria_met(rawDataSet.objHeader, header):
                    rawDataSet.flatFileNames.append(filename)
            elif (header['IMAGETYP'] == 'dark'):
                if dark_criteria_met(rawDataSet.objHeader, header):
                    rawDataSet.darkFileNames.append(filename)
        rawDataSet.flatFileNames.sort()
        rawDataSet.darkFileNames.sort()
             
    # remove data sets for which no flats are available
    for rawDataSet in rawDataSets:
        if len(rawDataSet.flatFileNames) == 0:
            logger.info('no flats for {}'.format(
                    rawDataSet.objFileName[rawDataSet.objFileName.rfind('/') + 1 :]))
            rawDataSets.remove(rawDataSet)
          
    return(rawDataSets)


def get_headers(in_dir):
    """
    Make a list of fits files found in in_dir, unzip them as needed.
    Return a dictionary of headers indexed by file name
    """
    cmnd = "find " + in_dir + " -name \*fits\* | sort"
    filenames, err = subprocess.Popen([cmnd], stdout=subprocess.PIPE, shell=True).communicate()
    filenames = filter(None, filenames.split('\n'))
    
    headers = dict()
    
    for filename in filenames:
        if filename.endswith('gz'):
            os.system('gunzip ' + filename)
            filename = filename.rstrip('.gz')
        headers[filename] = fits.getheader(filename)

    return headers

def obj_criteria_met(header):
    if header['IMAGETYP'].lower() != 'object':
        return False
    if header['DISPERS'].lower() != 'high':
        return False
    if header['NAXIS1'] != constants.N_COLS:
        return False
    if header['NAXIS2'] != constants.N_ROWS:
        return False
    return True
    
def flat_criteria_met(obj_header, flat_header):
    eq_kwds = ['disppos', 'echlpos', 'filname', 'slitname', 'dispers']
    for kwd in eq_kwds:
        if obj_header[kwd] != flat_header[kwd]:
            return False
    return True


def dark_criteria_met(obj_header, dark_header):
    eq_kwds = ['elaptime']
    for kwd in eq_kwds:
        if obj_header[kwd] != dark_header[kwd]:
            return False
    return True