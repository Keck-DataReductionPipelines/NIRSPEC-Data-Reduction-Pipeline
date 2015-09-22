import os
import subprocess
import logging
from astropy.io import fits

import RawDataSet
from __builtin__ import False

def create(in_dir):
    """
    Given an input directory path, creates and returns a list of raw data sets.
    """
    
    logger = logging.getLogger('main')
    logger.info("creating raw data sets from files in " + in_dir)
    
    # get list of fits files
    headers = get_headers(in_dir)
            
    if (len(headers) == 0):
        logger.critical('no fits files found')
        return
    
    logger.info(str(len(headers)) + ' fits files found')
    
    rawDataSets = []

    for filename, header in headers.items():
        if (header['IMAGETYP'] == 'object'):
            rawDataSets.append(RawDataSet.RawDataSet(filename, header))
    
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
                         
#     for dataSet in rawDataSets:
#         print "obj: " + dataSet.objFileName
#         print "darks: " + str(dataSet.darkFileNames)
#         print "flats: " + str(dataSet.flatFileNames)
          
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