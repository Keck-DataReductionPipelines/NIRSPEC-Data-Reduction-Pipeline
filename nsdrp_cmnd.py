import logging

def process_frame(fn1, fn2, out_dir):
    
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

    if obj_header['ECHLPOS'] > 100:
        print('ERROR: cannot reduce low-resolution image (ECHLPOS > 100')
        exit(1)
        
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
                
    logger = logging.getLogger('main')

    # generate reduced data set by reducing raw data set
    reducedDataSet = reduce_frame.reduce(rawDataSet, out_dir)
        
    # produce data products from reduced data set
    if config.params['no_products'] is True:
        logger.info('data product generation inhibited by command line switch')
    else:
        products.gen(reducedDataSet, out_dir)

    # if diagnostics mode is enabled, then produce diagnostic data products
    if config.params['dgn'] is True:
        logger.info('diagnostic mode enabled, generating diagnostic data products')
        dgn.gen(reducedDataSet, out_dir)
        
    return        