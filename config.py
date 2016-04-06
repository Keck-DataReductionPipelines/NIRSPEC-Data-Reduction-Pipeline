params = {}

# command line arguments
params['cmnd_line_mode']        = False
params['debug']                 = False     # comment
params['no_cosmic']             = False     # cosmic ray rejection inhibited if True
params['no_products']           = False     # data product generation inhibited if True
params['obj_window']            = 9         # comment
params['sky_window']            = 8         # comment
params['sky_separation']        = 2         # comment

params['oh_filename']           = './ir_ohlines.dat'
params['oh_envar_name']         = 'NSDRP_OH_FILENAME'
params['oh_envar_override']     = False     # if True then use params['oh_filename'] 
                                            # even if envar is set

params['int_c']                 = False
params['dgn']                   = False     # diagnostic data product generation enabled if True
params['npy']                   = False
params['verbose']               = False     # all log file messages printed to stdout if True
params['subdirs']               = False     # data products in per-frame subdirectory if True
params['lla']                   = 2         # sky line location algorithm
params['pipes']                 = False
params['shortsubdir']           = True
params['ut']                    = None
params['gunzip']                = False
params['out_dir']               = './nsdrp_out'       # used only in command line mode

# configuration and tuning parameters
params['max_n_flats']           = 8
params['max_n_darks']           = 8
params['max_spatial_trace_res'] = 1.0


params['long_slit_edge_margin'] = 6         # cut out margin in pixels
params['large_tilt_threshold']  = 20
params['large_tilt_extra_padding'] = 10
params['overscan_width']        = 10

# expected order number at bottom of detector
starting_order = {
        'NIRSPEC-1': 80, 
        'NIRSPEC-2': 70, 
        'NIRSPEC-3': 67, 
        'NIRSPEC-4': 61, 
        'NIRSPEC-5': 53, 
        'NIRSPEC-6': 49, 
        'NIRSPEC-7': 41 
}

def get_starting_order(filtername):
    return starting_order[filtername.upper()]

# order edge location error threshold
max_edge_location_error = {
        'NIRSPEC-1': 30, 
        'NIRSPEC-2': 50, 
        'NIRSPEC-3': 50, 
        'NIRSPEC-4': 20, 
        'NIRSPEC-5': 50, 
        'NIRSPEC-6': 20, 
        'NIRSPEC-7': 60 
}

def get_max_edge_location_error(filtername, slit):
    if '24' in slit:
        if 'NIRSPEC-7' in filtername.upper():
            return 35
        else:
            return 30
    else:
        return max_edge_location_error[filtername.upper()]
    
# order cutout padding
long_slit_cutout_padding = {
        'NIRSPEC-1': 0, 
        'NIRSPEC-2': 0, 
        'NIRSPEC-3': 0, 
        'NIRSPEC-4': 0, 
        'NIRSPEC-5': 0, 
        'NIRSPEC-6': 15, 
        'NIRSPEC-7': 30          
}
short_slit_cutout_padding = {
        'NIRSPEC-1': 0, 
        'NIRSPEC-2': 0, 
        'NIRSPEC-3': 10, 
        'NIRSPEC-4': 15, 
        'NIRSPEC-5': 15, 
        'NIRSPEC-6': 15, 
        'NIRSPEC-7': 30          
}

def get_cutout_padding(filtername, slit):
    if '24' in slit:
        return(long_slit_cutout_padding[filtername.upper()])
    else:
        return(short_slit_cutout_padding[filtername.upper()])
    