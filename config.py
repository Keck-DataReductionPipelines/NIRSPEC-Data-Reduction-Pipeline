params = {}

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
params['max_n_flats']           = 8
params['max_n_darks']           = 8
params['out_dir']               = './nsdrp_out'       # used only in command line mode
params['max_spatial_trace_res'] = 1.0

