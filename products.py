import matplotlib

#matplotlib.use('Agg')
    
import pylab as pl
import logging
import os
import errno
import warnings

try:
    import Image
except ImportError:
    from PIL import Image
    
import numpy as np
from astropy.io import fits
from skimage import exposure
import image_lib
import config
import nsdrp

warnings.filterwarnings('ignore')

main_logger = logging.getLogger('main')
obj_logger = logging.getLogger('obj')
file_count = [0]

#import ReducedDataSet
#from statsmodels.sandbox.regression.ols_anova_original import products

#def produceProfileFitsTable(name, order,profile):
#    continue

saveFitsTables = True
saveAsciiTables = True
savePlots = True
showPlots = False
saveJpgs = False

# This dictionary maps data product filename suffix (e.g. flux_tbl.fits)
# to output subdirectory (e.g. fitstbl/flux).
subdirs = dict([
                ('flux_tbl.fits',       'fitstbl/flux'      ),
                ('flux.txt',            'ascii/flux'        ),
                ('flux.fits',           'fits/flux'         ),
                ('flux.png',            'previews/flux'     ),
                ('profile_tbl.fits',    'fitstbl/profile'   ),
                ('profile.txt',         'ascii/profile'     ),
                ('profile.fits',        'fits/profile'      ),
                ('profile.png',         'previews/profile'  ),
                ('wavecal_tbl.fits',    'fitstbl/wavecal'   ),
                ('wavecal.txt',         'ascii/wavecal'     ),
                ('order.fits',          'fits/order'        ),
                ('order.png',           'previews/order'    ),
                ('sky.fits',            'fits/sky'          ),
                ('sky.png',             'previews/sky'      ),
                ('snr.fits',            'fits/snr'          ),
                ('snr.png',             'previews/snr'      ),
                ('trace.fits',          'fits/trace'        ),
                ('trace.png',           'previews/trace'    ),
#               ('spectra.png',     'previews/spectra'  ),
               ])


def constructFileName(outpath, base_name, order, fn_suffix):
    """
    Constructs data product filename, including directory path given:
    outpath - data product root directory path
    base_name - base object file name, e.g. NS.20000325.49894
    order - order number, or None if this is not a per-order file
    fn_suffix - filename suffix, e.g. flux_tbl.fits
    """
    fn = outpath + '/' + subdirs[fn_suffix] + '/' + base_name + '_' + fn_suffix
    if order is None:
        return fn
    else:
        return fn[:fn.rfind(base_name) + len(base_name) + 1] + str(order) + fn[fn.rfind(fn_suffix)-1:]
   
def log_fn(fn):  
        #obj_logger.info('saving {}'.format(fn))
        file_count[0] += 1 
        return

def gen(reduced, out_dir):
    """
    Given a ReducedDataSet object and a root output directory, generate all
    data products and store results in output directory and subdirectories.
    """
    obj_logger.info('generating data products for {}...'.format(reduced.getBaseName()))
    
    file_count[0] = 0
    
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    
    # make sub directories
    for k, v in subdirs.iteritems():
        try:
            os.makedirs(out_dir + '/' + v)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                main_logger.critical('failed to create output directory ' + v)
                obj_logger.critical('failed to create output directory ' + v)
                raise IOError('failed to create output directory ' + v)
            
    # prepare extended fits header
    header = reduced.header
    header['NSDRPVER'] = (nsdrp.VERSION, 'NSDRP Version')
    header['COMMENT'] = ('NSDRP', 'NSDRP')
    if reduced.frameCalRmsRes is not None: 
        header['WFITRMS'] = (round(reduced.frameCalRmsRes, 4), 
                'RMS per-frame wavelength fit residual')
    if reduced.frameCalCoeffs is not None:
        for i in range(6):
            header['WFIT{}'.format(i)] = (round(reduced.frameCalCoeffs[i], 6), 
                'wavelength fit coefficient {}'.format(i))
    header['DARK'] = (reduced.darkKOAId, 'KOAID of dark frame or none')
    for i in range(len(reduced.flatKOAIds)):
        header['FLAT' + str(i)] = (reduced.flatKOAIds[i], 'KOAID of flat {}'.format(i))
    
    #
    # produce per-frame data products
    # 
    
    # construct arrays for global wavelength cal table
    order_num = []
    col = []
    source = []
    wave_exp = []
    wave_fit = []
    res = []
    peak = []
    slope = []
    
    for order in reduced.orders:
        for line in order.lines:
            if line.frameFitOutlier == False:
                order_num.append(order.flatOrder.orderNum)
                col.append(line.centroid)
                source.append('sky')
                wave_exp.append(line.waveAccepted)
                wave_fit.append(line.frameWaveFit)
                res.append(abs(line.frameWaveFit - line.waveAccepted))
                peak.append(line.peak)
                slope.append(line.frameFitSlope)
                
    # global wavelength cal tables
    
    wavelengthCalAsciiTable(
            out_dir, reduced.getBaseName(), order_num, col, source, wave_exp, wave_fit, res, peak, slope)
    
    wavelengthCalFitsTable(
            out_dir, reduced.getBaseName(), order_num, col, source, 
            wave_exp, wave_fit, res, peak, slope)
 
    #       
    # produce per-order data products
    #
    
    SKYDIST_INFO = 'Distance between sky and object extraction windows'
    
    for order in reduced.orders:
            
        # extend header further with per-order data
        header['FLATSCAL'] = (round(order.flatOrder.mean, 5), 
                'flat field normalization scale factor')
        header['ECHLORD'] = (order.flatOrder.orderNum, 
                'Echelle order number')
        if order.isPair:
            obj_window_w = len(order.objWindow['AB'])
        else:
            obj_window_w = len(order.objWindow['A'])

        header['OBJEXTRW'] = (obj_window_w, 'width of object extraction window in pixels')
         
        if order.isPair:
            centroid = order.centroid['AB']
            header['ORDERSNR'] = (round(order.snr['AB'], 3), 
                    'signal-to-noise ratio for order')
        else:
            centroid = order.centroid['A']
            header['ORDERSNR'] = (round(order.snr['A'], 3), 
                    'sign-to-noise ratio for order')
            header['SKYEXTRW'] = (max(len(order.topSkyWindow), len(order.botSkyWindow)), 
                    'width of sky subtraction window in pixels')
            try:
                if len(order.topSkyWindow) > 0:
                    header['SKYDIST'] = (order.topSkyWindow['A'][0] - order.objWindow['A'][-1], SKYDIST_INFO)
                else:
                    header['SKYDIST'] = (order.objWindow['A'][0] - order.botSkyWindow['A'][-1], SKYDIST_INFO)
            except IndexError as e:
                header['SKYDIST'] = 'Undefined'
                      
                      
        header['PROFPEAK'] = (round(centroid, 3), 'fractional row number of profile peak')
            

        #
        # flux ASCII and FITS tables
        #
        if order.isPair:
            obj_spec = order.objSpec['AB']
            noise_spec = order.noiseSpec['AB']
            snr_spec = np.absolute(order.objSpec['AB'] / order.noiseSpec['AB'])
        else:
            obj_spec = order.objSpec['A']
            noise_spec = order.noiseSpec['A']
            snr_spec = np.absolute(order.objSpec['A'] / order.noiseSpec['A'])
            
        fluxAsciiTable(out_dir, reduced.getBaseName(), order.flatOrder.orderNum, order.waveScale, 
                obj_spec, order.skySpec['A'], order.synthesizedSkySpec, snr_spec,
                order.flatSpec, order.flatOrder.topEdgeTrace, order.flatOrder.botEdgeTrace, 
                order.flatOrder.avgEdgeTrace, order.flatOrder.smoothedSpatialTrace, 
                order.flatOrder.smoothedSpatialTrace - order.flatOrder.avgEdgeTrace)
              
        fluxFitsTable(out_dir, reduced.getBaseName(), order.flatOrder.orderNum, order.waveScale, 
                obj_spec, order.skySpec['A'], order.synthesizedSkySpec, noise_spec,
                order.flatSpec, order.flatOrder.topEdgeTrace, order.flatOrder.botEdgeTrace, 
                order.flatOrder.avgEdgeTrace, order.flatOrder.smoothedSpatialTrace, 
                order.flatOrder.smoothedSpatialTrace - order.flatOrder.avgEdgeTrace)

        #
        # spatial trace plots
        # 
        tracePlot(out_dir, reduced.getBaseName(), order.flatOrder.orderNum, order.flatOrder.avgEdgeTrace, 
                order.flatOrder.smoothedSpatialTrace, order.flatOrder.spatialTraceMask)
     
        traceFits(out_dir, reduced.getBaseName(), order.flatOrder.orderNum, order.flatOrder.avgEdgeTrace)
        
        #
        # spatial profile per order plot, tables and 1-d FITS file
        #
        for frame in order.frames:
 
            if frame == 'AB':
                topSkyWindow = None
                botSkyWindow = None
                topBgMean = None
                botBgMean = None
            else:
                topSkyWindow = order.topSkyWindow[frame]
                botSkyWindow = order.botSkyWindow[frame]
                topBgMean = order.topBgMean[frame]
                botBgMean = order.botBgMean[frame]
            
            profilePlot(out_dir, reduced.baseNames[frame], order.flatOrder.orderNum, 
                    order.spatialProfile[frame], order.peakLocation[frame], order.centroid[frame], 
                    order.objWindow[frame], topSkyWindow, botSkyWindow, topBgMean, botBgMean, 
                    order.gaussianParams[frame], order.snr[frame])
                
            profileAsciiTable(out_dir, reduced.baseNames[frame], order.flatOrder.orderNum, 
                    order.spatialProfile[frame])
             
            profileFitsTable(out_dir, reduced.baseNames[frame], order.flatOrder.orderNum, 
                    order.spatialProfile[frame])
            
            profileFits(out_dir, reduced.baseNames[frame], order.flatOrder.orderNum, 
                    order.spatialProfile[frame], header)
            

        #
        # flux spectrum plot and 1-d FITS file
        #   
        for frame in reduced.frames:
            
            spectrumPlot(out_dir, reduced.baseNames[frame], 'flux', order.flatOrder.orderNum, 
                'counts', order.objSpec[frame], order.waveScale, order.calMethod)
             
            fitsSpectrum(out_dir, reduced.baseNames[frame], 'flux', order.flatOrder.orderNum, 
                'counts', order.objSpec[frame], order.waveScale, header)

        #
        # sky spectrum plot and 1-d FITS file
        #
        if reduced.isPair:
            frames = ['A', 'B']
        else:
            frames = ['A']
        
        for frame in frames:
            spectrumPlot(out_dir, reduced.baseNames[frame], 'sky', order.flatOrder.orderNum, 
                'counts', order.skySpec[frame], order.waveScale, order.calMethod)
    
            fitsSpectrum(out_dir, reduced.baseNames[frame], 'sky', order.flatOrder.orderNum, 
                'counts', order.skySpec[frame], order.waveScale, header)
        
        #
        # noise spectrum
        #
        for frame in reduced.frames:

            spectrumPlot(out_dir, reduced.baseNames[frame], 'snr', order.flatOrder.orderNum, 
                '', np.absolute(order.objSpec[frame]/order.noiseSpec[frame]), order.waveScale, 
                order.calMethod)
              
            fitsSpectrum(out_dir, reduced.baseNames[frame], 'snr', order.flatOrder.orderNum, 
                '', np.absolute(order.objSpec[frame]/order.noiseSpec[frame]), order.waveScale, 
                header)


        #
        # rectified order plot 2-d image plot and FITS file
        #
        for frame in order.frames:
            twoDimOrderPlot(out_dir, order.baseNames[frame], 'rectified order image', 'order.png', 
                order.flatOrder.orderNum, order.ffObjImg[frame], order.waveScale, order.calMethod)
            twoDimOrderFits(out_dir, order.baseNames[frame], order.flatOrder.orderNum, 
                order.ffObjImg[frame], header)
        
        # end of for each order


    main_logger.info('n data products generated for {} = {}'.format(
            reduced.getBaseName(), str(file_count[0])))
    
    obj_logger.info('n data products generated for {} = {}'.format(
            reduced.getBaseName(), str(file_count[0])))
    return 
    
    
def tracePlot(outpath, base_name, order_num, raw, fit, mask):

    pl.figure("Trace Plot", figsize=(6, 5), facecolor='white')
    pl.title('trace, ' + base_name + ", order " + str(order_num), fontsize=14)
    pl.xlabel('column (pixels)')
    pl.ylabel('row (pixels)')
    
#     yrange = offraw.max() - offraw.min()
 
    x = np.arange(raw.shape[0])
    
    pl.plot(x[mask], raw[mask], "ko", mfc="none", ms=1.0, linewidth=1, label="derived")
    pl.plot(x, fit, "k-", mfc="none", ms=1.0, linewidth=1, label="fit")
        
    pl.plot(x[np.logical_not(mask)], raw[np.logical_not(mask)], 
        "ro", mfc="red", mec="red", ms=2.0, linewidth=2, label="ignored")

    
    rms = np.sqrt(np.mean(np.square(raw - fit)))
    pl.annotate('RMS residual = ' + "{:.3f}".format(rms), (0.3, 0.8), xycoords="figure fraction")
    
    pl.minorticks_on()
    pl.grid(True)
    pl.legend(loc='best', prop={'size': 8})

    fn = constructFileName(outpath, base_name, order_num, 'trace.png')
    savePreviewPlot(fn)
    pl.close()
    log_fn(fn)
    
    return

def traceFits(outpath, base_name, order_num, trace):
    hdu = fits.PrimaryHDU(trace)
    hdulist = fits.HDUList(hdu)
    fn = constructFileName(outpath, base_name, order_num, 'trace.fits')     
    hdulist.writeto(fn, clobber=True)
    log_fn(fn)
    return

# def profileAsciiTable(outpath, base_name, order_num, profile):
#     
#     names = ['row', 'mean_flux'] 
#     units = ['pixels', 'counts']
#     formats = ['.0f', '.0f']
#     widths = [6, 9] 
#     
#     buff = []
#     
#     line = []
#     for i, name in enumerate(names):
#         line.append('{:>{w}}'.format(name, w=widths[i]))
#     buff.append('| {} |'.format('| '.join(line)))
#     
#     line = []
#     for i, unit in enumerate(units):
#         line.append('{:>{w}}'.format(unit, w=widths[i]))
#     buff.append('| {} |'.format('| '.join(line)))
#         
#     if UNDERLINE:
#         line = []
#         for i, unit in enumerate(units):
#             line.append('{:->{w}}'.format('', w=widths[i]))
#         buff.append('--'.join(line))
#     
#     for row in range(profile.shape[0]):
#         data = [row, profile[row]]
#         line = []
#         for i, val in enumerate(data):
#             line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
#         buff.append('  {}  '.format('  '.join(line)))
#                 
#     fn = constructFileName(outpath, base_name, order_num, 'profile.txt')
#     fptr = open(fn, 'w')
#     fptr.write('\n'.join(buff))
#     fptr.close()
#     log_fn(fn)
#  
#     return


def profileAsciiTable(outpath, base_name, order_num, profile):
    
    names = ['row', 'mean_flux'] 
    units = ['pixels', 'counts']
    formats = ['.0f', '.0f']
    widths = [6, 9] 
    
    buff = []
    
    if config.params['pipes'] is True:
        p_char = '|'
    else:
        p_char = ' '
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('{} {} {}'.format(p_char, (p_char + ' ').join(line), p_char))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('{} {} {}'.format(p_char, (p_char + ' ').join(line), p_char))
    
    for row in range(profile.shape[0]):
        data = [row, profile[row]]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
                
    fn = constructFileName(outpath, base_name, order_num, 'profile.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
    log_fn(fn)
 
    return

 
def profileFitsTable(outpath, base_name, order_num, profile):  
    
    prihdr = fits.Header()
    prihdr['COMMENT'] = "profile table"
    prihdu = fits.PrimaryHDU(header=prihdr)
    tbhdu = fits.new_table(
        fits.ColDefs([
            fits.Column(name='row (pix)', format='1I', array=np.arange(profile.shape[0], dtype=int)),
            fits.Column(name='mean_flux (cnts)', format='1D', array=profile)]))
    thdulist = fits.HDUList([prihdu, tbhdu]) 
    fn = constructFileName(outpath, base_name, order_num, 'profile_tbl.fits')
    thdulist.writeto(fn, clobber=True)  
    log_fn(fn)
    return  

def profileFits(outpath, base_name, order_num, profile, header):
    hdu = fits.PrimaryHDU(profile)
    hdulist = fits.HDUList(hdu)
    hdr = hdulist[0].header
    for k, v in header.iteritems():
        try:
            hdr[k] = v
        except Exception as e:
            #obj_logger.warning(e.message)
            pass

            
    fn = constructFileName(outpath, base_name, order_num, 'profile.fits')     
    hdulist.writeto(fn, clobber=True)
    log_fn(fn)
    return
    
def profilePlot(outpath, base_name, order_num, profile, peak, centroid,
            ext_range, sky_range_top, sky_range_bot, top_bg_mean, bot_bg_mean, gaussian, snr):

    pl.figure('spatial profile', facecolor='white')
    pl.cla()
    pl.title('spatial profile, ' + base_name + ', order ' + str(order_num), fontsize=14)

    pl.xlabel('relative row (pixels)')
    pl.ylabel('flux (counts)')

    # set axes limits
    yrange = profile.max() - profile.min()
    ymin = profile.min() - (0.1 * yrange)
    ymax = profile.max() + (0.1 * yrange)
    pl.ylim(ymin, ymax)
    
    # set ticks and grid
    pl.minorticks_on()
    pl.grid(False)
    
    # plot profile
    pl.plot(profile, "ko-", mfc='none', ms=3.0, linewidth=1)
    
    # plot vertical line to indicate location of centroid
    pl.plot([centroid, centroid], [ymin, profile.max()], "g-", linewidth=0.5, label='peak')
   
    # draw extraction window
    wvlh = 0.01 * yrange;
    ewh = 0.05 * yrange
    
    pl.plot((ext_range[0], ext_range[-1]), (profile.max(), profile.max()), 
            'r', linewidth=0.5, label='extraction window')
    pl.plot((ext_range[0], ext_range[0]),
            (profile.max() - wvlh, profile.max() + wvlh), 
            'r', linewidth=1.5)
    pl.plot((ext_range[-1], ext_range[-1]),
            (profile.max() - wvlh, profile.max() + wvlh), 
            'r', linewidth=1.5)  
    
    # draw annotation showing centroid, width and SNR
    if gaussian is None:
        width = 'unknown'
    else:
        width = '{:.1f}'.format(abs(gaussian[2]))
    
    aStr = 'centroid = {:.1f} pixels'.format(centroid)
    aStr = aStr + '\nwidth = {} pixels'.format(width)
    if snr is not None:
        aStr = aStr + '\nSNR = {:.1f}'.format(snr)
    
    if peak > (len(profile) / 2):
        pl.annotate(aStr, (peak - (len(ext_range) / 2) - 20, ((ymax - ymin) * 3 / 5) + ymin))
    else:
        pl.annotate(aStr, (peak + (len(ext_range) / 2) + 5, ((ymax - ymin) * 3 / 5) + ymin))
        
    # draw sky windows
 
    if sky_range_top is not None and len(sky_range_top) > 0:

        pl.plot((sky_range_top[0], sky_range_top[-1]),
                (top_bg_mean, top_bg_mean), 
                'b', linewidth=0.5, label='sky window')  
        pl.plot((sky_range_top[0], sky_range_top[0]),
                (top_bg_mean - wvlh, top_bg_mean + wvlh), 
                'b', linewidth=1.5)
        pl.plot((sky_range_top[-1], sky_range_top[-1]),
                (top_bg_mean - wvlh, top_bg_mean + wvlh), 
                'b', linewidth=1.5)  
        
    if sky_range_bot is not None and len(sky_range_bot) > 0:
        
        pl.plot((sky_range_bot[0], sky_range_bot[-1]),
                (bot_bg_mean, bot_bg_mean), 
                'b', linewidth=0.5)   
        pl.plot((sky_range_bot[0], sky_range_bot[0]),
                (bot_bg_mean - wvlh, bot_bg_mean + wvlh), 
                'b', linewidth=1.5)
        pl.plot((sky_range_bot[-1], sky_range_bot[-1]),
                (bot_bg_mean - wvlh, bot_bg_mean + wvlh), 
                'b', linewidth=1.5)
        
    # draw best fit Gaussian
    if gaussian is not None:
        min = np.amin(profile)
        if min < 0:
            offset = 0
        else:
            offset = min
        pl.plot(image_lib.gaussian(range(len(profile)), gaussian[0], gaussian[1], gaussian[2]) + offset, 
                'k--', linewidth=0.5, label='Gaussian fit')
        
    pl.legend(loc='best', prop={'size': 8})
    
    fn = constructFileName(outpath, base_name, order_num, 'profile.png')
    savePreviewPlot(fn)
    pl.close()
    log_fn(fn)
    return

def fluxAsciiTable(outpath, base_name, order_num, wave, flux, sky, synth_sky, snr,
                     flat, trace_upper, trace_lower, trace_mean, trace_fit, fit_res):
     
    names = [   'col',         'wave',         'flux',         'sky', 
                'synth_sky',   'snr',          'flat',        'trace_upper', 
                'trace_lower', 'trace_mean',   'trace_fit',    'fit_res']
    units = [   'pixels',      'Angstroms',    'counts',       'counts', 
                '',            '',             '',             'pixel', 
                'pixels',      'pixels',       'pixels',       'pixels']
    formats = [ 'd',            '.6e',          '.3e',          '.3e',
                '.3e',          '.3e',          '.3e',          '.3e',
                '.3e',          '.3e',          '.3e',          '.3e']
    nominal_width = 10
    widths = []
    
    if config.params['pipes'] is True:
        p_char = '|'
    else:
        p_char = ' '
        
    if trace_lower is None:
        trace_lower = np.zeros(1024, dtype=float)
    if trace_mean is None:
        trace_mean = np.zeros(1024, dtype=float)
        
    for name in names:
        widths.append(max(len(name), nominal_width))
            
    for i in range(len(names)):
        widths[i] = max(len(units[i]), widths[i])
        
    widths[0] = 6
    widths[1] = 13
    
    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('{} {} {}'.format(p_char, (p_char + ' ').join(line), p_char))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('{} {} {}'.format(p_char, (p_char + ' ').join(line), p_char))
    
    for col in range(wave.shape[0]):
        data = [col, wave[col], flux[col], sky[col], synth_sky[col], snr[col], flat[col],
                trace_upper[col], trace_lower[col], trace_mean[col], trace_fit[col], fit_res[col]]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
        
    #print('\n'.join(buff))
        
    fn = constructFileName(outpath, base_name, order_num, 'flux.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
    log_fn(fn)
    return

    
def fluxFitsTable(outpath, base_name, order_num, wave, flux, sky, synth_sky, error,
                     flat, trace_upper, trace_lower, trace_mean, trace_fit, fit_res):
    # create binary FITS table
    prihdr = fits.Header()

    prihdu = fits.PrimaryHDU(header=prihdr)

    tbhdu = fits.BinTableHDU.from_columns([
                fits.Column(name='col', format='1I', array=np.arange(1024, dtype=int)),
                fits.Column(name='wave (A)', format='1D', array=wave),
                fits.Column(name='flux (cnts)', format='1D', array=flux),
                fits.Column(name='noise (cnts)', format='1D', array=error),
                fits.Column(name='sky (cnts)', format='1D', array=sky),
                fits.Column(name='synth_sky', format='1D', array=synth_sky),
                fits.Column(name='sig_to_noise', format='1D', array=np.absolute(flux/error)),
                fits.Column(name='flat (cnts)', format='1D', array=flat),
                fits.Column(name='trace_upper (pix)', format='1D', array=trace_upper),
                fits.Column(name='trace_lower (pix)', format='1D', array=trace_lower),
                fits.Column(name='trace_mean (pix)', format='1D', array=trace_mean),
                fits.Column(name='trace_fit (pix)', format='1D', array=trace_fit),
                fits.Column(name='fit_res (pix)', format='1D', array=trace_mean - trace_fit)])
    fn = constructFileName(outpath, base_name, order_num, 'flux_tbl.fits')   
    tbhdu.writeto(fn, clobber=True)
    log_fn(fn)
    return

def spectrumPlot(outpath, base_name, title, order_num, y_units, cont, wave, 
                 wave_note='unknown'):
    
    pl.figure(title, facecolor='white')
    pl.clf()
    pl.title(title + ', ' + base_name + ", order " + str(order_num), fontsize=12)
    pl.xlabel('wavelength ($\AA$) (' + wave_note + ')')
    if len(y_units) > 0:
        pl.ylabel(title + '(' + y_units + ')')
    else:
        pl.ylabel(title)
    pl.grid(True)
    pl.plot(wave[:1004], cont[:1004], "k-", mfc="none", ms=3.0, linewidth=1)
    
    ymin, ymax = pl.ylim()
    pl.plot(wave, cont, "k-", mfc="none", ms=3.0, linewidth=1)
    pl.ylim(ymin, ymax)
    
#     axes = pl.gca()
    #axes.set_ylim(0, 300)
    
    fn = constructFileName(outpath, base_name, order_num, title + '.png')
    savePreviewPlot(fn)
    log_fn(fn)
    return
    
def fitsSpectrum(outpath, base_name, title, order_num, y_units, cont, wave, header):
    
    hdu = fits.PrimaryHDU(cont)
    hdulist = fits.HDUList(hdu)
    hdr = hdulist[0].header
    for k, v in header.iteritems():
        try:
            hdr[k] = v
        except Exception as e:
            #obj_logger.warning(e.message)
            pass
            
    fn = constructFileName(outpath, base_name, order_num, title + '.fits')     
    hdulist.writeto(fn, clobber=True)
    log_fn(fn)
    return
    
    
def multiSpectrumPlot(outpath, base_name, order, y_units, cont, sky, noise, wave):
    
    title = 'spectrum'
    pl.figure(title, facecolor='white')
    pl.clf()
    pl.title(title + ', ' + base_name + ", order " + str(order), fontsize=12)
    pl.xlabel('wavelength ($\AA$)')
    pl.ylabel(title + '(' + y_units + ')')
    pl.grid(True)
    
    pl.plot(wave[:1004], cont[:1004], "k-", mfc="none", ms=3.0, linewidth=1, label='object')
    pl.plot(wave[:1004], sky[:1004], "b-", mfc="none", ms=3.0, linewidth=1, label='sky')
    pl.plot(wave[:1004], noise[:1004], "r-", mfc="none", ms=3.0, linewidth=1, label='noise (1 sigma)')

    pl.legend(loc='best', prop={'size': 8})

    fn = constructFileName(outpath, base_name, order, 'spectra.png')
        
    savePreviewPlot(fn)
    log_fn(fn)
    pl.close()
    return
    
#     hdu = fits.PrimaryHDU(cont)
#     hdulist = fits.HDUList(hdu)
#     hdulist[0].header['object'] = '511 Davida'
#     fits_fn = outpath + '/' + base_name + '_' + str(order) + '_' + title + '.fits'
#     logger.info('writing ' + title + ' plot for order ' + str(order) + ' to ' + fits_fn)
#     hdulist.writeto(fits_fn, clobber=True)
    
def wavelengthCalAsciiTable(outpath, base_name, order, col, source, wave_exp, wave_fit, res, peak, 
        slope):
    
    names = ['order', 'source', 'col', 'wave_exp', 'wave_fit', 'res', 'peak', 'disp'] 
    units = ['', '', 'pixels', 'Angstroms', 'Angstroms',  'Angstroms', 'counts', 'Angstroms/pixel']
    formats = [ 'd', '', '.3f', '.6e', '.6e', '.3f', 'd', '.3e']
    nominal_width = 10
    widths = []
    
    if config.params['pipes'] is True:
        p_char = '|'
    else:
        p_char = ' '
        
    for name in names:
        widths.append(max(len(name), nominal_width))
            
    for i in range(len(names)):
        widths[i] = max(len(units[i]), widths[i])
       
    widths[0] = 6 
    widths[1] = 6
    widths[2] = 8
    widths[3] = 13
    widths[4] = 13
    widths[5] = 9
    widths[6] = 6

    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('{} {} {}'.format(p_char, (p_char + ' ').join(line), p_char))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('{} {} {}'.format(p_char, (p_char + ' ').join(line), p_char))
    
    for i in range(len(order)):
        data = [order[i], source[i], col[i], wave_exp[i], wave_fit[i], res[i], (int)(peak[i]), slope[i]]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
        
    #print('\n'.join(buff))
        
    fn = constructFileName(outpath, base_name, None, 'wavecal.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
    log_fn(fn)
    return



def wavelengthCalFitsTable(outpath, base_name, order, col, source, wave_exp, wave_fit, res, peak, slope):
    prihdr = fits.Header()
    prihdr['COMMENT'] = "wavelength calibration table"
    prihdu = fits.PrimaryHDU(header=prihdr)
    tbhdu = fits.new_table(
            fits.ColDefs([
                fits.Column(name='order', format='1I', array=order),
                fits.Column(name='source', format='1A', array=source),
                fits.Column(name='col (pixels)', format='1D', array=col),
                fits.Column(name='wave_exp (Angstroms)', format='1D', array=wave_exp),
                fits.Column(name='wave_fit (Angstroms)', format='1D', array=wave_fit),
                fits.Column(name='res (Angstroms)', format='1D', array=res),
                fits.Column(name='peak (counts)', format='1D', array=peak),
                fits.Column(name='disp (Angstroms/pixel)', format='1D', array=slope)]))
    thdulist = fits.HDUList([prihdu, tbhdu]) 
    fn = constructFileName(outpath, base_name, None, 'wavecal_tbl.fits')   
    thdulist.writeto(fn, clobber=True)         
    log_fn(fn)
    return
    
def twoDimOrderPlot(outpath, base_name, title, base_filename, order_num, data, x_scale,
            wave_note='unknown'):
    """
    Produces a generic 2-d image plot.
    
    Arguments:
        output: Directory path of root products directory.
        base_name: Base name of object frame.
        title: Title of plot, e.g. rectified order image.
        base_filename:
        order_num:
        data:
        x_scale:
    """
    pl.figure('2d order image', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.title(title + ', ' + base_name + ", order " + str(order_num), fontsize=14)
    pl.xlabel('wavelength($\AA$) (' + wave_note + ')', fontsize=12)
    pl.ylabel('row (pixel)', fontsize=12)
    
    pl.imshow(exposure.equalize_hist(data), origin='lower', 
                  extent=[x_scale[0], x_scale[-1], 0, data.shape[0]], aspect='auto')      
       
#    pl.colorbar()
#    pl.set_cmap('jet')
    pl.set_cmap('gray')
#     pl.set_cmap('Blues_r')

    fn = constructFileName(outpath, base_name, order_num, base_filename)
    savePreviewPlot(fn)
    log_fn(fn)
    pl.close()
        
    return
    

def twoDimOrderFits(outpath, base_name, order_num, data, header):     
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList(hdu)
    hdr = hdulist[0].header
    for k, v in header.iteritems():
        try:
            hdr[k] = v
        except Exception as e:
            #obj_logger.warning(e.message)
            pass
            
    fn = constructFileName(outpath, base_name, order_num, 'order.fits')
    hdulist.writeto(fn, clobber=True)
    log_fn(fn)
    return


def savePreviewPlot(fn):
    """
    Saves a preview plot in either PNG or JPG format, depending on configuration parameter.
    
    If config.params['jpg'] is false then the plot is saved in PNG format which is the default.
    If config.params['jpg'] is true then the PNG file is converted to JPG format and the PNG file
    is deleted.
    
    Args:
        fn: filename for PNG format, e.g. out/previews/flux/KOAID.png
        
    """
    pl.savefig(fn)
    if (config.params['jpg']):
        outfn = fn[:fn.find('.png')] + '.jpg'
        Image.open(fn).save(outfn, 'JPEG')
        os.remove(fn)
    
        
