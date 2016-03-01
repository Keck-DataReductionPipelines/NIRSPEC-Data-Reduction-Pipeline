import logging
import os
import errno
import warnings
import numpy as np
import pylab as pl
from astropy.io import fits
from astropy.table import Table, Column
from astropy.io import ascii
from skimage import exposure
#import Image
import image_lib
#import PIL
import config
import warnings

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

# This dictionary maps data product filename suffix (e.g. flux.tbl)
# to output subdirectory (e.g. fitstbl/flux).
subdirs = dict([
                ('flux.tbl',        'fitstbl/flux'      ),
                ('flux.txt',        'ascii/flux'        ),
                ('flux.fits',       'fits/flux'         ),
                ('flux.png',        'previews/flux'     ),
                ('profile.tbl',     'fitstbl/profile'   ),
                ('profile.txt',     'ascii/profile'     ),
                ('profile.fits',    'fits/profile'      ),
                ('profile.png',     'previews/profile'  ),
                ('calids.tbl',      'fitstbl/cal'       ),
                ('calids.txt',      'ascii/cal'         ),
                ('order.fits',      'fits/order'        ),
                ('order.png',       'previews/order'    ),
                ('sky.fits',        'fits/sky'          ),
                ('sky.png',         'previews/sky'      ),
                ('noise.fits',      'fits/noise'        ),
                ('noise.png',       'previews/noise'    ),
                ('trace.fits',      'fits/trace'        ),
                ('trace.png',       'previews/trace'    ),
#                 ('spectra.png',     'previews/spectra'  ),
               ])


def constructFileName(outpath, base_name, order, fn_suffix):
    """
    Constructs data product filename, including directory path given:
    outpath - data product root directory path
    base_name - base object file name, e.g. NS.20000325.49894
    order - order number, or None if this is not a per-order file
    fn_suffix - filename suffix, e.g. flux.tbl
    """
    fn = outpath + '/' + subdirs[fn_suffix] + '/' + base_name + '_' + fn_suffix
    if order is None:
        return fn
    else:
        return fn[:fn.rfind('_') + 1] + str(order) + fn[fn.rfind('_'):]
   
def log_fn(fn):  
        #obj_logger.info('saving {}'.format(fn))
        file_count[0] += 1 
        return

def gen(reduced, out_dir):
    """
    Given a ReducedDataSet object and a root output directory, generate all
    data products and store results in output directory and subdirectories.
    """
    obj_logger.info('generating data products...')
    
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
    header['COMMENT'] = ('NSDRP', 'NSDRP')
    if reduced.rmsFitRes is not None: 
        header['WFITRMS'] = (round(reduced.rmsFitRes, 4), 'RMS wavelength fit residual')
    if reduced.coeffs is not None:
        for i in range(6):
            header['WFIT{}'.format(i)] = (round(reduced.coeffs[i], 6), 'wavelength fit coefficient {}'.format(i))
    header['DARK'] = (reduced.darkKOAId, 'KOAID of dark frame or none')
    for i in range(len(reduced.flatKOAIds)):
        header['FLAT' + str(i)] = (reduced.flatKOAIds[i], 'KOAID of flat {}'.format(i))
    
    
    # per frame data products
    # -----------------------
    
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
            if line.usedInGlobalFit:
                order_num.append(order.orderNum)
                col.append(line.centroid)
                source.append('sky')
                wave_exp.append(line.acceptedWavelength)
                wave_fit.append(line.globalFitWavelength)
                res.append(abs(line.globalFitWavelength - line.acceptedWavelength))
                peak.append(line.peak)
                slope.append(line.globalFitSlope)
                
    # global wavelength cal tables
    
    wavelengthCalAsciiTable(
            out_dir, reduced.baseName, order_num, col, source, wave_exp, wave_fit, res, peak, slope)
    
    wavelengthCalFitsTable(
            out_dir, reduced.baseName, order_num, col, source, 
            wave_exp, wave_fit, res, peak, slope)
 
            
    # per order data products
    
    for order in reduced.orders:
        
#         if reduced.coeffs is None:
#             wavelength_scale = order.wavelengthScaleMeas
#         else:
#             wavelength_scale = order.wavelengthScaleCalc
            
        # extend header further with per-order data
        header['FLATSCAL'] = (round(order.flatMean, 5), 'flat field normalization scale factor')
        header['ECHLORD'] = (order.orderNum, 'Echelle order number')
        header['OBJEXTRW'] = (len(order.objWindow), 'width of object extraction window in pixels')
        header['SKYEXTRW'] = (max(len(order.topSkyWindow), len(order.botSkyWindow)), 
                              'width of sky subtraction window in pixels')

        if len(order.topSkyWindow) > 0:
            header['SKYDIST'] = (order.topSkyWindow[0] - order.objWindow[-1], )
        else:
            header['SKYDIST'] = (order.objWindow[0] - order.botSkyWindow[-1], )
            
        header['PROFPEAK'] = (round(order.centroid, 3), 'fractional row number of profile peak')
        header['ORDERSNR'] = (round(order.snr, 3), 'sign-to-noise ratio for order')

        fluxAsciiTable(out_dir, reduced.baseName, order.orderNum, order.wavelengthScaleMeas, 
                order.objSpec, order.skySpec, order.synthesizedSkySpec, order.noiseSpec,
                order.flatSpec, order.topTrace, order.botTrace, order.avgTrace, 
                order.smoothedTrace, order.smoothedTrace - order.avgTrace)
             
        fluxFitsTable(out_dir, reduced.baseName, order.orderNum, order.wavelengthScaleMeas, 
                order.objSpec, order.skySpec, order.synthesizedSkySpec, order.noiseSpec,
                order.flatSpec, order.topTrace, order.botTrace, order.avgTrace, 
                order.smoothedTrace, order.smoothedTrace - order.avgTrace)
                
        tracePlot(out_dir, reduced.baseName, order.orderNum, order.avgTrace, 
                order.smoothedTrace, order.traceMask, order.botMeas + order.padding)
     
        traceFits(out_dir, reduced.baseName, order.orderNum, order.avgTrace)
        
        profilePlot(out_dir, reduced.baseName, order.orderNum, order.spatialProfile, 
            order.peakLocation, order.centroid, order.objWindow, order.topSkyWindow, 
            order.botSkyWindow, order.topBgMean, order.botBgMean, order.gaussianParams, order.snr)
         
        profileAsciiTable(out_dir, reduced.baseName, order.orderNum, order.spatialProfile)
         
        profileFitsTable(out_dir, reduced.baseName, order.orderNum, order.spatialProfile)
        
        profileFits(out_dir, reduced.baseName, order.orderNum, order.spatialProfile, header)

        spectrumPlot(out_dir, reduced.baseName, 'flux', order.orderNum, 
            'counts', order.objSpec, order.wavelengthScaleMeas)
         
        fitsSpectrum(out_dir, reduced.baseName, 'flux', order.orderNum, 
            'counts', order.objSpec, order.wavelengthScaleMeas, header)

        spectrumPlot(out_dir, reduced.baseName, 'sky', order.orderNum, 
            'counts', order.skySpec, order.wavelengthScaleMeas)
        
#         skyLinesPlot(out_dir, reduced.baseName, order)
#         skyLinesAsciiTable(out_dir, reduced.baseName, order)

        fitsSpectrum(out_dir, reduced.baseName, 'sky', order.orderNum, 
            'counts', order.skySpec, order.wavelengthScaleMeas, header)

        spectrumPlot(out_dir, reduced.baseName, 'noise', order.orderNum, 
            'counts', order.noiseSpec, order.wavelengthScaleMeas)
         
        fitsSpectrum(out_dir, reduced.baseName, 'noise', order.orderNum, 
            'counts', order.noiseSpec, order.wavelengthScaleMeas, header)

#         multiSpectrumPlot(out_dir, reduced.baseName, order.orderNum, 
#             'counts', order.objSpec, order.skySpec, order.noiseSpec, wavelength_scale)


        twoDimOrderFits(out_dir, reduced.baseName, order.orderNum, order.objImg, header)

        if len(order.wavelengthScaleMeas) > 0:
            twoDimOrderPlot(out_dir, reduced.baseName, 'rectified order image', 
                    reduced.getObjectName(), 'order.png', order.orderNum, order.objImg, 
                    order.wavelengthScaleMeas)

        else:
            twoDimOrderPlot(out_dir, reduced.baseName, 'rectified order image *', 
                    reduced.getObjectName(), 'order.png', order.orderNum, order.objImg, 
                    order.wavelengthScaleMeas)

    main_logger.info('n data products generated for {} = {}'.format(
            reduced.baseName, str(file_count[0])))
    return 
    
    
def tracePlot(outpath, base_name, order_num, raw, fit, mask, shift_offset):
        
    offraw = raw + shift_offset
    offfit = fit + shift_offset
    
    pl.figure("Trace Plot", figsize=(6, 5), facecolor='white')
    pl.title('trace, ' + base_name + ", order " + str(order_num), fontsize=14)
    pl.xlabel('column (pixels)')
    pl.ylabel('row (pixels)')
    
#     yrange = offraw.max() - offraw.min()
 
    x = np.arange(offraw.shape[0])
    
    pl.plot(x[mask], offraw[mask], "ko", mfc="none", ms=1.0, linewidth=1, label="derived")
    pl.plot(x, offfit, "k-", mfc="none", ms=1.0, linewidth=1, label="fit")
        
    pl.plot(x[np.logical_not(mask)], offraw[np.logical_not(mask)], 
        "ro", mfc="red", mec="red", ms=2.0, linewidth=2, label="ignored")

    
    rms = np.sqrt(np.mean(np.square(offraw - offfit)))
    pl.annotate('RMS residual = ' + "{:.3f}".format(rms), (0.3, 0.8), xycoords="figure fraction")
    
    pl.minorticks_on()
    pl.grid(True)
    pl.legend(loc='best', prop={'size': 8})

    fn = constructFileName(outpath, base_name, order_num, 'trace.png')
    pl.savefig(fn)
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
    fn = constructFileName(outpath, base_name, order_num, 'profile.tbl')
    thdulist.writeto(fn, clobber=True)  
    log_fn(fn)
    return  

def profileFits(outpath, base_name, order_num, profile, header):
    hdu = fits.PrimaryHDU(profile)
    hdulist = fits.HDUList(hdu)
    hdr = hdulist[0].header
    for k, v in header.iteritems():
        hdr[k] = v
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
    #pl.xlim(0, profile.shape[0])
    yrange = profile.max() - profile.min()
    ymin = profile.min() - (0.1 * yrange)
    ymax = profile.max() + (0.1 * yrange)
    pl.ylim(ymin, ymax)
    
    pl.minorticks_on()
    pl.grid(False)
    
    pl.plot(profile, "ko-", mfc='none', ms=3.0, linewidth=1)
    
#     pl.plot([centroid, centroid], [ymin, profile[peak]], "g-", linewidth=0.5, label='peak')
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
    
    if gaussian is None:
        width = 'unknown'
    else:
        width = '{:.1f}'.format(abs(gaussian[2]))
    if peak > (len(profile) / 2):
        pl.annotate('centroid = {:.1f} pixels\nwidth = {} pixels\nSNR = {:.1f}'.format(
                centroid, width, snr), 
                    (peak - (len(ext_range) / 2) - 20, ((ymax - ymin) * 3 / 5) + ymin))
    else:
        pl.annotate('centroid = {:.1f} pixels\nwidth = {} pixels\nSNR = {:.1f}'.format(
                centroid, width, snr), 
                    (peak + (len(ext_range) / 2) + 5, ((ymax - ymin) * 3 / 5) + ymin))
        
    # draw sky windows
 
    if (sky_range_top):

        pl.plot((sky_range_top[0], sky_range_top[-1]),
                (top_bg_mean, top_bg_mean), 
                'b', linewidth=0.5, label='sky window')  
        pl.plot((sky_range_top[0], sky_range_top[0]),
                (top_bg_mean - wvlh, top_bg_mean + wvlh), 
                'b', linewidth=1.5)
        pl.plot((sky_range_top[-1], sky_range_top[-1]),
                (top_bg_mean - wvlh, top_bg_mean + wvlh), 
                'b', linewidth=1.5)  
        
    if (sky_range_bot):
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
        pl.plot(image_lib.gaussian(range(len(profile)), gaussian[0], gaussian[1], gaussian[2]) + np.amin(profile), 
                'k--', linewidth=0.5, label='Gaussian fit')
        
    pl.legend(loc='best', prop={'size': 8})
    
    fn = constructFileName(outpath, base_name, order_num, 'profile.png')
    pl.savefig(fn)
    pl.close()
    log_fn(fn)
    return

def fluxAsciiTable(outpath, base_name, order_num, wave, flux, sky, synth_sky, error,
                     flat, trace_upper, trace_lower, trace_mean, trace_fit, fit_res):
     
    names = [   'col',         'wave',         'flux',         'sky', 
                'synth_sky',   'error',        'flat',         'trace_upper', 
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
        data = [col, wave[col], flux[col], sky[col], synth_sky[col], error[col], flat[col],
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
                fits.Column(name='error (cnts)', format='1D', array=error),
                fits.Column(name='sky (cnts)', format='1D', array=sky),
                fits.Column(name='synth_sky', format='1D', array=synth_sky),
                fits.Column(name='sig_to_noise', format='1D', array=np.arange(1024, dtype=float)),
                fits.Column(name='flat (cnts)', format='1D', array=flat),
                fits.Column(name='trace_upper (pix)', format='1D', array=trace_upper),
                fits.Column(name='trace_lower (pix)', format='1D', array=trace_lower),
                fits.Column(name='trace_mean (pix)', format='1D', array=trace_mean),
                fits.Column(name='trace_fit (pix)', format='1D', array=trace_fit),
                fits.Column(name='fit_res (pix)', format='1D', array=trace_mean - trace_fit)])
    fn = constructFileName(outpath, base_name, order_num, 'flux.tbl')   
    tbhdu.writeto(fn, clobber=True)
    log_fn(fn)
    return

def spectrumPlot(outpath, base_name, title, order_num, y_units, cont, wave):
    
    pl.figure(title, facecolor='white')
    pl.clf()
    pl.title(title + ', ' + base_name + ", order " + str(order_num), fontsize=12)
    pl.xlabel('wavelength ($\AA$)')
    pl.ylabel(title + '(' + y_units + ')')
    pl.grid(True)
    pl.plot(wave[:1004], cont[:1004], "k-", mfc="none", ms=3.0, linewidth=1)
    
    ymin, ymax = pl.ylim()
    pl.plot(wave, cont, "k-", mfc="none", ms=3.0, linewidth=1)
    pl.ylim(ymin, ymax)
    
#     axes = pl.gca()
    #axes.set_ylim(0, 300)
    
    fn = constructFileName(outpath, base_name, order_num, title + '.png')
    pl.savefig(fn)
    log_fn(fn)
    return
    
def fitsSpectrum(outpath, base_name, title, order_num, y_units, cont, wave, header):
    
    hdu = fits.PrimaryHDU(cont)
    hdulist = fits.HDUList(hdu)
    hdr = hdulist[0].header
    for k, v in header.iteritems():
        hdr[k] = v
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
        
    pl.savefig(fn)
    log_fn(fn)
    pl.close()
    return
    
#     hdu = fits.PrimaryHDU(cont)
#     hdulist = fits.HDUList(hdu)
#     hdulist[0].header['object'] = '511 Davida'
#     fits_fn = outpath + '/' + base_name + '_' + str(order) + '_' + title + '.fits'
#     logger.info('writing ' + title + ' plot for order ' + str(order) + ' to ' + fits_fn)
#     hdulist.writeto(fits_fn, clobber=True)
    
def wavelengthCalAsciiTable(outpath, base_name, order, col, source, wave_exp, wave_fit, res, peak, slope):
    
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
        
    fn = constructFileName(outpath, base_name, None, 'calids.txt')
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
    fn = constructFileName(outpath, base_name, None, 'calids.tbl')   
    thdulist.writeto(fn, clobber=True)         
    log_fn(fn)
    return
    
def twoDimOrderPlot(outpath, base_name, title, obj_name, base_filename, order_num, data, x):
    pl.figure('2d order image', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.title(title + ', ' + base_name + ", order " + str(order_num), fontsize=14)
    pl.xlabel('wavelength($\AA$)', fontsize=12)
    pl.ylabel('row (pixel)', fontsize=12)
    #pl.imshow(img, aspect='auto')
    #pl.imshow(data, vmin=0, vmax=1024, aspect='auto')
    
    pl.imshow(exposure.equalize_hist(data), origin='lower', 
                  extent=[x[0], x[-1], 0, data.shape[0]], aspect='auto')      
#     from matplotlib import colors
#     norm = colors.LogNorm(data.mean() + 0.5 * data.std(), data.max(), clip='True')
#     pl.imshow(data, norm=norm, origin='lower',
#                   extent=[x[0], x[-1], 0, data.shape[0]], aspect='auto')               
    pl.colorbar()
    pl.set_cmap('jet')
#     pl.set_cmap('Blues_r')
    fn = constructFileName(outpath, base_name, order_num, base_filename)
    pl.savefig(fn)
    log_fn(fn)
    pl.close()
    
#     np.save(fn[:fn.rfind('.')], data)
    
    return
    

def twoDimOrderFits(outpath, base_name, order_num, data, header):     
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList(hdu)
    hdr = hdulist[0].header
    for k, v in header.iteritems():
        hdr[k] = v
    fn = constructFileName(outpath, base_name, order_num, 'order.fits')
    hdulist.writeto(fn, clobber=True)
    log_fn(fn)
    return
    
        
