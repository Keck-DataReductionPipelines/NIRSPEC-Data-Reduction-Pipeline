import logging
import os
import errno
import numpy as np
import pylab as pl
from astropy.io import fits
from astropy.table import Table, Column
from astropy.io import ascii
import Image

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
UNDERLINE = False

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
                ('localcalids.txt', 'ascii/cal'         ),
                ('order.fits',      'fits/order'        ),
                ('order.png',       'previews/order'    ),
                ('sky.fits',        'fits/sky'          ),
                ('sky.png',         'previews/sky'      ),
                ('noise.fits',      'fits/noise'        ),
                ('noise.png',       'previews/noise'    ),
                ('trace.fits',      'fits/trace'        ),
                ('trace.png',       'previews/trace'    ),
                ('spectra.png',     'previews/spectra'  ),
                ('skylines.png',    'previews/skylines' ),
                ('skylines.txt',    'ascii/skylines'    )
                ])

def constructFileName(outpath, base_name, order, fn_suffix):
    fn = outpath + '/' + subdirs[fn_suffix] + '/' + base_name + '_' + fn_suffix
    if order is None:
        return fn
    else:
        return fn[:fn.rfind('_') + 1] + str(order) + fn[fn.rfind('_'):]
   
def log_fn(fn):  
        obj_logger.info('saving {}'.format(fn))
        file_count[0] += 1 
        return

def gen(reduced, out_dir):
            
    # make sub directories
    for k, v in subdirs.iteritems():
        try:
            os.makedirs(out_dir + '/' + v)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                main_logger.critical('failed to create output directory ' + v)
                obj_logger.critical('failed to create output directory ' + v)
                raise IOError('failed to create output directory ' + v)
    
    # per frame data products
    
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
                col.append(line.col)
                source.append('sky')
                wave_exp.append(line.acceptedWavelength)
                wave_fit.append(line.globalFitWavelength)
                res.append(abs(line.globalFitWavelength - line.acceptedWavelength))
                peak.append(line.peak)
                slope.append(line.globalFitSlope)
                
    wavelengthCalAsciiTable(
            out_dir, reduced.baseName, order_num, col, source, wave_exp, wave_fit, res, peak, slope)
    
    order_num = []
    col = []
    source = []
    wave_exp = []
    wave_fit = []
    res = []
    peak = []
    slope = []
    
    for order in reduced.orders:
        if order.perOrderSlope > 0.95 and order.perOrderSlope < 1.05:
            for line in order.lines:
                order_num.append(order.orderNum)
                col.append(line.col)
                source.append('sky')
                wave_exp.append(line.acceptedWavelength)
                wave_fit.append(line.localFitWavelength)
                res.append(line.localFitResidual)
                peak.append(line.peak)
                slope.append(line.localFitSlope)
                
    perOrderWavelengthCalAsciiTable(
            out_dir, reduced.baseName, order_num, col, source, wave_exp, wave_fit, res, peak, slope)
            
    # per order data products
    for order in reduced.orders:
        
        if reduced.coeffs is None:
            wavelength_scale = order.wavelengthScaleMeas
        else:
            wavelength_scale = order.wavelengthScaleCalc
            
        fluxAsciiTable(out_dir, reduced.baseName, order.orderNum, wavelength_scale, 
                order.objSpec, order.skySpec, order.synthesizedSkySpec, order.noiseSpec,
                np.arange(1024), order.topTrace, order.botTrace, order.avgTrace, 
                order.smoothedTrace, order.smoothedTrace - order.avgTrace)
             
        tracePlot(out_dir, reduced.baseName, order.orderNum, order.avgTrace, 
                order.smoothedTrace, order.traceMask, order.botMeas + order.padding)
     
        profilePlot(out_dir, reduced.baseName, order.orderNum, order.spatialProfile, 
            order.peakLocation, order.centroid, order.objWindow, order.topSkyWindow, 
            order.botSkyWindow)
         
        profileAsciiTable(out_dir, reduced.baseName, order.orderNum, order.spatialProfile)
#         
#         profileFitsTable(out_dir, reduced.baseName, order.orderNum, order.spatialProfile)
#         file_count += 1

        spectrumPlot(out_dir, reduced.baseName, 'flux', order.orderNum, 
            'counts', order.objSpec, order.wavelengthScaleMeas)
#         
#         fitsSpectrum(out_dir, reduced.baseName, 'flux', order.orderNum, 
#             'counts', order.objSpec, order.wavelengthScaleMeas)
#         file_count += 1

#         
        spectrumPlot(out_dir, reduced.baseName, 'sky', order.orderNum, 
            'counts', order.skySpec, order.wavelengthScaleMeas)
        
        skyLinesPlot(out_dir, reduced.baseName, order)
        skyLinesAsciiTable(out_dir, reduced.baseName, order)

#         
#         fitsSpectrum(out_dir, reduced.baseName, 'sky', order.orderNum, 
#             'counts', order.skySpec, order.wavelengthScaleMeas)
#         file_count += 1

#         
        spectrumPlot(out_dir, reduced.baseName, 'noise', order.orderNum, 
            'counts', order.noiseSpec, order.wavelengthScaleMeas)

#         
#         fitsSpectrum(out_dir, reduced.baseName, 'noise', order.orderNum, 
#             'counts', order.noiseSpec, order.wavelengthScaleMeas)

         
        multiSpectrumPlot(out_dir, reduced.baseName, order.orderNum, 
            'counts', order.objSpec, order.skySpec, order.noiseSpec, wavelength_scale)


        twoDimOrderFits(out_dir, reduced.baseName, order.orderNum, order.objImg)

        if len(order.wavelengthScaleMeas) > 0:
            twoDimOrderPlot(out_dir, reduced.baseName, 'rectified order image', 
                    reduced.getObjectName(), order.orderNum, order.objImg, 
                    order.wavelengthScaleMeas)

        else:
            twoDimOrderPlot(out_dir, reduced.baseName, 'rectified order image *', 
                    reduced.getObjectName(), order.orderNum, order.objImg, 
                    order.wavelengthScaleCalc)

    main_logger.info(str(file_count[0]) + ' product files saved for ' + reduced.baseName)  
    return 
    
    
def tracePlot(outpath, base_name, order_num, raw, fit, mask, shift_offset):
        
    offraw = raw + shift_offset
    offfit = fit + shift_offset
    
    pl.figure("Trace Plot", figsize=(6, 5), facecolor='white')
    pl.title('trace, ' + base_name + ", order " + str(order_num), fontsize=14)
    pl.xlabel('column (pixels)')
    pl.ylabel('row (pixels)')
    
    yrange = offraw.max() - offraw.min()
#     ymin = offraw.min() - (0.1 * yrange)
#     ymax = offraw.max() + (0.1 * yrange)
 
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

def profileAsciiTable(outpath, base_name, order_num, profile):
    
    names = ['row', 'mean_flux'] 
    units = ['pixels', 'counts']
    formats = ['.0f', '.0f']
    widths = [6, 9] 
    
    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
        
    if UNDERLINE:
        line = []
        for i, unit in enumerate(units):
            line.append('{:->{w}}'.format('', w=widths[i]))
        buff.append('--'.join(line))
    
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
    
def profilePlot(outpath, base_name, order_num, profile, peak, centroid,
            ext_range, sky_range_top, sky_range_bot):

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
    pl.grid(True)
    
    pl.plot(profile, "ko-", mfc='none', ms=3.0, linewidth=1)
    
    pl.plot([centroid, centroid], [ymin, profile[peak]], "g-", linewidth=0.5, label='peak')

    
    # draw extraction window
    wvlh = 0.01 * yrange;
    ewh = 0.05 * yrange
    
    pl.plot((peak + ext_range[0], peak + ext_range[-1]), (ymax - ewh, ymax - ewh), 
            'r', linewidth=0.5, label='extraction window')
    pl.plot((peak + ext_range[0], peak + ext_range[0]),
            (ymax - ewh - wvlh, ymax - ewh + wvlh), 
            'r', linewidth=1.5)
    pl.plot((peak + ext_range[-1], peak + ext_range[-1]),
            (ymax - ewh - wvlh, ymax - ewh + wvlh), 
            'r', linewidth=1.5)  
    
    # indicate centroid location
    pl.annotate('centroid = ' + "{:.1f} pixels".format(centroid), 
                (peak + ext_range[-1] + 5, (ymax - ymin) * 4 / 5))
        
    # draw sky windows
    swh = 0.2 * yrange
    if (sky_range_top):
        pl.plot((peak + sky_range_top[0], peak + sky_range_top[-1]),
                (profile.min() + swh, profile.min() + swh), 
                'b', linewidth=0.5, label='sky window')  
        pl.plot((peak + sky_range_top[0], peak + sky_range_top[0]),
                (profile.min() + swh - wvlh, profile.min() + swh + wvlh), 
                'b', linewidth=1.5)
        pl.plot((peak + sky_range_top[-1], peak + sky_range_top[-1]),
                (profile.min() + swh - wvlh, profile.min() + swh + wvlh), 
                'b', linewidth=1.5)  
        
    if (sky_range_bot):
        pl.plot((peak + sky_range_bot[0], peak + sky_range_bot[-1]),
                (profile.min() + swh, profile.min() + swh), 
                'b', linewidth=0.5)   
        pl.plot((peak + sky_range_bot[0], peak + sky_range_bot[0]),
                (profile.min() + swh - wvlh, profile.min() + swh + wvlh), 
                'b', linewidth=1.5)
        pl.plot((peak + sky_range_bot[-1], peak + sky_range_bot[-1]),
                (profile.min() + swh - wvlh, profile.min() + swh + wvlh), 
                'b', linewidth=1.5)
        
    pl.legend(loc='best', prop={'size': 8})
    
    fn = constructFileName(outpath, base_name, order_num, 'profile.png')
    pl.savefig(fn)
    pl.close()
    log_fn(fn)
    return

def fluxAsciiTable(outpath, base_name, order_num, wave, flux, sky, synth_sky, sigma,
                     flat, trace_upper, trace_lower, trace_mean, trace_fit, fit_res):
     
    names = [   'col',         'wave',         'flux',         'sky', 
                'synth_sky',   'sigma',        'flat',         'trace_upper', 
                'trace_lower', 'trace_mean',   'trace_fit',    'fit_res']
    units = [   'pixels',      'Angstroms',    'counts',       'counts', 
                '',            '',             '',             'pixel', 
                'pixels',      'pixels',       'pixels',       'pixels']
    formats = [ 'd',            '.6e',          '.3e',          '.3e',
                '.3e',          '.3e',          '.3e',          '.3e',
                '.3e',          '.3e',          '.3e',          '.3e']
    nominal_width = 10
    widths = []
    
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
    buff.append('| {} |'.format('| '.join(line)))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
        
    if UNDERLINE:
        line = []
        for i, unit in enumerate(units):
            line.append('| {} |'.format('{:->{w}}'.format('', w=widths[i])))
        buff.append('--'.join(line))
    
    for col in range(wave.shape[0]):
        data = [col, wave[col], flux[col], sky[col], synth_sky[col], sigma[col], flat[col],
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


    
def fluxFitsTable(outpath, base_name, order_num, wave, flux, sky, synth_sky,
                     trace_upper, trace_lower, trace_mean, trace_fit):
    # create binary FITS table
    prihdr = fits.Header()
    prihdr['COMMENT'] = "flux table"
    prihdu = fits.PrimaryHDU(header=prihdr)
    tbhdu = fits.new_table(
            fits.ColDefs([fits.Column(name='col', format='1I', array=np.arange(1024, dtype=int)),
                         fits.Column(name='wave (A)', format='1D', array=wave),
                         fits.Column(name='flux (cnts)', format='1D', array=flux),
                         fits.Column(name='error (cnts)', format='1D', array=np.arange(1024, dtype=float)),
                         fits.Column(name='sky (cnts)', format='1D', array=sky),
                         fits.Column(name='synth_sky', format='1D', array=synth_sky),
                         fits.Column(name='sig_to_noise', format='1D', array=np.arange(1024, dtype=float)),
                         fits.Column(name='flat (cnts)', format='1D', array=np.arange(1024, dtype=float)),
                         fits.Column(name='trace_upper (pix)', format='1D', array=trace_upper),
                         fits.Column(name='trace_lower (pix)', format='1D', array=trace_lower),
                         fits.Column(name='trace_mean (pix)', format='1D', array=trace_mean),
                         fits.Column(name='trace_fit (pix)', format='1D', array=trace_fit),
                         fits.Column(name='fit_res (pix)', format='1D', array=trace_mean - trace_fit)]))
    thdulist = fits.HDUList([prihdu, tbhdu]) 
    fn = constructFileName(outpath, base_name, order_num, 'flux.tbl')   
    thdulist.writeto(fn, clobber=True)         
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
    
    axes = pl.gca()
    #axes.set_ylim(0, 300)
    
    fn = constructFileName(outpath, base_name, order_num, title + '.png')
    pl.savefig(fn)
    log_fn(fn)
    return
    
def fitsSpectrum(outpath, base_name, title, order_num, y_units, cont, wave):
    
    hdu = fits.PrimaryHDU(cont)
    hdulist = fits.HDUList(hdu)
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
    pl.plot(wave[:1004], noise[:1004], "r-", mfc="none", ms=3.0, linewidth=1, label='noise (1 sigma')

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
    formats = [ 'd', '', 'd', '.6e', '.6e', '.3f', 'd', '.3e']
    nominal_width = 10
    widths = []
    
    for name in names:
        widths.append(max(len(name), nominal_width))
            
    for i in range(len(names)):
        widths[i] = max(len(units[i]), widths[i])
       
    widths[0] = 6 
    widths[1] = 6
    widths[2] = 6
    widths[3] = 13
    widths[4] = 13
    widths[5] = 9
    widths[6] = 6

    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
      
    if UNDERLINE:  
        line = []
        for i, unit in enumerate(units):
            line.append('{:->{w}}'.format('', w=widths[i]))
        buff.append('| {} |'.format('--'.join(line)))
    
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

def perOrderWavelengthCalAsciiTable(outpath, base_name, order, col, source, wave_exp, wave_fit, res, peak, slope):
    
    names = ['order', 'source', 'col', 'wave_exp', 'wave_fit', 'res', 'peak', 'disp'] 
    units = ['', '', 'pixels', 'Angstroms', 'Angstroms',  'Angstroms', 'counts', 'Angstroms/pixel']
    formats = [ 'd', '', 'd', '.6e', '.6e', '.3f', 'd', '.3e']
    nominal_width = 10
    widths = []
    
    for name in names:
        widths.append(max(len(name), nominal_width))
            
    for i in range(len(names)):
        widths[i] = max(len(units[i]), widths[i])
       
    widths[0] = 6 
    widths[1] = 6
    widths[2] = 6
    widths[3] = 13
    widths[4] = 13
    widths[5] = 9
    widths[6] = 6

    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
      
    if UNDERLINE:  
        line = []
        for i, unit in enumerate(units):
            line.append('{:->{w}}'.format('', w=widths[i]))
        buff.append('| {} |'.format('--'.join(line)))
    
    for i in range(len(order)):
        data = [order[i], source[i], col[i], wave_exp[i], wave_fit[i], res[i], (int)(peak[i]), slope[i]]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
        
    #print('\n'.join(buff))
        
    fn = constructFileName(outpath, base_name, None, 'localcalids.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
    log_fn(fn)
    return
    
def wavelengthCalFitsTable(outpath, base_name, order, col, source, wave_exp, wave_fit, peak, slope):

    # create binary FITS table
    prihdr = fits.Header()
    prihdr['COMMENT'] = "wavelength calibration table"
    prihdu = fits.PrimaryHDU(header=prihdr)
    tbhdu = fits.new_table(
            fits.ColDefs([fits.Column(name='order', format='B', array=order),
                         fits.Column(name='source', format='16A', array=source),
                         fits.Column(name='col (pix)', format='I', array=col),
                         fits.Column(name='wave_exp (A)', format='D', array=wave_exp),
                         fits.Column(name='wave_fit (A)', format='D', array=wave_fit),
                         fits.Column(name='peak (cnts)', format='D', array=peak),
                         fits.Column(name='disp (A/pix)', format='D', array=slope)]))
    thdulist = fits.HDUList([prihdu, tbhdu])   
    fn = constructFileName(outpath, base_name, order, 'calid_tbl.fits')
    thdulist.writeto(fn, clobber=True)                
    log_fn(fn)
    return 

def skyLinesPlot(outpath, base_name, order):
    pl.figure('sky lines', facecolor='white', figsize=(8,5))
    pl.cla()
    pl.title("sky lines" + ', ' + base_name + ", order " + str(order.orderNum), fontsize=14)
    pl.xlabel('x (pixels)', fontsize=12)
    #pl.ylabel('row (pixel)', fontsize=12)
    
    synmax = np.amax(order.synthesizedSkySpec);
    skymax = np.amax(order.skySpec);
    

    if synmax > skymax:
        syn = order.synthesizedSkySpec
        sky = order.skySpec * (synmax / skymax)
    else:
        syn = order.synthesizedSkySpec * (skymax / synmax)
        sky = order.skySpec

    pl.plot((syn * 0.4) + (max(synmax, skymax) / 2), 'g-', linewidth=1, label='synthesized sky')
    pl.plot(sky * 0.4, 'b-', linewidth=1, label='sky')

    ymin, ymax = pl.ylim()
    ymax = max(synmax, skymax) * 1.1
    ymin = -200
    pl.ylim(ymin, ymax)
    pl.xlim(0, 1024)

    pl.legend(loc='best', prop={'size': 8})
    
    pl.annotate('shift = ' + "{:.3f}".format(order.wavelengthShift), 
                (0.3, 0.8), xycoords="figure fraction")
    
    dy = 0
    for line in order.lines:
        pl.plot([line.col, line.col], [ymin, ymax], "k--", linewidth=0.5)
        pl.annotate(str(line.acceptedWavelength), (line.col, ((ymax - ymin) / 2) + dy), size=8)
        pl.annotate(str(line.col) + ', ' + '{:.3f}'.format(order.wavelengthScaleCalc[line.col]), (line.col, ((ymax - ymin) / 3) - dy), size=8)
        dy += 300
        if dy > 1500:
            dy = 0

    
    #pl.minorticks_on()
    pl.grid(True)

    fn = constructFileName(outpath, base_name, order.orderNum, 'skylines.png')
        
    pl.savefig(fn)
    log_fn(fn)
    pl.close()
    return
    

def skyLinesAsciiTable(outpath, base_name, order):
    
    names = ['order', 'col', 'calc', 'accepted'] 
    units = [' ', 'pixels', '$\AA$', '$\AA$']
    formats = ['.0f', '.0f', '.1f', '.1f']
    widths = [6, 6, 6, 9] 
    
    buff = []
    
    line = []
    for i, name in enumerate(names):
        line.append('{:>{w}}'.format(name, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
    
    line = []
    for i, unit in enumerate(units):
        line.append('{:>{w}}'.format(unit, w=widths[i]))
    buff.append('| {} |'.format('| '.join(line)))
        
    if UNDERLINE:
        line = []
        for i, unit in enumerate(units):
            line.append('{:->{w}}'.format('', w=widths[i]))
        buff.append('--'.join(line))
    
    for l in order.lines:
        data = [order.orderNum, l.col, order.wavelengthScaleCalc[l.col], l.acceptedWavelength]
        line = []
        for i, val in enumerate(data):
            line.append('{:>{w}{f}}'.format(val, w=widths[i], f=formats[i]))
        buff.append('  {}  '.format('  '.join(line)))
                
    fn = constructFileName(outpath, base_name, order.orderNum, 'skylines.txt')
    fptr = open(fn, 'w')
    fptr.write('\n'.join(buff))
    fptr.close()
    log_fn(fn)
 
    return
    
def twoDimOrderPlot(outpath, base_name, title, obj_name, order_num, data, x):
    pl.figure('2d order image', facecolor='white', figsize=(8, 5))
    pl.cla()
    pl.title(title + ', ' + base_name + ", order " + str(order_num), fontsize=14)
    pl.xlabel('wavelength($\AA$)', fontsize=12)
    pl.ylabel('row (pixel)', fontsize=12)
    #pl.imshow(img, aspect='auto')
    #pl.imshow(data, vmin=0, vmax=1024, aspect='auto')
    
    pl.imshow(data, origin='lower', vmin=0, vmax=1024, 
                  extent=[x[0], x[-1], 0, data.shape[0]], aspect='auto')               
    pl.colorbar()
    fn = constructFileName(outpath, base_name, order_num, 'order.png')
    pl.savefig(fn)
    log_fn(fn)
    pl.close()
    
    np.save(fn[:fn.rfind('.')], data)
    
    return
    

def twoDimOrderFits(outpath, base_name, order_num, data):     
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList(hdu)
#     hdulist[0].header['object'] = '511 Davida'
    fn = constructFileName(outpath, base_name, order_num, 'order.fits')
    hdulist.writeto(fn, clobber=True)
    log_fn(fn)
    return
    
        
def quickImagePlot(img, title, x_label, y_label):
    #return;
    pl.figure('quick image plot', facecolor='white', figsize=(8, 6))
    pl.cla()
    pl.title(title, fontsize=14)
    pl.xlabel(x_label, fontsize=12)
    pl.ylabel(y_label, fontsize=12)
    #pl.imshow(img, aspect='auto')
    pl.imshow(img, vmin=0, vmax=256, aspect='auto', cmap="gray")
#     pl.colorbar()
    #pl.set_cmap('spectral')
    pl.show()
    
    
def quickPlot(data, title):
    pl.figure('quick')
    pl.cla()
    pl.title(title)
    pl.plot(data, "ro")
    pl.show()
    