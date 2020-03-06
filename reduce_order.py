import logging
import numpy as np
import scipy.stats
import scipy.optimize
# import scipy.ndimage

import config
import image_lib
import nirspec_lib
import wavelength_utils
import Line
import DrpException

logger = logging.getLogger('obj')


def reduce_order(order):

    # flatten object images for this order
    __flatten(order)

    # rectify obj and flattened obj in spatial dimension
    __rectify_spatial(order)

    # trim rectified order
    __trim(order)

    # save rectified images before spectral rectify for diagnostics
    order.srNormFlatImg = order.flatOrder.rectFlatImg
    for frame in order.frames:
        order.srFfObjImg[frame] = order.ffObjImg[frame]

    # find spatial profile and peak
    __find_spatial_profile_and_peak(order)

    # characterize spatial profile by fitting to Gaussian
    __characterize_spatial_profile(order)

    # Find and smooth spectral trace, always use frame A
    try:
        order.spectralTrace = nirspec_lib.smooth_spectral_trace(
            nirspec_lib.find_spectral_trace(
                order.ffObjImg['A']), order.ffObjImg['A'].shape[0])
    except Exception as e:
        logger.warning('not rectifying order {} in spectral dimension'.format(
            order.flatOrder.orderNum))

    else:
        order.flatOrder.rectFlatImg = image_lib.rectify_spectral(
            order.flatOrder.rectFlatImg, order.spectralTrace)
        __rectify_spectral(order)
        order.spectralRectified = True

    # compute noise image
    order.noiseImg = nirspec_lib.calc_noise_img(
        order.objImg['A'], order.flatOrder.rectFlatImg, order.integrationTime)

    # extract spectra
    __extract_spectra(order)

    # calculate approximate SNR
    __calc_approximate_snr(order)

    # find and identify sky lines
    line_pairs = None  # line_pairs are (column number, accepted wavelength)
    try:
        oh_wavelengths, oh_intensities = wavelength_utils.get_oh_lines()
    except IOError as e:
        logger.critical('cannot read OH line file: ' + str(e))
        raise

    try:
        # synthesize sky spectrum and store in order object
        order.synthesizedSkySpec = wavelength_utils.synthesize_sky(
            oh_wavelengths, oh_intensities, order.flatOrder.gratingEqWaveScale)

        # identify lines and return list of (column number, accepted
        # wavelength) tuples
        line_pairs = wavelength_utils.line_id(
            order, oh_wavelengths, oh_intensities)

    except (IOError, ValueError) as e:
        logger.warning('sky line matching failed: ' + str(e))

    if line_pairs is not None:

        logger.info(str(len(line_pairs)) + ' matched sky lines found in order')

        # add line pairs to Order object as Line objects
        for line_pair in line_pairs:
            col, waveAccepted = line_pair
            peak = order.skySpec['A'][col]
            cent = image_lib.centroid(order.skySpec['A'], 1024, 5, col)
            line = Line.Line(order.baseNames['A'], order.flatOrder.orderNum,
                             waveAccepted, col, cent, peak)
            order.lines.append(line)

        if len(order.lines) >= 3:
            # do local wavelength fit
            measured = []
            accepted = []
            for line in order.lines:
                measured.append(order.flatOrder.gratingEqWaveScale[line.col])
                accepted.append(line.waveAccepted)
            (order.orderCalSlope, order.orderCalIncpt, order.orderCalCorrCoeff, p, e) = \
                scipy.stats.linregress(np.array(measured), np.array(accepted))
            order.orderCalNLines = len(order.lines)
            logger.info('per order wavelength fit: n = {}, a = {:.6f}, b = {:.6f}, r = {:.6f}'.format(
                len(order.lines), order.orderCalIncpt, order.orderCalSlope,
                order.orderCalCorrCoeff))

            for line in order.lines:
                line.orderWaveFit = order.orderCalIncpt + \
                    (order.orderCalSlope *
                     order.flatOrder.gratingEqWaveScale[line.col])
                line.orderFitRes = abs(line.orderWaveFit - line.waveAccepted)
                line.orderFitSlope = (order.orderCalSlope *
                                      (order.flatOrder.gratingEqWaveScale[1023] -
                                       order.flatOrder.gratingEqWaveScale[0])) / 1024.0
    else:
        logger.warning('no matched sky lines in order ' +
                       str(order.flatOrder.orderNum))
        order.orderCal = False

    return


def __flatten(order):
    """Flat field object image[s] but keep originals for noise calculation.
    """

    for frame in order.frames:

        order.objImg[frame] = np.array(order.objCutout[frame])

        order.ffObjImg[frame] = np.array(
            order.objCutout[frame] /
            order.flatOrder.normFlatImg)

        if frame != 'AB':
            if np.amin(order.ffObjImg[frame]) < 0:
                order.ffObjImg[frame] -= np.amin(order.ffObjImg[frame])

    order.flattened = True
    logger.info('order has been flat fielded')
    return


def __rectify_spatial(order):
    """
    """
    for frame in order.frames:
        order.objImg[frame] = image_lib.rectify_spatial(
            order.objImg[frame], order.flatOrder.smoothedSpatialTrace)
        order.ffObjImg[frame] = image_lib.rectify_spatial(
            order.ffObjImg[frame], order.flatOrder.smoothedSpatialTrace)

    order.spatialRectified = True
    logger.info('order has been rectified in the spatial dimension')

    return


def __trim(order):
    """
    """
    for frame in order.frames:
        order.objImg[frame] = \
            order.objImg[frame][order.flatOrder.botTrim:order.flatOrder.topTrim, :]
        order.ffObjImg[frame] = \
            order.ffObjImg[frame][order.flatOrder.botTrim:order.flatOrder.topTrim, :]

    return


def __rectify_spectral(order):
    """
    """
    for frame in order.frames:
        order.objImg[frame] = image_lib.rectify_spectral(
            order.objImg[frame], order.spectralTrace)
        order.ffObjImg[frame] = image_lib.rectify_spectral(
            order.ffObjImg[frame], order.spectralTrace)

    return


def __extract_spectra(order):

    if order.isPair:
        # get object extraction range for AB
        order.objWindow['AB'], _, _ = \
            image_lib.get_extraction_ranges(order.objImg['AB'].shape[0],
                                            order.peakLocation['AB'], config.params['obj_window'], None, None)

        logger.info('frame AB extraction window width = {}'.format(
            str(len(order.objWindow['AB']))))

        # extract object spectrum from AB
        order.objSpec['AB'] = np.sum(order.ffObjImg['AB'][i, :]
                                     for i in order.objWindow['AB'])

        frames = ['A', 'B']
    else:
        frames = ['A']

    for frame in frames:

        # get sky extraction ranges for A or A and B
        order.objWindow[frame], order.topSkyWindow[frame], order.botSkyWindow[frame] = \
            image_lib.get_extraction_ranges(order.objImg[frame].shape[0],
                                            order.peakLocation[frame], config.params['obj_window'],
                                            config.params['sky_window'], config.params['sky_separation'])

        logger.info('frame {} extraction window width = {}'.format(
            frame, str(len(order.objWindow[frame]))))
        logger.info('frame {} top background window width = {}'.format(
            frame, str(len(order.topSkyWindow[frame]))))
        if len(order.topSkyWindow[frame]) > 0:
            logger.info('frame {} top background window separation = {}'.format(
                frame, str(order.topSkyWindow[frame][0] - order.objWindow[frame][-1])))
        logger.info('frame {} bottom background window width = {}'.format(
            frame, str(len(order.botSkyWindow[frame]))))
        if len(order.botSkyWindow[frame]) > 0:
            logger.info('frame {} bottom background window separation = {}'.format(
                frame, str(order.objWindow[frame][0] - order.botSkyWindow[frame][-1])))

        # extract object, sky and noise spectra for A and B and flat spectrum
        order.objSpec[frame], order.flatSpec, order.skySpec[frame], order.noiseSpec[frame], \
            order.topBgMean[frame], order.botBgMean[frame] = image_lib.extract_spectra(
            order.ffObjImg[frame], order.flatOrder.rectFlatImg, order.noiseImg,
            order.objWindow[frame], order.topSkyWindow[frame],
            order.botSkyWindow[frame])

    if order.isPair:
        order.noiseSpec['AB'] = order.noiseSpec['A'] + order.noiseSpec['B']

    return


def __calc_approximate_snr(order):

    if order.isPair:
        frames = ['A', 'B']
    else:
        frames = ['A']

    for frame in frames:

        bg = 0.0

        if order.topBgMean[frame] is not None:
            bg += order.topBgMean[frame]
        if order.botBgMean[frame] is not None:
            bg += order.botBgMean[frame]
        if order.topBgMean[frame] is not None and order.botBgMean[frame] is not None:
            bg /= 2

        order.snr[frame] = np.absolute(
            np.mean(order.ffObjImg[frame]
                    [order.peakLocation[frame]: order.peakLocation[frame] + 1, :]) / bg)

        logger.info(
            'frame {} signal-to-noise ratio = {:.1f}'.format(frame, order.snr[frame]))

    if order.isPair:
        order.snr['AB'] = 0.0
        if order.snr['A'] is not None:
            order.snr['AB'] = order.snr['A']
        if order.snr['B'] is not None:
            order.snr['AB'] = order.snr['AB'] + order.snr['B']
        if order.snr['A'] is not None and order.snr['B'] is not None:
            order.snr['AB'] = order.snr['AB'] / 2.0

    return


def __characterize_spatial_profile(order):

    for frame in order.frames:
        try:
            for w in range(10, 30, 10):
                logger.debug('gaussian window width = {}'.format(2 * w))
                x0 = max(0, order.peakLocation[frame] - w)
                x1 = min(
                    len(order.spatialProfile[frame]) - 1, order.peakLocation[frame] + w)
                x = list(range(x1 - x0))
                order.gaussianParams[frame], pcov = scipy.optimize.curve_fit(
                    image_lib.gaussian, x, order.spatialProfile[frame][x0:x1] -
                    np.amin(order.spatialProfile[frame][x0:x1]))
                order.gaussianParams[frame][1] += x0
                if order.gaussianParams[frame][2] > 1.0:
                    break
        except Exception as e:
            logger.warning(
                'cannot fit frame {} spatial profile to Gaussian'.format(frame))
            order.gaussianParams[frame] = None
        else:
            logger.info('frame {} spatial peak width = {:.1f} pixels'.format(
                frame, abs(order.gaussianParams[frame][2])))

    return


def __find_spatial_profile_and_peak(order):
    """
    """

    MARGIN = 5

    for frame in order.frames:

        # find spatial profile(s)
        order.spatialProfile[frame] = order.ffObjImg[frame].mean(axis=1)
        if len(order.spatialProfile[frame]) < (2 * MARGIN) + 2:
            raise DrpException.DrpException(
                'cannot find spatial profile for frame {} order {}'.format(
                    frame, order.flatOrder.orderNum))

        # find peak locations
        order.peakLocation[frame] = np.argmax(
            order.spatialProfile[frame][MARGIN:-MARGIN]) + MARGIN
        logger.info('frame {} spatial profile peak intensity row {:d}'.format(
            frame, order.peakLocation[frame]))

        # fit peak to Gaussian, save Gaussian parameters and real centroid
        # location
        p0 = order.peakLocation[frame] - (config.params['obj_window'] // 2)
        p1 = order.peakLocation[frame] + (config.params['obj_window'] // 2)
        order.centroid[frame] = (scipy.ndimage.measurements.center_of_mass(
            order.spatialProfile[frame][p0:p1]))[0] + p0
        logger.info('frame {} spatial profile peak centroid row {:.1f}'.format(
            frame, float(order.centroid[frame])))

    return


def __calculate_SNR(order):

    bg = 0.0

    if order.topBgMean is not None:
        bg += order.topBgMean
    if order.botBgMean is not None:
        bg += order.botBgMean
    if order.topBgMean is not None and order.botBgMean is not None:
        bg /= 2
    order.snr = np.absolute(np.mean(
        order.ffObjImg['A'][order.peakLocation['A']:order.peakLocation['A'] + 1, :]) / bg)
