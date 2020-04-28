import numpy as np
import scipy.ndimage


def trace_edge(data, start, searchWidth, bgWidth, jumpThresh):

    # initialize trace array
    trace = np.zeros(data.shape[1])

    # nJumps is the number of times successive centroids differed by more than
    # threshold
    nJumps = 0

    # first centroid assumed to be at start
    trace[0] = start

    # find centroids for the rest of the columns in data

    for i in range(1, data.shape[1]):

        # define search window
        ymin = int(trace[i - 1] - searchWidth)
        ymax = int(trace[i - 1] + searchWidth)

        # clip search window at top and bottom of column
        if abs(ymax) > data.shape[0]:
            ymax = int(data.shape[0])

        if ymin < 1:  # don't let it trace the bottom of the detector
            ymin = 1
        if ymax <= 0:
            # how can this happen?
            ymax = int(trace[i] + searchWidth) + 1

        if bgWidth <= 0:
            bgMean = 0.0

        else:
            # If bgWidth > 0 then we will subtract average of pixel values at two locations
            # from each pixel value in search window.  Two locations are previous centroid
            # plus and minus search width.

            bgMin = trace[i - 1] - bgWidth
            if bgMin < 0:
                bgMin = 0
            bgMax = trace[i - 1] + bgWidth
            if bgMax > data.shape[0]:
                bgMax = data.shape[0] - 1

            try:
                bgMean = (data[bgMin, i] + data[bgMax, i]) / 2.0
            except BaseException:
                bgMean = 0.0

        trace[i] = scipy.ndimage.measurements.center_of_mass(
            data[int(ymin):int(ymax) + 1, i] - bgMean)[0] + ymin

#         import pylab as pl
#         x0 = max(0, int(trace[i - 1]) - 50)
#         x1 = min(1023, int(trace[i-1]) + 50)
#         pl.figure()
#         pl.cla()
#         pl.plot(np.arange(x0, x1), data[x0:x1, i])
#         pl.plot([trace[i], trace[i]], pl.ylim())


#         pl.plot(data[x0:x1, i], 'ro')
#         pl.plot(data[int(ymin):int(ymax) + 1, i] - bgMean, 'go')
#         pl.plot([trace[i], trace[i]], [0, pl.ylim()[0]], 'g-')
#         print(trace[max(0, i-10):i])


#         pl.show()

        if trace[i] is np.inf or trace[i] is -np.inf:
            # went off array
            print('went off array')
            return None, None

        # centroid jumped more than traceDelta
        if np.abs(trace[i] - trace[i - 1]) > jumpThresh:
            nJumps += 1
            if i > 4:
                # jump is past beginning, use past three centroids
                trace[i] = trace[i - 3:i - 1].mean()
            elif i > 1:
                # average as many traces as we have gone through
                trace[i] = trace[i - 2:i - 1].mean()
            else:
                # use the first one found
                trace[i] = trace[i - 1]

    return trace, nJumps
