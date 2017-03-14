
class Line:
    """
    """
    
    def __init__(self, frame, order, wAccepted, col, centroid, peak):
        
        # line identification
        self.frame = frame              # frame name or KOAID from which line extracted
        self.col = col                  # integer column number of line, pixels
        self.centroid = centroid        # centroid column of line, fractional pixels
        self.order = order              # order number in which this line was found
        self.waveAccepted = wAccepted   # accepted wavelength from atlas, Angstroms
        self.peak = peak                # peak intensity, counts


        # per-order fit results
        self.orderWaveFit = 0.0         # wavelength of line based on per order fit (Angstroms)
        self.orderFitRes = 0.0          # fit residual (Angstroms)
        self.orderFitSlope = 0.0        # slope of per-order WL equation at this wavelength
                                        # (Angstroms/pixel)
        self.orderFitOutlier = False    # true if line discarded from per-order fit as outlier
 
 
        # per-frame fit results
        self.frameWaveFit = 0.0         # wavelength of line based on per-frame fit (Angstroms)
        self.frameFitRes = 0.0          # fit residual (Angstroms)
        self.frameFitSlope = 0.0        # slope of per-frame WL equation at this WL
                                        # (Angstroms/pixel)
        self.frameFitOutlier = True     # true if line discarded from per-frame fit as outlier
        

        # multi-frame fit results
        self.multiFrames = []           # names or KOAIDs of other frames in multi-frame fit
        self.mfWaveFit = 0.0            # wavelength of line based on multi-frame fit (Angstroms)
        self.mfFitRes = 0.0             # fit residual (Angstroms)
        self.mfFitSlope = 0.0           # slope of per-frame WL equation at this wavelength 
                                        # (Angstroms/pixel)
        self.multiFrameOutlier = False  # true if line discarded from fit as outlier
        