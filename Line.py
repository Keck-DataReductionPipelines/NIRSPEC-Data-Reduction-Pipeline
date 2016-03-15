
class Line:
    
    def __init__(self):
        
        self.col = -1                   # integer column number of line, pixels
        self.centroid = None            # centroid on wavelength axis, fractional pixels
        self.order = 0                  # order number in which this line was found
        self.acceptedWavelength = 0.0   # accepted wavelength from atlas, Angstroms
        self.peak = 0.0                 # peak intensity, counts

        self.localFitWavelength = 0.0
        self.localFitResidual = 0.0
        self.localFitSlope = 0.0
        
        self.globalFitWavelength = 0.0
        self.globalFitResidual = 0.0
        self.globalFitSlope = 0.0
        self.usedInGlobalFit = False
        
