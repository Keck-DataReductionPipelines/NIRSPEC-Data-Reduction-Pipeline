
class Line:
    
    def __init__(self):
        
        self.col = -1 
        self.order = 0
        self.acceptedWavelength = 0.0
        self.peak = 0.0

        self.localFitWavelength = 0.0
        self.localFitResidual = 0.0
        self.localFitSlope = 0.0
        
        self.globalFitWavelength = 0.0
        self.globalFitResidual = 0.0
        self.globalFitSlope = 0.0
        self.usedInGlobalFit = False
        
