


class WaveCalTable:
    
    def __init__(self, rds):
        self.rds = rds
        self.names = []
        self.units = []
        self.formats = []
        self.widths = []
        
        self.defCols()
        self.collateData()
        
        
    def genAsciiTable(self, fn):
        pass
        
    def defCols(self):
        self.defCol('frame',    '',     '',     10)
        self.defCol('order',    '',     'd',    6)
        self.defCol('source',   '',     '',     '')
        self.defCol('col',      '',     '',     '')
        self.defCol('peak',     '',     '',     '')
        self.defCol('wave_ge',  '',     '',     '')
        self.defCol('wave_exp', '',     '',     '')
        self.defCol('used_o',   '',     '',     '')
        self.defCol('wave_o',   '',     '',     '')
        self.defCol('res_o',    '',     '',     '')
        self.defCol('disp_o',   '',     '',     '')
        self.defCol('used_f',   '',     '',     '')
        self.defCol('save_f',   '',     '',     '')
        self.defCol('res_f',    '',     '',     '')
        self.defCol('disp_f',   '',     '',     '')
        self.defCol('used_mf',  '',     '',     '')
        self.defCol('wave_mf',  '',     '',     '')
        self.defCol('res_mf',   '',     '',     '')
        self.defCol('disp_mf',  '',     '',     '')

    def defCol(self, name, units, format, width):
        self.names.append(name)
        self.units.append(units)
        self.formats.append(format)
        return
    
    def collateData(self):
        
        
        
    def genFitsTable(self, fn):