# from __builtin__ import False


class NirspecConfig:

    def __init__(self, header):
        self.header = header

    def isTheSame(self, header):
        for kwd in ['disppos', 'echlpos', 'filname', 'slitname']:
            if self.header[kwd] != header[kwd]:
                return False
        return True

    def toString(self):
        return 'disppos={}, echlpos={}, filname={}, slitname={}'.format(
            self.header['disppos'], self.header['echlpos'],
            self.header['filname'], self.header['slitname'])
