import imp

MODULES = [
    'os',
    'numpy',
    'math',
    'subprocess',
    'fnmatch',
    'logging',
    'pylab',
    'errno',
    'datetime',
    'warnings',
    'astropy',
    'scipy',
    'argparse',
    'statsmodels',
    'PIL']

missingModules = []


def is_missing():
    for m in MODULES:
        try:
            imp.find_module(m)
        except ImportError:
            missingModules.append(m)

    return missingModules
