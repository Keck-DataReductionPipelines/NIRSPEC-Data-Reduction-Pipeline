import os

modules = ['numpy', 'math', 'subprocess', 'fnmatch', 'logging', 'pylab', 'errno', 'datetime', 'warnings', 'astropy', 'scipy', 'argparse', 'statsmodels']

for m in modules:
    if m not in os.sys.modules:
        print 'Module %s not installed' % m
        
