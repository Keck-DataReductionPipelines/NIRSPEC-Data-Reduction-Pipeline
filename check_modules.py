import imp

modules = ['os', 'numpy', 'math', 'subprocess', 'fnmatch', 'logging', 'pylab', 'errno', 'datetime', 'warnings', 'astropy', 'scipy', 'argparse', 'statsmodels', 'PIL']

for m in modules:
    try:
        imp.find_module(m)
    except ImportError:
        print 'Module %s not installed' % m
