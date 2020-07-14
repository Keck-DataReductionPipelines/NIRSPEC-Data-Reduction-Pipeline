import subprocess as sp
import sys
import importlib


MODULES = [
    'os',
    'numpy',
    'math',
    'subprocess',
    'fnmatch',
    'logging',
    'matplotlib',
    'errno',
    'datetime',
    'warnings',
    'astropy',
    'scipy',
    'argparse',
    'statsmodels',
    'PIL',
    'astroscrappy',
    'skimage']

missingModules = []
installedModules = []
failedToInstall = []


def is_missing():
    for module in MODULES:
        try:
            importlib.import_module(module)
        except ImportError:
            missingModules.append(module)

    return missingModules


def install_packages(packages):
    print(('Missing modules: ' + ', '.join(packages)))

    if 'PIL' in packages:
        idx = packages.index('PIL')
        packages[idx] = 'Pillow'

    if 'skimage' in packages:
        idx = packages.index('skimage')
        packages[idx] = 'scikit-image'

    for package in packages:
        try:
            if sys.version_info[0]  == 2:
                sp.check_call([sys.executable, '-m', 'pip', 'install', '-q', package], stderr=sp.PIPE)
            else:
                sp.check_call([sys.executable, '-m', 'pip', 'install', '-q', package], stderr=sp.DEVNULL)

            installedModules.append(package)

        except sp.CalledProcessError:
            failedToInstall.append(package)

        if package == 'Pillow':
            package = 'PIL'

        if package == 'scikit-image':
            package = 'skimage'

        #globals()[package] = importlib.import_module(package)

    return installedModules, failedToInstall
