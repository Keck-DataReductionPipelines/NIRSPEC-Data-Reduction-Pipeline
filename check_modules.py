import subprocess
import sys
import importlib
import pip


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


def is_missing():
    for module in MODULES:
        try:
            importlib.import_module(module)
        except ImportError:
            missingModules.append(module)

    return missingModules


def install_packages(packages):
    print(('Missing modules: ' + ', '.join(packages)))
    print('XXXXXXX' * 9)

    if 'PIL' in packages:
        idx = packages.index('PIL')
        packages[idx] = 'Pillow'

    if 'skimage' in packages:
        idx = packages.index('skimage')
        packages[idx] = 'scikit-image'

    for package in packages:
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', package])
            print('Installed ' + package)
        except:
            if hasattr(pip, 'main'):
                pip.main(['install', package])
                print('Installed ' + package)
            else:
                pip._internal.main(['install', package])
                print('Installed ' + package)

        if package == 'Pillow':
            package = 'PIL'

        if package == 'scikit-image':
            package = 'skimage'

        globals()[package] = importlib.import_module(package)
