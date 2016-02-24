PREREQUISITES:

1)    Anaconda Python 2.7 is required to run the NIRSpec DRP. It can be downloaded and installed following         the directions found at https://www.continuum.io/downloads.

2)    Verify correct installation by typing python on the command line. The terminal will display the              version number and the distribution, i.e. Anaconda.

3)    The NIRSpec DRP software is maintained on GitHub. Download the code from https://github.com/2ichard/nirspec_drp.git and install by decompressing the folder in the directory of your choice.

EXECUTING NIRSPEC DRP:

1)    Simple command line usage: python nirspec_drp.py in_dir out_dir

      Note that you must provide the full directory path of the python file if your current working directory      is not the source code.  

      OH sky lines are used to refine the theoretical wavelength scale. These lines are listed in the              ir_ohlines.dat file which is included in the download. In case you would like to use a different file,       it can be specified by using an additional command line argument descibed in advanced usage.  

      in_dir: Path to directory which contains the raw FITS files to be reduced. This must include at least 1              object frame and 1 flat image.
      out_dir: Path to the root of the output directory to used for storage of generated data products

2)    Advanced Usage
    
      Additional command line arguments are available. These are listed and described after executing: python      nirspec_drp.py -h

