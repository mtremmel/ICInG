from distutils.core import setup

DESCRIPTION = 'Tools for creating initial conditions for an isolated cosmological collapse simulation '
LONG_DESCRIPTION = open('README.md').read()
NAME = 'ICInG'
VERSION = '1.0'
AUTHOR = 'Michael Tremmel'
AUTHOR_EMAIL = 'm.tremmel6@gmail.com'
MAINTAINER = 'Michael Tremmel'
MAINTAINER_EMAIL = 'm.tremmel6@gmail.com'
URL = 'http://github.com/mtremmel/ICInG'
DOWNLOAD_URL = 'http://github.com/mtremmel/ICInG'
LICENSE = 'BSD'



setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR, 
      author_email= AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      package_dir={'ICInG/':''},
      packages=['ICInG'], 
      package_data={'ICInG':['data/snapshot_005.hdf5']},
      classifiers = ["Development Status :: 3 - Alpha",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
                     "Programming Language :: Python :: 2",
                     "Topic :: Scientific/Engineering :: Astronomy"])
