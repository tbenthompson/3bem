from distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
import numpy.distutils.command.build_ext as _build_ext
from tbempy.build_ext import tbempyBuildExt
import urllib
import shutil
import os
import sys
import copy
import subprocess

def download_libs():
    if os.path.exists('lib'):
        print('')
        print('Not downloading libraries. If libraries should be re-downloaded, delete the lib directory')
        print('')
        return

    catch_url = 'https://raw.githubusercontent.com/philsquared/Catch/develop/single_include/catch.hpp'
    boost_url = 'http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.bz2/download'
    boost_numpy_url = 'https://github.com/ndarray/Boost.NumPy/archive/master.zip'

    # Delete the lib tree and recreate an empty directory
    os.makedirs('lib')

    print('Downloading Catch unit testing framework')
    urllib.urlretrieve(catch_url, os.path.join('lib', '_catch.hpp'))

    print('Downloading Boost for building C++ <--> python wrappers')
    urllib.urlretrieve(boost_url, 'boost.archive')
    cmd = ['tar', '-xvf', 'boost.archive']
    proc = subprocess.Popen(cmd)
    proc.wait()
    os.remove('boost.archive')
    boost_directory_name = [f for f in os.listdir(os.curdir) if f.startswith('boost')][0]
    shutil.move(boost_directory_name, os.path.join('lib', 'boost'))

    print('Download Boost.NumPy for clean C++ <--> python array transfer')
    urllib.urlretrieve(boost_numpy_url, 'boost.numpy.archive')
    cmd = ['unzip', 'boost.numpy.archive']
    proc = subprocess.Popen(cmd)
    proc.wait()
    os.remove('boost.numpy.archive')
    numpy_directory_name = [f for f in os.listdir(os.curdir) if f.startswith('Boost')][0]
    shutil.move(numpy_directory_name, os.path.join('lib', 'boost_numpy'))

def configuration(parent_package='', top_path = None):
    config = Configuration(None, parent_package, top_path)
    config.add_subpackage('tbempy')
    return config

def setup_package():
    download_libs()

    metadata = dict(
        name = 'tbempy',
        version = '0.1',
        description = 'The black box boundary element method.',
        author = 'T. Ben Thompson',
        author_email = 't.ben.thompson@gmail.com'
    )
    metadata['cmdclass'] = dict(build_ext = tbempyBuildExt)
    metadata['configuration'] = configuration

    if len(sys.argv) == 1:
        sys.argv = ['', 'build_ext', '--inplace']
    setup(**metadata)

if __name__ == '__main__':
    setup_package()
