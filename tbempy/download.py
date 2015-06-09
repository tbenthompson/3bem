import urllib
import shutil
import os
import sys
import copy
import subprocess

def download_libs():
    if not os.path.exists('lib'):
        os.makedirs('lib')
    print(os.listdir('lib'))
    print(os.listdir('lib'))
    print(os.listdir('lib'))
    print(os.listdir('lib'))
    if not os.path.exists('lib/_catch.hpp'):
        download_catch()
    else:
        print('Not downloading Catch.')
    if not os.path.exists('lib/boost'):
        download_boost()
    else:
        print('Not downloading boost.')
    if not os.path.exists('lib/boost_numpy'):
        download_boost_numpy()
    else:
        print('Not downloading boost_numpy.')
    if not os.path.exists('lib/gte'):
        download_gte()
    else:
        print('Not downloading Geometric Tools Engine.')
    print('To redownload a library, delete the lib directory')

def download_catch():
    catch_url = 'https://raw.githubusercontent.com/philsquared/Catch/develop/single_include/catch.hpp'
    print('Downloading Catch unit testing framework')
    urllib.urlretrieve(catch_url, os.path.join('lib', '_catch.hpp'))

def download_boost():
    boost_url = 'http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.bz2/download'
    print('Downloading Boost for building C++ <--> python wrappers')
    urllib.urlretrieve(boost_url, 'boost.archive')
    cmd = ['tar', '-xvf', 'boost.archive']
    proc = subprocess.Popen(cmd)
    proc.wait()
    os.remove('boost.archive')
    boost_directory_name = [
        f for f in os.listdir(os.curdir) if f.startswith('boost')
    ][0]
    shutil.move(boost_directory_name, os.path.join('lib', 'boost'))

def download_boost_numpy():
    boost_numpy_url = 'https://github.com/ndarray/Boost.NumPy/archive/master.zip'
    print('Download Boost.NumPy for clean C++ <--> python array transfer')
    urllib.urlretrieve(boost_numpy_url, 'boost.numpy.archive')
    cmd = ['unzip', 'boost.numpy.archive']
    proc = subprocess.Popen(cmd)
    proc.wait()
    os.remove('boost.numpy.archive')
    numpy_directory_name = [
        f for f in os.listdir(os.curdir) if f.startswith('Boost')
    ][0]
    shutil.move(numpy_directory_name, os.path.join('lib', 'boost_numpy'))

def download_gte():
    gte_url = 'http://www.geometrictools.com/Downloads/GeometricToolsEngine1p14.zip'
    print('Downloading Geometric Tools Engine')
    urllib.urlretrieve(gte_url, 'gte.archive')
    cmd = ['unzip', 'gte.archive']
    proc = subprocess.Popen(cmd)
    proc.wait()
    os.remove('gte.archive')
    shutil.move('GeometricTools/GTEngine', os.path.join('lib', 'gte'))
    os.rmdir('GeometricTools')

