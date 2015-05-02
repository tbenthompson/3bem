import os
import shutil
import urllib
import subprocess

catch_url = 'https://raw.githubusercontent.com/philsquared/Catch/develop/single_include/catch.hpp'
boost_url = 'http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.bz2/download'
boost_numpy_url = 'https://github.com/ndarray/Boost.NumPy/archive/master.zip'
eigen_url = 'http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz'

# Delete the lib tree and recreate an empty directory
shutil.rmtree('lib')
os.makedirs('lib')

print("Downloading Catch unit testing framework")
urllib.urlretrieve(catch_url, "lib/_catch.hpp")

print("Downloading Boost for building C++ <--> python wrappers")
urllib.urlretrieve(boost_url, "boost.archive")
cmd = ['tar', '-xvf', 'boost.archive']
proc = subprocess.Popen(cmd)
proc.wait()
os.remove('boost.archive')
boost_directory_name = [f for f in os.listdir('./') if f.startswith('boost')][0]
shutil.move(boost_directory_name, 'lib/boost')

print("Download Boost.NumPy for clean C++ <--> python array transfer")
urllib.urlretrieve(boost_numpy_url, "boost.numpy.archive")
cmd = ['unzip', 'boost.numpy.archive']
proc = subprocess.Popen(cmd)
proc.wait()
os.remove('boost.numpy.archive')
numpy_directory_name = [f for f in os.listdir('./') if f.startswith('Boost')][0]
shutil.move(numpy_directory_name, 'lib/boost_numpy')

print("Downloading Eigen for linear algebra")
urllib.urlretrieve(eigen_url, "eigen.archive")
cmd = ['tar', '-xvf', 'eigen.archive']
proc = subprocess.Popen(cmd)
proc.wait()
os.remove('eigen.archive')
eigen_directory_name = [f for f in os.listdir('./') if f.startswith('eigen')][0]
shutil.move(eigen_directory_name, 'lib/eigen')
