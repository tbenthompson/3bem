import os
import shutil
import urllib
import subprocess

catch_url = 'https://raw.githubusercontent.com/philsquared/Catch/develop/single_include/catch.hpp'
boost_url = 'http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.bz2/download'
eigen_url = 'http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz'

# Delete the lib tree and recreate an empty directory
shutil.rmtree('lib')
os.makedirs('lib')

print("Downloading Catch unit testing framework")
urllib.urlretrieve(catch_url, "lib/_catch.hpp")

print("Downloading Boost for building C++ --> python wrappers")
urllib.urlretrieve(boost_url, "boost.archive")
cmd = ['tar', '-xvf', 'boost.archive']
proc = subprocess.Popen(cmd)
proc.wait()
os.remove('boost.archive')
boost_directory_name = [f for f in os.listdir('./') if f.startswith('boost')][0]
print(boost_directory_name)
shutil.move(boost_directory_name, 'lib/boost')

print("Downloading Eigen for linear algebra")
urllib.urlretrieve(eigen_url, "eigen.archive")
cmd = ['tar', '-xvf', 'eigen.archive']
proc = subprocess.Popen(cmd)
proc.wait()
os.remove('eigen.archive')
eigen_directory_name = [f for f in os.listdir('./') if f.startswith('eigen')][0]
shutil.move(eigen_directory_name, 'lib/eigen')
