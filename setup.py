from distutils.core import setup
import build_tools.entry
import os
import sys
import copy

# I want to call fabricate with sys.argv empty, but still pass sys.argv to
# the python packaging tools
stored_argv = copy.copy(sys.argv)
sys.argv = []

# Fabricate calls system.exit when it is done compiling. Let's prevent it
# from actually exiting by catching the SystemExit exception
try:
    # The root of the build tree should be the containing folder for this script
    dir = os.path.dirname(os.path.realpath(__file__))
    build_tools.entry.run_fabricate('build', dir)
except SystemExit as e:
    pass

# Check the build by running unit tests
print('Checking the build by running the unit tests')
returncode = build_tools.entry.tests()
if returncode != 0:
    print('')
    print('')
    print('')
    print('**** Tests are failing. The build or code have problems. Quitting. **** ')
    print('')
    print('')
    print('')
    sys.exit(returncode)

# Restore sys.argv to be passed to python packager
sys.argv = stored_argv

# Call python packager
setup(
    name = 'tbempy',
    version = '0.1',
    description = 'The black box boundary element method.',
    author = 'T. Ben Thompson',
    author_email = 't.ben.thompson@gmail.com',
    packages = ['tbempy'],
    package_data = dict(tbempy = ['tbempy.so'])
)
