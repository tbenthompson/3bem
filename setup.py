from distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
import numpy.distutils.command.build_ext as _build_ext
from tbempy.build_ext import tbempyBuildExt
from tbempy.download import download_libs
import sys
import os
import shutil

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
        tbempy_path = os.path.join('tbempy', '_tbempy.so')
        if os.path.exists(tbempy_path):
            os.remove(tbempy_path)
        sys.argv.extend(['build_ext', '--inplace'])
    elif sys.argv[1] == 'reinstall':
        if os.path.exists('build'):
            shutil.rmtree('build')
        sys.argv[1] = 'install'
    setup(**metadata)

if __name__ == '__main__':
    setup_package()
