import sys
import os
dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(dir, os.pardir))

from tbempy.setup import get_extension_config, setup_parallel_compile, files_in_dir
import numpy.distutils.command.build_ext as _build_ext
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

def get_test_srces():
    return files_in_dir('unit_tests', 'cpp')

def configuration(parent_package='',top_path=None):
    config = Configuration('unit_tests', parent_package, top_path)

    ext_config = get_extension_config()
    ext_config['sources'] += get_test_srces()

    config.add_extension('test_runner', **ext_config)
    return config

# distutils is good at building shared libraries. Here, I hack the system
# to make it produce an executable
class TestBuildExt(_build_ext.build_ext):
    def get_ext_filename(self, base):
        so_name = _build_ext.build_ext.get_ext_filename(self, base)
        return os.path.splitext(so_name)[0]

    def build_extension(self, ext):
        self._cxx_compiler.linker_so.remove('-shared')
        _build_ext.build_ext.build_extension(self, ext)

def setup_package():
    setup_parallel_compile()

    metadata = dict(name = 'tbempy')
    metadata['configuration'] = configuration
    metadata['cmdclass'] = dict(build_ext = TestBuildExt)

    # Default to building in place, the most common use case during development.
    if len(sys.argv) == 1:
        sys.argv = ['', 'build_ext', '--inplace']
    setup(**metadata)

if __name__ == "__main__":
    setup_package()
