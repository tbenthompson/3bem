import sys
import os
# Prepending the parent directory to the path is necessary to force the
# tbempy.* imports to import from this directory rather than a globally installed
# copy of tbempy
this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(this_dir, os.pardir))

from tbempy.download import download_libs
from tbempy.setup import get_extension_config, setup_parallel_compile, \
    files_in_dir
from tbempy.build_ext import tbempyBuildExt
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

def get_test_srces():
    srces = files_in_dir('unit_tests', 'cpp')
    return srces

def configuration(parent_package='',top_path=None):
    config = Configuration('unit_tests', parent_package, top_path)

    ext_config = get_extension_config()
    ext_config['sources'] += get_test_srces()

    config.add_extension('test_runner', **ext_config)
    return config

# distutils is good at building shared libraries. Here, hack the system
# to make it produce an executable
class TestBuildExt(tbempyBuildExt):
    def get_ext_filename(self, base):
        so_name = tbempyBuildExt.get_ext_filename(self, base)
        return os.path.splitext(so_name)[0]

    def build_extension(self, ext):
        def remove_param(name):
            if name in self._cxx_compiler.linker_so:
                self._cxx_compiler.linker_so.remove(name)

        # Shared libraries are built with "-bundle" on mac. Don't want that!
        remove_param('-bundle')
        # And "-shared" on linux. Don't want that!
        remove_param('-shared')
        tbempyBuildExt.build_extension(self, ext)

def setup_package():
    download_libs()
    setup_parallel_compile()

    metadata = dict(name = 'tbempy')
    metadata['configuration'] = configuration
    metadata['cmdclass'] = dict(build_ext = TestBuildExt)

    # Default to building in place, the most common use case during development.
    # Removing the existing executable is necessary in order to force a rebuild
    # a rebuild is desirable because distutils does a poor job detecting changes
    if len(sys.argv) == 1:
        executable_path = os.path.join('unit_tests', 'test_runner')
        if os.path.exists(executable_path):
            os.remove(executable_path)
        sys.argv.extend(['build_ext', '--inplace'])
    setup(**metadata)

if __name__ == "__main__":
    setup_package()
