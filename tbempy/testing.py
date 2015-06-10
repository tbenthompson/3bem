from build_ext import tbempyBuildExt
from setup import setup_parallel_compile, files_in_dir, get_extension_config
from download import download_libs
import os
import sys
from numpy.distutils.core import setup

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

def setup_tests(executable_path, config_fnc):
    download_libs()
    setup_parallel_compile()

    metadata = dict(name = 'tbempy')
    metadata['configuration'] = config_fnc
    metadata['cmdclass'] = dict(build_ext = TestBuildExt)

    # Default to building in place, the most common use case during development.
    # Removing the existing executable is necessary in order to force a rebuild
    # a rebuild is desirable because distutils does a poor job detecting changes
    if os.path.exists(executable_path):
        os.remove(executable_path)
    sys.argv = ['__', 'build_ext', '--inplace']
    setup(**metadata)
