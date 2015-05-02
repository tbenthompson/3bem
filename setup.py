from distutils.core import setup, Extension
from tools.config import boost_include_dirs, boost_lib, boost_sources, include_dirs,\
    base_cpp_flags,cpp_flag_types, wrapper_srces, lib_srces
from tools.util import files_in_dir
import os

# Setting OPT to '' prevents distutils from appending the -Wstrict-prototypes flag
# which does not apply to c++ compilation
if 'OPT' not in os.environ:
    os.environ['OPT'] = ''

tbempy_srces = lib_srces + wrapper_srces + boost_sources
tbempy_srces = [s + '.cpp' for s in tbempy_srces]
no_warnings = True
cpp_flags = base_cpp_flags + cpp_flag_types['release']
if no_warnings:
    cpp_flags.append('-w')
includes = include_dirs + boost_include_dirs
libs = [
]
link_args = [
    '-fopenmp'
]

tbempy = Extension(
    'tbempy',
    sources = tbempy_srces,
    include_dirs = includes,
    libraries = libs,
    extra_compile_args = cpp_flags,
    extra_link_args = link_args,
    language = 'c++'
)
print(dir(tbempy))

# monkey-patch for parallel compilation
def parallelCCompile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    import multiprocessing.pool
    import multiprocessing
    N = 1
    try:
        N = multiprocessing.cpu_count()
    except e:
        pass
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects
import distutils.ccompiler
distutils.ccompiler.CCompiler.compile=parallelCCompile

setup(
    name = 'tbempy',
    version = '0.1',
    description = 'The black box boundary element method.',
    author = 'T. Ben Thompson',
    author_email = 't.ben.thompson@gmail.com',
    ext_modules = [tbempy]
)
