from distutils.core import setup, Extension
from tools.config import includes, base_cpp_flags, cpp_flag_types
from tools.util import files_in_dir

lib_srces = files_in_dir('cpp', 'cpp')
wrapper_srces = files_in_dir('python_wrapper', 'cpp')
tbempy_srces = lib_srces + wrapper_srces
tbempy_srces = [s + '.cpp' for s in tbempy_srces]
includes = ['./cpp', './lib/eigen']
cpp_flags = base_cpp_flags + cpp_flag_types['release']
libs = [
    'boost_python'
]
link_args = [
    '-fopenmp'
]
tbempy = Extension(
    'tbempy',
    sources = tbempy_srces,
    include_dirs = includes,
    libraries = libs,
    extra_compile_args = base_cpp_flags,
    extra_link_args = link_args
)

setup(
    name = 'tbempy',
    version = '1.0',
    description = 'This is a demo package',
    ext_modules = [tbempy]
)
