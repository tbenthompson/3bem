from distutils.core import setup, Extension
from Cython.Build import cythonize
from tools.config import include_dirs, base_cpp_flags, cpp_flag_types
# from tools.util import files_in_dir

# lib_srces = files_in_dir('cpp', 'cpp')
# wrapper_srces = files_in_dir('python_wrapper', 'cpp')
# tbempy_srces = lib_srces + wrapper_srces
# tbempy_srces = [s + '.cpp' for s in tbempy_srces]
includes = ['./cpp', './lib/eigen']
cpp_flags = base_cpp_flags + cpp_flag_types['release']
link_args = [
    '-fopenmp', '-fPIC'
]
# tbempy = Extension(
#     'tbempy',
#     sources = tbempy_srces,
#     include_dirs = includes,
#     libraries = libs,
#     extra_compile_args = cpp_flags,
#     extra_link_args = link_args,
#     language = 'c++'
# )
tbempy = cythonize(Extension(
    'mesh',
    sources = ['cython_wrapper/mesh.pyx'],
    include_dirs = includes,
    extra_compile_args = cpp_flags,
    extra_link_args = link_args,
    language = 'c++'
))

setup(
    name = 'Name',
    ext_modules = tbempy
)
