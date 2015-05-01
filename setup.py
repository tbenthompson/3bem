from distutils.core import setup, Extension
from tools.config import boost_include_dir, boost_lib, boost_sources, include_dirs,\
    base_cpp_flags,cpp_flag_types
from tools.util import files_in_dir

lib_srces = files_in_dir('cpp', 'cpp')
wrapper_srces = files_in_dir('python_wrapper', 'cpp')
tbempy_srces = lib_srces + wrapper_srces + boost_sources
tbempy_srces = [s + '.cpp' for s in tbempy_srces]
cpp_flags = base_cpp_flags + cpp_flag_types['release']
cpp_flags += ['-Wno-unused-variable']
includes = include_dirs + [boost_include_dir]
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

setup(
    name = 'tbempy',
    version = '0.1',
    description = 'The black box boundary element method.',
    author = 'T. Ben Thompson',
    author_email = 't.ben.thompson@gmail.com',
    ext_modules = [tbempy]
)
