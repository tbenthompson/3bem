from __future__ import print_function
from tools.util import files_in_dir
from tools.testing import unit_testing_targets
from numpy.distutils.system_info import get_info as np_config_info
import distutils.sysconfig
import warnings
import os

def include(dirs):
    return ['-I' + d for d in dirs]

src_dir = 'cpp'
test_dir = 'test'
py_wrap_dir = 'python_wrapper'

include_dirs = ['./' + str(src_dir), './lib', './lib/eigen']
include_flags = include(include_dirs)

base_cpp_flags = ['-Wextra', '-std=c++11', '-fopenmp', '-DDEBUG=1']
cpp_flag_types = dict()
cpp_flag_types['debug'] = ['-g', '-Og']
cpp_flag_types['release'] = ['-Ofast']
cpp_flag_types['profile'] = cpp_flag_types['release'] + ['-g']
cpp_flag_types['coverage'] = ['--coverage'] + cpp_flag_types['debug']

base_link_flags = ['-fopenmp']
link_flag_types = dict()
link_flag_types['coverage'] = ['--coverage']

boost_python_dir = ''
boost_python_lib = 'boost_python'
python_include_dir = distutils.sysconfig.get_python_inc()
python_lib_dir = distutils.sysconfig.get_python_lib(standard_lib = True)
python_lib = distutils.sysconfig.get_config_var('LDLIBRARY')

def get_config(command_params):
    build_type = 'release'
    if '-d' in command_params:
        build_type = 'debug'
    if '-p' in command_params:
        build_type = 'profile'
    printer = lambda x: None
    if '-v' in command_params:
        printer = print

    cpp_flags = base_cpp_flags + include_flags + cpp_flag_types[build_type]

    link_flags = base_link_flags + link_flag_types.get(build_type, [])

    lib = dict()
    lib['cpp_flags'] = cpp_flags + ['-fPIC']
    lib['link_flags'] = link_flags + ['-shared']
    lib['sources'] = files_in_dir(src_dir, 'cpp')
    lib['linked_sources'] = []
    lib['binary_name'] = 'lib3bem.so'
    lib['priority'] = 0

    python_wrapper_cpp_flags = lib['cpp_flags'] + include([python_include_dir])
    python_wrapper_link_flags = lib['link_flags'] +\
        ['-L' + python_lib_dir] +\
        ['-l:' + python_lib] +\
        ['-l' + boost_python_lib]
    python_wrapper = dict()
    python_wrapper['cpp_flags'] = python_wrapper_cpp_flags
    python_wrapper['link_flags'] = python_wrapper_link_flags
    python_wrapper['sources'] = files_in_dir(py_wrap_dir, 'cpp')
    python_wrapper['linked_sources'] = lib['sources']
    python_wrapper['binary_name'] = 'tbempy.so'
    python_wrapper['priority'] = 2

    build_dir = os.path.join('build', str(build_type))
    lib_dep_flags = ['-Wl,-rpath=./' + build_dir, '-L./' + build_dir, '-l3bem']

    c = dict()
    c['build_dir'] = build_dir
    c['subdirs'] = dict(
        src_dir = src_dir,
        test_dir = test_dir,
        py_wrap_dir = py_wrap_dir
    )
    c['compiler'] = 'g++'
    c['command_params'] = command_params
    c['printer'] = printer
    c['lib_dep_flags'] = lib_dep_flags
    c['cpp_flags'] = cpp_flags
    c['link_flags'] = link_flags
    c['targets'] = dict()
    c['targets']['lib'] = lib
    c['targets']['python_wrapper'] = python_wrapper
    c['targets'].update(unit_testing_targets(c))
    return c

