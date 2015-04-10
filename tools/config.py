from __future__ import print_function
from tools.util import files_in_dir
from tools.testing import unit_testing_targets
from numpy.distutils.system_info import get_info as np_config_info
import warnings
import os

def get_config(command_params):
    build_type = 'debug'
    if '-r' in command_params:
        build_type = 'release'
    if '-p' in command_params:
        build_type = 'profile'
    printer = lambda x: None
    if '-v' in command_params:
        printer = print
    src_dir = 'cpp'
    test_dir = 'test'
    py_wrap_dir = 'python_wrapper'

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        blas_lapack_info = np_config_info('lapack_opt', 0)

    includes = [
        './' + str(src_dir),
        '../lib/unittest-cpp/UnitTest++',
        '../lib/autocheck/include'
    ]

    base_cpp_flags = [
        '-Wextra',
        '-std=c++11',
        '-fopenmp',
        '-DDEBUG=1'
    ] + ['-I' + loc for loc in includes]

    flag_types = dict()
    flag_types['debug'] = ['-g', '-Og']
    flag_types['release'] = ['-Ofast']
    flag_types['profile'] = flag_types['release'] + ['-g']
    flag_types['coverage'] = ['--coverage'] + flag_types['debug']

    cpp_flags = base_cpp_flags + flag_types[build_type]

    link_flags = ['--coverage', '-fopenmp']
    link_flags.extend(['-L' + l_dir for l_dir in blas_lapack_info['library_dirs']])
    link_flags.extend(['-l' + library for library in blas_lapack_info['libraries']])

    lib = dict()
    lib['cpp_flags'] = cpp_flags + ['-fPIC']
    lib['link_flags'] = link_flags + ['-shared']
    lib['sources'] = files_in_dir(src_dir, 'cpp')
    lib['linked_sources'] = []
    lib['binary_name'] = 'lib3bem.so'
    lib['priority'] = 0

    python_wrapper = dict()
    python_wrapper['cpp_flags'] = lib['cpp_flags'] + ['-I/usr/include/python2.7']
    python_wrapper['link_flags'] = lib['link_flags'] + ['-lpython2.7', '-lboost_python']
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

