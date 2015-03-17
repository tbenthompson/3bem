from __future__ import print_function
from tools.build.util import files_in_dir
from tools.build.testing import unit_testing_targets, acceptance_testing_targets
import os

def get_config(command_params):
    build_type = 'debug'
    if '-r' in command_params:
        build_type = 'release'
    printer = lambda x: None
    if '-v' in command_params:
        printer = print
    petsc_dir = os.environ['PETSC_DIR']
    petsc_arch = os.environ['PETSC_ARCH']


    includes = [
        './3bem',
        '../lib/',
        '../lib/unittest-cpp/UnitTest++',
        '../lib/autocheck/include',
        petsc_dir + '/' + petsc_arch + '/include',
        petsc_dir + '/include'
    ]

    base_cpp_flags = [
        '-Wall',
        '-std=c++11',
        '-fopenmp',
        '-DDEBUG=1'
    ] + ['-I' + loc for loc in includes]

    flag_types = dict()
    flag_types['debug'] = ['-g', '-Og']
    flag_types['release'] = ['-Ofast','-ffast-math','-funroll-loops']
    flag_types['coverage'] = ['--coverage'] + flag_types['debug']

    cpp_flags = base_cpp_flags + flag_types[build_type]

    link_flags = [
        '--coverage',
        '-fopenmp',
        '-lhdf5',
        '-larmadillo',
        '-Wl,-rpath=' + petsc_dir + '/' + petsc_arch + '/lib',
        '-L' + petsc_dir + '/' + petsc_arch + '/lib',
        '-lpetsc'
        ]

    lib = dict()
    lib['cpp_flags'] = cpp_flags + ['-fPIC']
    lib['link_flags'] = link_flags + ['-shared']
    lib['sources'] = files_in_dir('3bem', 'cpp')
    lib['linked_sources'] = []
    lib['binary_name'] = 'lib3bem.so'
    lib['priority'] = 0

    build_dir = 'build_' + str(build_type)
    lib_dep_flags = ['-Wl,-rpath=./' + build_dir, '-L./' + build_dir, '-l3bem']

    c = dict()
    c['build_dir'] = build_dir
    c['subdirs'] = ['3bem', 'test', 'acctests']
    c['compiler'] = 'mpic++'
    c['command_params'] = command_params
    c['printer'] = printer
    c['lib_dep_flags'] = lib_dep_flags
    c['cpp_flags'] = cpp_flags
    c['link_flags'] = link_flags
    c['targets'] = dict()
    c['targets']['lib'] = lib
    c['targets'].update(unit_testing_targets(c))
    c['targets'].update(acceptance_testing_targets(c))
    return c

