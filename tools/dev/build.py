"""
This is the build script for 3bem. It uses the fabricate.py build tool.
Think "make", but way better. It automatically determines dependencies
between build steps so that if one file changes, only the minimal number of
other files are recompiled or relinked.

Also, it's nice to write a build script in python instead of some arcane
domain specific language like "make".

More info:
http://code.google.com/p/fabricate/
"""
from __future__ import print_function
from tools.fabricate import *

import os
import sys
import shutil
import subprocess

from tools.dev.testing import run_fast_tests, run_slow_tests
from tools.dev.util import files_in_dir, oname

def get_config():
    build_type = 'debug'
    if '-r' in command_params:
        build_type = 'release'
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

    lib_dep_flags = ['-Wl,-rpath=./build', '-L./build', '-l3bem']

    # c['test_link_flags'] = link_flags +\
    #     lib_dep_flags +\
    #     ['-L../lib/unittest-cpp/builds', '-lUnitTest++']

    lib = dict()
    lib['cpp_flags'] = cpp_flags + ['-fPIC']
    lib['link_flags'] = link_flags + ['-shared']
    lib['sources'] = files_in_dir('3bem', 'cpp')
    lib['linked_sources'] = []
    lib['binary_name'] = 'lib3bem.so'

    c = dict()
    c['build_dir'] = 'build_' + str(build_type)
    c['subdirs'] = ['3bem', 'test', 'inttest']
    c['compiler'] = 'mpic++'
    c['targets'] = [lib]
    return c


def just_test():
    t = test_info[command_params[0]]
    compile_flags = base_cpp_flags + flag_sets['debug_flags']
    compile_runner(t['src'], compile_flags)
    for s in t['lib_srcs']:
        compile_runner(s, base_cpp_flags + flag_sets['debug_flags'])
    after()
    link_runner([t['src']] + t['lib_srcs'], oname(t['src']), test_link_flags)

def build():
    c = get_config()
    setup_tree(c)
    for t in c['targets']:
        compile(c, t)
        link(c, t)

def setup_tree(c):
    def setup_dir(dirname):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    setup_dir(c['build_dir'])
    for d in c['subdirs']:
        setup_dir(os.path.join(c['build_dir'], d))

def compile(c, target):
    for source in target['sources']:
        compile_runner(c['build_dir'], c['compiler'], source, target['cpp_flags'])
    after()

def compile_runner(build_dir, compiler, source, flags):
    run(
        compiler,
        '-c',
        source + '.cpp',
        '-o',
        oname(build_dir, source + '.o'),
        flags
    )

def link(c, t):
    objs = [oname(c['build_dir'], s + '.o')
        for s in t['sources'] + t['linked_sources']]
    binary_path = oname(c['build_dir'], t['binary_name'])
    link_runner(c['compiler'], binary_path, objs, t['link_flags'])
    after()

def link_runner(compiler, binary_path, objs, flags):
    run(compiler, '-o', binary_path, objs, flags)

def fast_tests():
    run_fast_tests(get_build_dir())

def slow_tests():
    run_slow_tests()

def lcov():
    coverage_file = oname('coverage.info')
    lcov_outdir = oname('lcov_out')
    after()
    run('lcov', '--capture', '--directory', build_dir, '--output-file', coverage_file)
    after()
    run('genhtml', coverage_file, '--output-directory', lcov_outdir)

def clean():
    autoclean()

def rebuild():
    clean()
    build()

command_params = []
def save_parameters():
    if len(sys.argv) > 2:
        command_params.extend(sys.argv[2:])
        del sys.argv[2:]

def entrypoint(dir):
    save_parameters()
    main(parallel_ok = True, build_dir = dir, jobs = 12)

if __name__ == "__main__":
    entrypoint()
