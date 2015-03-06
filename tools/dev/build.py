"""
This is the build script for 3bem. It uses the fabricate.py build tool.
Think "make", but way better. It automatically determines dependencies
between build steps so that if one file changes, only the minimal number of
other files are recompiled or relinked.

Also, it's nice to write a build script in python instead of some arcane
domain specific language like "make".

More info:
http://code.google.com/p/fabricate/

This build script separates the different executables and libraries to be
constructed into "targets", which are compiled individually.

Some basic command line parsing is done to determine which target to compile
and whether to do a debug, release, etc build. At some point, the
sophistication of the command line parsing may need to improve.
"""
from __future__ import print_function
from tools.fabricate import *

import os
import sys
import shutil
import subprocess

from tools.dev.testing import run_fast_tests, run_slow_tests, testing_targets
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
    c['subdirs'] = ['3bem', 'test', 'inttest']
    c['compiler'] = 'mpic++'
    c['targets'] = dict()
    c['targets']['lib'] = lib
    c['targets'].update(testing_targets(cpp_flags, link_flags, lib_dep_flags))
    return c

def determine_targets(c):
    targets = dict()
    if len(command_params) > 0:
        for entry in command_params:
            if not entry.startswith('-'):
                targets[entry] = c['targets'][entry]
    if len(targets) > 0:
        return targets
    else:
        return c['targets']

def priority_groupings(targets):
    buckets = dict()
    for t in targets.values():
        p = t['priority']
        if p in buckets:
            buckets[p].append(t)
        else:
            buckets[p] = [t]
    return buckets

def build():
    c = get_config()
    targets = determine_targets(c)
    buckets = priority_groupings(targets)
    setup_tree(c)
    for b in sorted(buckets.keys()):
        for t in buckets[b]:
            print('\nCompiling target: ' + t['binary_name'])
            compile(c, t)
    after()
    for b in sorted(buckets.keys()):
        for t in buckets[b]:
            print('\nLinking target: ' + t['binary_name'])
            link(c, t)
        after()

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

if __name__ == '__main__':
    entrypoint()
