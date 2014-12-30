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
from build.fabricate import *

import os
import subprocess

def files_in_dir(directory, ext):
    ret = []
    for file in os.listdir(directory):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.' + ext:
            ret.append(os.path.join(directory, file_name))
    return ret

dirs = ['3bem', 'examples', 'test', 'inttest']
lib_srces = files_in_dir("3bem", "cpp")
examples = files_in_dir("examples", "cpp")
tests = files_in_dir("test", "cpp")
inttests = files_in_dir("inttest", "cpp")

compiler = 'mpic++'

petsc_dir = os.environ['PETSC_DIR']
petsc_arch = os.environ['PETSC_ARCH']
includes = [
    './3bem',
    '../lib/',
    '../lib/unittest-cpp/src',
    '../lib/autocheck/include',
    petsc_dir + '/' + petsc_arch + '/include',
    petsc_dir + '/include'
]

cpp_flags = '-Wall -std=c++11 -fopenmp'.split()
cpp_flags.extend(['-I' + loc for loc in includes])

debug_flags = '-g -Og -DDEBUG=1'.split()
release_flags = '-DNDEBUG=1 -Ofast -ffast-math -funroll-loops'.split()
profile_flags = release_flags + ['-g']
test_coverage_flags = ['--coverage']
test_coverage_flags.extend(debug_flags)

# cpp_flags.extend(test_coverage_flags)
# cpp_flags.extend(debug_flags)
# cpp_flags.extend(release_flags)
cpp_flags.extend(profile_flags)

lib_cpp_flags = ['-fPIC']
lib_cpp_flags.extend(cpp_flags)

link_flags = '--coverage -fopenmp -lhdf5'.split()
link_flags.append('-Wl,-rpath=' + petsc_dir + '/' + petsc_arch + '/lib')
link_flags.append('-L' + petsc_dir + '/' + petsc_arch + '/lib')
link_flags.append('-lpetsc')

lib_link_flags = ['-shared']
lib_link_flags.extend(link_flags)

lib_dep_flags = ['-Wl,-rpath=./build', '-L./build', '-l3bem']

test_link_flags = ['-L../lib/unittest-cpp']
test_link_flags.append('-lUnitTest++')
test_link_flags.extend(lib_dep_flags)
test_link_flags.extend(link_flags)

example_link_flags = lib_dep_flags
example_link_flags.extend(link_flags)
example_link_flags.extend('-lglut -lGL -lGLU -lGLEW'.split())

build_dir = 'build'

def setup_dir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def build():
    for d in dirs:
        setup_dir(os.path.join(build_dir, d))
    compile()
    link()

def oname(filename):
    return os.path.join(build_dir, filename)

def compile():
    for source in lib_srces:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), lib_cpp_flags)
    for source in examples:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), cpp_flags)
    for source in tests:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), cpp_flags)
    for source in inttests:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), cpp_flags)

def link():
    lib_objs = [oname(s + '.o') for s in lib_srces]
    after()
    run(compiler, '-o', oname('lib3bem.so'), lib_objs, lib_link_flags)
    after()
    for source in examples:
        run(compiler, '-o', oname(source), oname(source + '.o'), example_link_flags)
    for source in tests:
        run(compiler, '-o', oname(source), oname(source + '.o'), test_link_flags)
    for source in inttests:
        run(compiler, '-o', oname(source), oname(source + '.o'), test_link_flags)

def run_test_set(test_names):
    for test_file in test_names:
        # I go outside fabricate's runner and use subprocess.call instead because
        # I see no reason to track the running of tests. Running all of the tests
        # each time is desirable. Using fabricate would only run the tests that
        # had changed since the last run.
        print("\nRunning test set: " + test_file)
        subprocess.call(oname(test_file))

def run_tests():
    run_test_set(tests)

def run_integration_tests():
    run_test_set(inttests)

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

main(parallel_ok = True, jobs = 12)
