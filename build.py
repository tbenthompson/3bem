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
from tools.fabricate import *

import os
import sys
import shutil
import subprocess

def files_in_dir(directory, ext):
    ret = []
    for file in os.listdir(directory):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.' + ext:
            ret.append(os.path.join(directory, file_name))
    return ret

dirs = ['3bem', 'test', 'inttest']
lib_srces = files_in_dir("3bem", "cpp")
tests = files_in_dir("test", "cpp")
inttests_cpp = files_in_dir("inttest", "cpp")
inttests_exec = files_in_dir("inttest", "cpp")
inttests_exec.remove('inttest/laplace')

compiler = 'mpic++'

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

cpp_flags = '-Wall -std=c++11 -fopenmp'.split()
cpp_flags.extend(['-I' + loc for loc in includes])

debug_flags = '-g -Og -DDEBUG=1'.split()
release_flags = '-DNDEBUG=1 -Ofast -ffast-math -funroll-loops'.split()
profile_flags = release_flags + ['-g']
test_coverage_flags = ['--coverage']
test_coverage_flags.extend(debug_flags)

flag_sets = dict()
flag_sets['test_coverage_flags'] = test_coverage_flags
flag_sets['debug_flags'] = debug_flags
flag_sets['release_flags'] = release_flags
flag_sets['profile_flags'] = profile_flags
cpp_flags.extend(flag_sets['debug_flags'])

lib_cpp_flags = ['-fPIC']
lib_cpp_flags.extend(cpp_flags)

link_flags = '--coverage -fopenmp -lhdf5'.split()
link_flags.append('-Wl,-rpath=' + petsc_dir + '/' + petsc_arch + '/lib')
link_flags.append('-L' + petsc_dir + '/' + petsc_arch + '/lib')
link_flags.append('-lpetsc')
link_flags.append('-larmadillo')

lib_link_flags = ['-shared']
lib_link_flags.extend(link_flags)

lib_dep_flags = ['-Wl,-rpath=./build', '-L./build', '-l3bem']

test_link_flags = ['-L../lib/unittest-cpp/builds']
test_link_flags.append('-lUnitTest++')
test_link_flags.extend(lib_dep_flags)
test_link_flags.extend(link_flags)

build_dir = 'build'

command_params = []

def entrypoint():
    save_parameters()
    main(parallel_ok = True, jobs = 12)

def save_parameters():
    if len(sys.argv) > 2:
        command_params.extend(sys.argv[2:])
        del sys.argv[2:]

def just_test():
    test_filename = command_params[0]
    test_rootname, _ = os.path.splitext(test_filename)
    setup_tree()
    compile_lib()
    compile_runner(test_rootname, cpp_flags)
    link_lib()
    link_runner(test_rootname, test_link_flags)
    run(oname(test_rootname))


def build():
    setup_tree()
    compile()
    link()

def setup_dir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def setup_tree():
    for d in dirs:
        setup_dir(os.path.join(build_dir, d))

def compile():
    compile_lib()
    compile_tests()
    compile_inttests()

def compile_lib():
    for source in lib_srces:
        compile_runner(source, lib_cpp_flags)

def compile_tests():
    for source in tests:
        compile_runner(source, cpp_flags)

def compile_inttests():
    for source in inttests_cpp:
        compile_runner(source, cpp_flags)

def compile_runner(source, flags):
    run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), flags)

def oname(filename):
    return os.path.join(build_dir, filename)

def link():
    link_lib()
    link_tests()
    link_inttests()

def link_lib():
    lib_objs = [oname(s + '.o') for s in lib_srces]
    after()
    run(compiler, '-o', oname('lib3bem.so'), lib_objs, lib_link_flags)
    after()

def link_tests():
    for source in tests:
        link_runner(source, test_link_flags)

def link_inttests():
    for source in inttests_exec:
        link_runner(source, test_link_flags,
                additional_objs = [oname('inttest/laplace.o')])

def link_runner(source, flags, additional_objs = []):
    objs = [oname(source + '.o')] + additional_objs
    run(compiler, '-o', oname(source), objs, flags)

def fast_tests():
    run_test_set(tests)

def slow_tests():
    run_test_set(inttests_exec)
    check_slow_tests()

def check_slow_tests():
    if os.path.exists('tools/__pycache__'):
        shutil.rmtree('tools/__pycache__')
    subprocess.call('\
        py.test -s \
        tools/check_planestrain.py \
        tools/check_antiplane.py\
    ', shell = True)

def run_test_set(test_names):
    for test_file in test_names:
        # I go outside fabricate's executor and use subprocess.call instead because
        # I see no reason to track the running of tests. Running all of the tests
        # each time is desirable. Using fabricate would only run the tests that
        # had changed since the last run.
        print("\nRunning test set: " + test_file)
        subprocess.call(oname(test_file))

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

if __name__ == "__main__":
    entrypoint()
