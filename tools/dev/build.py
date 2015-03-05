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

def files_in_dir(directory, ext):
    ret = []
    for file in os.listdir(directory):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.' + ext:
            ret.append(os.path.join(directory, file_name))
    return ret

def oname(filename):
    return os.path.join(build_dir, filename)


dirs = ['3bem', 'test', 'inttest']
lib_srces = files_in_dir("3bem", "cpp")
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

base_cpp_flags = '-Wall -std=c++11 -fopenmp'.split()
base_cpp_flags.extend(['-I' + loc for loc in includes])

debug_flags = '-g -Og -DDEBUG=1'.split()
release_flags = '-DDEBUG=1 -Ofast -ffast-math -funroll-loops'.split()
profile_flags = release_flags + ['-g']
test_coverage_flags = ['--coverage']
test_coverage_flags.extend(debug_flags)

flag_sets = dict()
flag_sets['test_coverage_flags'] = test_coverage_flags
flag_sets['debug_flags'] = debug_flags
flag_sets['release_flags'] = release_flags
flag_sets['profile_flags'] = profile_flags
cpp_flags = base_cpp_flags + flag_sets['release_flags']

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

def tests():
    return files_in_dir("test", "cpp")

def build_executables():
    executables = []
    for t in tests():
        test_data = dict()
        test_data['source_files'] = [t]
        test_data['cpp_flags'] = cpp_flags + flag_sets['release_flags']
        test_data['exec_name'] = oname(t)
        test_data['link_flags'] = test_link_flags
        if 'regression_021515' in t:
            test_data['cpp_flags'] = base_cpp_flags + ['-DDEBUG=1', '-O3']
        executables.append(test_data)
    return executables

executables = build_executables()
command_params = []

test_info = dict()
test_info['function'] = dict()
test_info['function']['src'] = 'test/test_function'
test_info['function']['lib_srcs'] = ['3bem/vectorx']
test_info['petsc'] = dict()
test_info['petsc']['src'] = 'test/test_petsc'
test_info['petsc']['lib_srcs'] = ['3bem/petsc_facade']
test_info['matrix_free'] = dict()
test_info['matrix_free']['src'] = 'test/test_matrix_free_builder'
test_info['matrix_free']['lib_srcs'] = []
test_info['quadrature'] = dict()
test_info['quadrature']['src'] = 'test/test_quadrature'
test_info['quadrature']['lib_srcs'] = ['3bem/quadrature']

def just_test():
    t = test_info[command_params[0]]
    compile_flags = base_cpp_flags + flag_sets['debug_flags']
    compile_runner(t['src'], compile_flags)
    for s in t['lib_srcs']:
        compile_runner(s, base_cpp_flags + flag_sets['debug_flags'])
    after()
    link_runner([t['src']] + t['lib_srcs'], oname(t['src']), test_link_flags)

def entrypoint(dir):
    save_parameters()
    main(parallel_ok = True, build_dir = dir, jobs = 12)

def save_parameters():
    if len(sys.argv) > 2:
        command_params.extend(sys.argv[2:])
        del sys.argv[2:]

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
    for e in executables:
        for s in e['source_files']:
            compile_runner(s, e['cpp_flags'])

def compile_inttests():
    for source in inttests_cpp:
        compile_runner(source, cpp_flags)

def compile_runner(source, flags):
    run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), flags)

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
    for e in executables:
        link_runner(e['source_files'], e['exec_name'], e['link_flags'])

def link_inttests():
    for source in inttests_exec:
        link_runner([source], oname(source), test_link_flags,
                additional_objs = [oname('inttest/laplace.o')])

def link_runner(sources, exec_name, flags, additional_objs = []):
    objs = [oname(s + '.o') for s in sources] + additional_objs
    run(compiler, '-o', exec_name, objs, flags)

def fast_tests():
    run_test_set(tests())

def slow_tests():
    run_test_set(inttests_exec, True)
    check_slow_tests()

def check_slow_tests():
    if os.path.exists('tools/__pycache__'):
        shutil.rmtree('tools/__pycache__')
    subprocess.call('\
        py.test -s \
        tools/check_planestrain.py \
        tools/check_antiplane.py\
    ', shell = True)

def run_test_set(test_names, print_stdout = False):
    stderr = ''
    n_success = 0
    n_failure = 0
    for test_file in sorted(test_names):
        sys.stdout.flush()
        p = subprocess.Popen(
            oname(test_file),
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        p.wait()
        stdout = p.stdout.read()
        if print_stdout:
            print(stdout)
        import re
        m = re.search('Success: ([0-9]+) tests passed', stdout)
        if m is not None:
            s = int(m.group(1))
            n_success += s
            print('.' * s, end='')
            continue
        stderr += p.stderr.read()
        m = re.search('FAILURE: ([0-9]+) out of ([0-9]+) tests failed', stdout)
        if m is not None:
            f = int(m.group(1))
            s = int(m.group(2))
            n_success += s
            n_failure += f
            print('F' * f + '.' * s, end='')
    print('')
    print(stderr, end='')
    print("Tests successfully passed: " + str(n_success))
    print("Tests failed: " + str(n_failure))

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
