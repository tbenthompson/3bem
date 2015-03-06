from __future__ import print_function
from tools.build.util import files_in_dir, oname
import sys
import subprocess
import re
import os

test_info = dict()
test_info['function'] = dict()
test_info['function']['src'] = 'test_function'
test_info['function']['lib_srcs'] = ['vectorx']
test_info['petsc'] = dict()
test_info['petsc']['src'] = 'test_petsc'
test_info['petsc']['lib_srcs'] = ['petsc_facade', 'vectorx', 'util']
test_info['matrix_free'] = dict()
test_info['matrix_free']['src'] = 'test_matrix_free_builder'
test_info['matrix_free']['lib_srcs'] = []
test_info['matrix_free']['link_lib'] = True
test_info['quadrature'] = dict()
test_info['quadrature']['src'] = 'test_quadrature'
test_info['quadrature']['lib_srcs'] = ['quadrature']
test_info['closest_pt'] = dict()
test_info['closest_pt']['src'] = 'test_closest_pt'
test_info['closest_pt']['lib_srcs'] = ['closest_pt']
test_info['integral_term'] = dict()
test_info['integral_term']['src'] = 'test_integral_term'
test_info['integral_term']['lib_srcs'] = []

acctest_info = dict()

def testing_targets(cpp_flags, link_flags, lib_dep_flags):
    test_link_flags = link_flags + ['-L../lib/unittest-cpp/builds', '-lUnitTest++']
    ts = dict()
    for test_name, test_data in test_info.iteritems():
        target_link_flags = test_link_flags
        if test_data.get('link_lib', False):
            target_link_flags += lib_dep_flags
        target = dict()
        target['cpp_flags'] = cpp_flags
        target['link_flags'] = test_link_flags
        target['sources'] = [os.path.join('test', test_data['src'])]
        target['linked_sources'] =\
            [os.path.join('3bem', s) for s in test_data['lib_srcs']]
        target['binary_name'] = test_data['src']
        target['priority'] = 1000
        ts[test_name] = target
    return ts

def tests():
    return files_in_dir("test", "cpp")

def run_fast_tests(build_dir):
    run_test_set(build_dir, tests())

def run_slow_tests(build_dir):
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

def run_test(build_dir, test_file):
    p = subprocess.Popen(
        oname(build_dir, test_file),
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    p.wait()
    return p

def interpret_success(stdout, stderr):
    m = re.search('Success: ([0-9]+) tests passed', stdout)
    s = int(m.group(1))
    print('.' * s, end='')
    return s, 0, ''

def interpret_failure(stdout, stderr):
    m = re.search('FAILURE: ([0-9]+) out of ([0-9]+) tests failed', stdout)
    f = int(m.group(1))
    s = int(m.group(2))
    print('F' * f + '.' * s, end='')
    return s, f, stderr

def interpret_test_results(p, build_dir, test_file, logger):
    stdout = p.stdout.read()
    stderr = p.stderr.read()
    logger(stdout)
    if 'Success' in stdout:
        return interpret_success(stdout, stderr)
    if 'FAILURE' in stdout:
        return interpret_failure(stdout, stderr)

def run_test_set(build_dir, test_names, print_stdout = False):
    stderr = ''
    n_success = 0
    n_failure = 0

    logger = lambda x: None
    if print_stdout:
        logger = lambda x: print(x)

    for test_file in sorted(test_names):
        logger("Running tests: " + str(test_file))
        sys.stdout.flush()
        p = run_test(build_dir, test_file)
        s, f, errors = interpret_test_results(p, build_dir, test_file, logger)
        stderr += errors
        n_success += s
        n_failure += f

    print('')
    print(stderr, end='')
    print("Tests successfully passed: " + str(n_success))
    print("Tests failed: " + str(n_failure))
