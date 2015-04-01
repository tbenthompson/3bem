from __future__ import print_function
from tools.build.util import files_in_dir, oname
from tools.build.test_info import unit_test_info, acceptance_test_info
import sys
import subprocess
import re
import os
from copy import copy

def testing_targets(test_info, loc, c):
    test_link_flags = c['link_flags'] + ['-L../lib/unittest-cpp/builds', '-lUnitTest++']
    ts = dict()
    for test_name, test_data in test_info.iteritems():
        target_link_flags = copy(test_link_flags)
        if test_data.get('link_lib', False):
            target_link_flags += c['lib_dep_flags']
        target = dict()
        target['cpp_flags'] = c['cpp_flags']
        target['link_flags'] = target_link_flags
        target['sources'] = [os.path.join(loc, test_data['src'])]
        target['sources'] += test_data.get('other_srces', [])
        target['linked_sources'] = [s for s in test_data['lib_srcs']]
        target['linked_sources_flags'] = c['targets']['lib']['cpp_flags']
        target['binary_name'] = test_data['src']
        target['priority'] = 1
        ts[test_name] = target
    return ts

def unit_testing_targets(c):
    return testing_targets(unit_test_info, 'test', c)

def acceptance_testing_targets(c):
    return testing_targets(acceptance_test_info, 'acctests', c)

def tests():
    return files_in_dir("test", "cpp")

def run_unit_tests(c):
    test_names = [unit_test_info[k]['src'] for k in unit_test_info]
    run_test_set(c, test_names)

def run_acceptance_tests(build_dir):
    run_test_set(build_dir, acceptance_test_info)
    check_acceptance_tests()

def check_acceptance_tests():
    if os.path.exists('tools/__pycache__'):
        shutil.rmtree('tools/__pycache__')
    subprocess.call('\
        py.test -s \
        tools/check_planestrain.py \
        tools/check_antiplane.py\
    ', shell = True)

def run_test(build_dir, test_file, n):
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
    logger(stdout.rstrip('\n'))
    if 'Success' in stdout:
        return interpret_success(stdout, stderr)
    if 'FAILURE' in stdout:
        return interpret_failure(stdout, stderr)
    return 1, 0, ''

def run_test_set(c, test_names, print_stdout = False):
    stderr = ''
    n_success = 0
    n_failure = 0

    for test_file in sorted(test_names):
        c['printer']("\nRunning tests: " + str(test_file))
        sys.stdout.flush()
        p = run_test(c['build_dir'], test_file, 123)
        s, f, errors = interpret_test_results(p, c['build_dir'], test_file, c['printer'])
        stderr += errors
        n_success += s
        n_failure += f
        c['printer']('')

    print('')
    print(stderr, end='')
    print("Tests successfully passed: " + str(n_success))
    print("Tests failed: " + str(n_failure))
