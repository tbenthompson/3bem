from build_tools.config import get_config
from build_tools.build import build_target
from build_tools.fabricate import main, autoclean
import pprint
import sys
import os
import subprocess
import shutil

def tests():
    c = get_config(command_params)
    test_runner = os.path.join(c['build_dir'], c['targets']['tests']['binary_name'])
    p = subprocess.Popen([test_runner] + command_params)
    p.wait()
    return p.returncode

def basic_config():
    c = get_config(command_params)
    del c['targets']
    pprint.pprint(c)

def full_config():
    c = get_config(command_params)
    pprint.pprint(c)

def build():
    build_tests()
    build_python_wrapper()

def build_tests():
    c = get_config(command_params)
    build_target(c, c['targets']['tests'])

def build_python_wrapper():
    c = get_config(command_params)
    build_target(c, c['targets']['python_wrapper'])
    binary_path = os.path.join(
        c['build_dir'], c['targets']['python_wrapper']['binary_name']
    )
    shutil.copy(binary_path, c['subdirs']['tbempy_dir'])

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

def run_fabricate(dir):
    import multiprocessing

    n_cores = 2
    try:
        n_cores = multiprocessing.cpu_count()
    except NotImplementedError as e:
        n_cores = 2
    print('Building with ' + str(n_cores) + ' cores')
    main(parallel_ok = True, build_dir = dir, jobs = 12)

def entrypoint(dir):
    save_parameters()
    run_fabricate(dir)

if __name__ == '__main__':
    entrypoint()
