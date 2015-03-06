from tools.dev.config import get_config
from tools.dev.build import run_build
from tools.dev.testing import run_fast_tests, run_slow_tests
from tools.dev.fabricate import main
import sys

def fast_tests():
    run_fast_tests(get_config())

def slow_tests():
    run_slow_tests()

def lcov():
    coverage_file = oname('coverage.info')
    lcov_outdir = oname('lcov_out')
    after()
    run('lcov', '--capture', '--directory', build_dir, '--output-file', coverage_file)
    after()
    run('genhtml', coverage_file, '--output-directory', lcov_outdir)

def build():
    run_build(get_config(command_params))

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
