from tools.build.config import get_config
from tools.build.build import run_build
from tools.build.testing import run_unit_tests, run_acceptance_tests
from tools.build.fabricate import main
from codegen.main import generate
import sys

def unit_tests():
    run_unit_tests(get_config(command_params))

def acceptance_tests():
    run_acceptance_tests(get_config(command_params))

def lcov():
    coverage_file = oname('coverage.info')
    lcov_outdir = oname('lcov_out')
    after()
    run('lcov', '--capture', '--directory', build_dir, '--output-file', coverage_file)
    after()
    run('genhtml', coverage_file, '--output-directory', lcov_outdir)

def codegen():
    generate()

def build():
    codegen()
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
