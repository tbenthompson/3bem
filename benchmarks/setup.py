import sys
import os
# Prepending the parent directory to the path is necessary to force the
# tbempy.* imports to import from this directory rather than a globally installed
# copy of tbempy
this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(this_dir, os.pardir))

from tbempy.testing import setup_tests
from tbempy.setup import files_in_dir, get_extension_config
from numpy.distutils.misc_util import Configuration
import subprocess
import os

benchmark_path = os.path.join('lib', 'benchmark')
benchmark_src_path = os.path.join(benchmark_path, 'src')
benchmark_include_path = os.path.join(benchmark_path, 'include')

def build_google_benchmark():
    cur_dir = os.getcwd()
    os.chdir(benchmark_path)
    subprocess.Popen(['cmake', '.']).wait()
    subprocess.Popen(['make']).wait()
    os.chdir(cur_dir)

def make_config(path, prefix):
    def configuration(parent_package='',top_path=None):
        build_google_benchmark()

        config = Configuration(path, parent_package, top_path)

        ext_config = get_extension_config()
        ext_config['sources'] += [
            f for f in files_in_dir(path, 'cpp')
            if f.startswith(os.path.join(path, prefix))
        ]
        ext_config['include_dirs'].append(benchmark_include_path)
        ext_config['extra_info']['libraries'].append('benchmark')
        ext_config['extra_info']['library_dirs'].append(benchmark_src_path)

        config.add_extension('runner', **ext_config)
        return config
    return configuration

if __name__ == "__main__":

    setup_tests(
        os.path.join('benchmarks', 'runner'),
        make_config('benchmarks', 'bench_')
    )
