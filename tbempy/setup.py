import os
import sys
import multiprocessing.pool
import multiprocessing
import distutils.ccompiler
from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import get_info as np_config_info
import warnings

# Patch for parallel compilation with distutils
# From: http://stackoverflow.com/questions/11013851/speeding-up-build-process-with-distutils
def parallelCCompile(self, sources, output_dir = None, macros = None,
        include_dirs = None, debug = 0, extra_preargs = None, extra_postargs = None,
        depends = None):

    # these lines are copied directly from distutils.ccompiler.CCompiler
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs
    )
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    # Determine the number of compilation threads. Unless there are special
    # circumstances, this is the number of cores on the machine
    N = 1
    try:
        N = multiprocessing.cpu_count()
    except e:
        pass

    def _single_compile(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return
        # import time
        # start = time.time()
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
        # end = time.time()
        # print("took " + str(end - start) + " to compile " + str(obj))

    # imap is evaluated on demand, converting to list() forces execution
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects

def setup_parallel_compile():
    distutils.ccompiler.CCompiler.compile = parallelCCompile

def files_in_dir(directory, ext):
    ret = []
    for file in os.listdir(directory):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.' + ext:
            ret.append(os.path.join(directory, file))
    return ret

def get_tbempy_srces():
    return files_in_dir('cpp', 'cpp')

def get_wrapper_srces():
    return files_in_dir('python_wrapper', 'cpp')

class BoostConfig(object):
    def __init__(self):
        self.root = os.path.join('lib', 'boost')
        self.numpy_root = os.path.join('lib', 'boost_numpy')
        self.include_dirs = [self.root, self.numpy_root]
        self.src_dirs = [
            os.path.join(self.root, 'libs', 'python', 'src'),
            os.path.join(self.numpy_root, 'libs', 'numpy', 'src')
        ]
        self.src_dirs += [os.path.join(self.src_dirs[0], 'converter')]
        self.src_dirs += [os.path.join(self.src_dirs[0], 'object')]
        self.sources = []
        for d in self.src_dirs:
            self.sources.extend(files_in_dir(d, 'cpp'))

def get_blas_flags():
    # Turn off warnings and ignore stdout while we grab the blas/lapack info
    # from numpy
    class RedirectStdStreams(object):
        def __init__(self, stdout=None, stderr=None):
            self._stdout = stdout or sys.stdout
            self._stderr = stderr or sys.stderr

        def __enter__(self):
            self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
            self.old_stdout.flush(); self.old_stderr.flush()
            sys.stdout, sys.stderr = self._stdout, self._stderr

        def __exit__(self, exc_type, exc_value, traceback):
            self._stdout.flush(); self._stderr.flush()
            sys.stdout = self.old_stdout
            sys.stderr = self.old_stderr

    devnull = open(os.devnull, 'w')
    with warnings.catch_warnings():
        with RedirectStdStreams(stdout=devnull, stderr=devnull):
            warnings.simplefilter("ignore")
            blas_lapack_info = np_config_info('lapack_opt', 0)
            return blas_lapack_info

def get_extension_config():
    # -UNDEBUG and -DDEBUG=1 ensure that asserts are turned on
    compile_args = [
        '-std=c++11',
        '-fopenmp',
        '-UNDEBUG',
        '-DDEBUG=1',
        '-O3',
        '-Wall',
        '-Wextra'
    ]
    link_args = [
        '-fopenmp'
    ]
    includes = ['cpp', 'lib']
    return dict(
        sources = get_tbempy_srces(),
        include_dirs = includes,
        extra_compile_args = compile_args,
        extra_link_args = link_args,
        extra_info = get_blas_flags()
    )

def configuration(parent_package='',top_path=None):
    setup_parallel_compile()

    config = Configuration('tbempy', parent_package, top_path)

    # Building the python wrappers requires building with boost.python and
    # boost.numpy, set up the configurations for these.
    boost_config = BoostConfig()

    ext_config = get_extension_config()
    ext_config['sources'] += get_wrapper_srces()
    ext_config['sources'] += boost_config.sources
    ext_config['include_dirs'] += boost_config.include_dirs
    # Warnings that we care about will be caught building the tests. The remaining
    # warnings are in the boost/numpy code and are annoying. Ignore warnings for
    # this build.
    ext_config['extra_compile_args'] += ['-w']

    config.add_extension('_tbempy', **ext_config)
    return config
