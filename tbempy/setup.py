import os
import multiprocessing.pool
import multiprocessing
import distutils.ccompiler
from numpy.distutils.misc_util import Configuration

# Patch for parallel compilation with distutils
# From: http://stackoverflow.com/questions/11013851/speeding-up-build-process-with-distutils
def parallelCCompile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N = 1
    try:
        N = multiprocessing.cpu_count()
    except e:
        pass
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects

def setup_parallel_compile():
    # Set distutils to use the parallel compiler from above
    distutils.ccompiler.CCompiler.compile=parallelCCompile

def files_in_dir(directory, ext):
    ret = []
    for file in os.listdir(directory):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.' + ext:
            ret.append(os.path.join(directory, file_name))
    return ret

def get_tbempy_srces():
    return files_in_dir('cpp', 'cpp')

def get_wrapper_srces():
    return files_in_dir('python_wrapper', 'cpp')

def get_test_srces():
    return files_in_dir('unit_tests', 'cpp')

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

compile_args = ['-std=c++11', '-fopenmp', '-DDEBUG=1']
link_args = ['-fopenmp']

includes = ['cpp', 'lib', os.path.join('lib', 'eigen')]

def configuration(parent_package='',top_path=None):

    setup_parallel_compile()

    config = Configuration('tbempy', parent_package, top_path)

    boost_config = BoostConfig()
    lib_includes = includes + boost_config.include_dirs
    sources = get_tbempy_srces() + get_wrapper_srces() + boost_config.sources
    sources = [s + '.cpp' for s in sources]
    config.add_extension(
        'tbempy',
        include_dirs = lib_includes,
        sources = sources,
        extra_compile_args = compile_args,
        extra_link_args = link_args
    )
    return config
