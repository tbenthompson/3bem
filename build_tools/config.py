from __future__ import print_function
from numpy.distutils.system_info import get_info as np_config_info
import distutils.sysconfig
import warnings
import os

# Code to get the compiler distutils will use
def get_compiler():
    import distutils.sysconfig
    import distutils.ccompiler
    compiler = distutils.ccompiler.new_compiler()
    distutils.sysconfig.customize_compiler(compiler)
    return compiler.compiler_cxx

def files_in_dir(directory, ext):
    ret = []
    for file in os.listdir(directory):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.' + ext:
            ret.append(os.path.join(directory, file_name))
    return ret

def include(dirs):
    return ['-I' + d for d in dirs]

src_dir = 'cpp'
test_dir = 'unit_tests'
py_wrap_dir = 'python_wrapper'
tbempy_dir = 'tbempy'

lib_srces = files_in_dir(src_dir, 'cpp')
wrapper_srces = files_in_dir(py_wrap_dir, 'cpp')
test_srces = files_in_dir(test_dir, 'cpp')

include_dirs = [str(src_dir), 'lib', os.path.join('lib', 'eigen')]
include_flags = include(include_dirs)

base_cpp_flags = ['-Wextra', '-std=c++11', '-fopenmp', '-DDEBUG=1', '-fPIC']
cpp_flag_types = dict()
cpp_flag_types['debug'] = ['-g', '-Og']
cpp_flag_types['release'] = ['-Ofast']
cpp_flag_types['profile'] = cpp_flag_types['release'] + ['-g']

base_link_flags = ['-fopenmp']
link_flag_types = dict()

boost_root = os.path.join('lib', 'boost')
boost_numpy_root = os.path.join('lib', 'boost_numpy')
boost_include_dirs = [boost_root, boost_numpy_root]
boost_src_dirs = [
    os.path.join(boost_root, 'libs', 'python', 'src'),
    os.path.join(boost_numpy_root, 'libs', 'numpy', 'src')
]
boost_src_dirs += [os.path.join(boost_src_dirs[0], 'converter')]
boost_src_dirs += [os.path.join(boost_src_dirs[0], 'object')]
boost_sources = reduce(lambda fs1, fs2: fs1 + fs2,
                       map(lambda d: files_in_dir(d, 'cpp'), boost_src_dirs))
python_include_dir = distutils.sysconfig.get_python_inc()
python_lib_dir = distutils.sysconfig.get_python_lib(standard_lib = True)
python_lib = distutils.sysconfig.get_config_var('LDLIBRARY')

def determine_build_type(command_params):
    build_type = 'release'
    if '-d' in command_params:
        build_type = 'debug'
    if '-p' in command_params:
        build_type = 'profile'
    return build_type

def determine_printer(command_params):
    printer = lambda x: None
    if '-v' in command_params:
        printer = print
    return printer

def get_config(command_params):
    build_type = determine_build_type(command_params);
    build_dir = os.path.join('build', str(build_type))
    cpp_flags = base_cpp_flags + include_flags + cpp_flag_types[build_type]
    link_flags = base_link_flags + link_flag_types.get(build_type, [])

    printer = determine_printer(command_params)

    tests = dict()
    tests['cpp_flags'] = cpp_flags
    tests['link_flags'] = link_flags
    tests['sources'] = test_srces + lib_srces
    tests['linked_sources'] = []
    tests['binary_name'] = 'test_runner'
    tests['priority'] = 1

    python_wrapper_cpp_flags = cpp_flags +\
            include([python_include_dir] + boost_include_dirs)
    python_wrapper_link_flags = link_flags + ['-shared'] +\
        ['-L' + python_lib_dir] +\
        ['-l:' + python_lib]
    python_wrapper = dict()
    python_wrapper['cpp_flags'] = python_wrapper_cpp_flags
    python_wrapper['link_flags'] = python_wrapper_link_flags
    python_wrapper['sources'] = wrapper_srces + boost_sources
    python_wrapper['linked_sources'] = lib_srces
    python_wrapper['binary_name'] = 'tbempy.so'
    python_wrapper['priority'] = 2

    c = dict()
    c['build_dir'] = build_dir
    c['subdirs'] = dict(
        src_dir = src_dir,
        test_dir = test_dir,
        py_wrap_dir = py_wrap_dir,
        tbempy_dir = tbempy_dir
    )
    c['subdirs'].update({d:d for d in boost_src_dirs})
    c['compiler'] = get_compiler()
    c['command_params'] = command_params
    c['printer'] = printer
    c['targets'] = dict()
    c['targets']['python_wrapper'] = python_wrapper
    c['targets']['tests'] = tests
    return c

