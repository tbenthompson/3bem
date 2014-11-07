from fabricate import *

import os

def files_in_dir(directory, ext):
    ret = []
    for file in os.listdir(directory):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.' + ext:
            ret.append(os.path.join(directory, file_name))
    return ret

dirs = ['src', 'examples', 'test', 'inttest']
lib_srces = files_in_dir("src", "cpp")
examples = files_in_dir("examples", "cpp")
tests = files_in_dir("test", "cpp")
inttests = files_in_dir("inttest", "cpp")

compiler = 'mpic++'

petsc_dir = os.environ['PETSC_DIR']
petsc_arch = os.environ['PETSC_ARCH']
includes = [
    './src',
    './lib/',
    './lib/unittest-cpp/src',
    './lib/actor-framework/libcaf_core',
    './lib/actor-framework/libcaf_opencl',
    './lib/autocheck/include',
    petsc_dir + '/' + petsc_arch + '/include',
    petsc_dir + '/include'
]

cpp_flags = '-Wall -std=c++11 -fopenmp -mavx'.split()
cpp_flags.extend(['-I' + loc for loc in includes])

debug_flags = '-g -Og -DDEBUG=1'.split()
release_flags = '-DNDEBUG=1 -O3 -funroll-loops'.split()
profile_flags = release_flags + ['-g']
# cpp_flags.extend(debug_flags)
# cpp_flags.extend(release_flags)
cpp_flags.extend(profile_flags)

lib_cpp_flags = ['-fPIC']
lib_cpp_flags.extend(cpp_flags)

link_flags = '-fopenmp -lOpenCL -lhdf5'.split()
link_flags.append('-Wl,-rpath=./lib/actor-framework/build/lib')
link_flags.append('-L./lib/actor-framework/build/lib')
link_flags.append('-lcaf_core')
link_flags.append('-lcaf_opencl')
link_flags.append('-Wl,-rpath=' + petsc_dir + '/' + petsc_arch + '/lib')
link_flags.append('-L' + petsc_dir + '/' + petsc_arch + '/lib')
link_flags.append('-lpetsc')

lib_link_flags = ['-shared']
lib_link_flags.extend(link_flags)

lib_dep_flags = ['-Wl,-rpath=./build', '-L./build', '-l3bem']

test_link_flags = ['-L./lib/unittest-cpp']
test_link_flags.append('-lUnitTest++')
test_link_flags.extend(lib_dep_flags)
test_link_flags.extend(link_flags)

example_link_flags = lib_dep_flags
example_link_flags.extend(link_flags)
example_link_flags.extend('-lglut -lGL -lGLU -lGLEW'.split())

build_dir = 'build'

def setup_dir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
def build():
    for d in dirs:
        setup_dir(os.path.join(build_dir, d))
    compile()
    link()

def oname(filename):
    return os.path.join(build_dir, filename)

def compile():
    for source in lib_srces:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), lib_cpp_flags)
    for source in examples:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), cpp_flags)
    for source in tests:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), cpp_flags)
    for source in inttests:
        run(compiler, '-c', source + '.cpp', '-o', oname(source + '.o'), cpp_flags)

def link():
    lib_objs = [oname(s + '.o') for s in lib_srces]
    after()
    run(compiler, '-o', oname('lib3bem.so'), lib_objs, lib_link_flags)
    after()
    for source in examples:
        run(compiler, '-o', oname(source), oname(source + '.o'), example_link_flags)
    for source in tests:
        run(compiler, '-o', oname(source), oname(source + '.o'), test_link_flags)
    for source in inttests:
        run(compiler, '-o', oname(source), oname(source + '.o'), test_link_flags)

def clean():
    autoclean()

def rebuild():
    clean()
    build()

main(parallel_ok = True, jobs = 12)
