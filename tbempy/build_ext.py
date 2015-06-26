from __future__ import print_function
import numpy.distutils.command.build_ext as _build_ext
import subprocess
import sys

def remove_arch(flags):
    new_flags = []
    # One parameter after '-arch' must also be removed since arch params
    # are specified like '-arch i386'
    skip = 0
    for i in range(len(flags)):
        if skip > 0:
            skip -= 1
            continue
        p = flags[i]
        if p == '-arch':
            skip = 1
            continue
        new_flags.append(p)
    return new_flags

def check_compiler(compiler):
    check_compiler_version_cmd = [compiler, '--version']
    p = subprocess.Popen(
        check_compiler_version_cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    p.wait()
    out = p.stdout.read().lower().decode('utf-8')
    if 'llvm' in out:
        return 'clang++'
    if 'Free Software Foundation' in out:
        return 'g++'
    if 'Intel' in out:
        return 'icpc'
    return

class tbempyBuildExt(_build_ext.build_ext):
    def warn_and_remove_openmp_if_clang(self, compiler, ext):
        if compiler != 'clang++':
            return

        print('')
        print('')
        print('')
        print('You are compiling with clang. clang does not yet support OpenMP '
              ' so OpenMP will be turned off. \nTo use a different compiler, try: '
              '\nCXX=compilername python ' + ' '.join(sys.argv))
        print('Continue (y/n):', end='')
        # response = raw_input()
        # if (response.lower() == 'n'):
        #     sys.exit()
        ext.extra_compile_args.remove('-fopenmp')
        ext.extra_link_args.remove('-fopenmp')

    def ignore_unused_arguments_if_clang(self, compiler, ext):
        if compiler != 'clang++':
            return

        ext.extra_compile_args.append('-Qunused-arguments')
        ext.extra_link_args.append('-Qunused-arguments')

    def build_extension(self, ext):
        # Some warning parameters that are set in the distutils.sysconfig
        # for clang on Mac OS X do not work on gcc. A simple solution is to
        # just remove all warning parameters and add back in my own later.
        self._cxx_compiler.compiler_so = [
            p for p in self._cxx_compiler.compiler_so
            if not p.startswith('-W')
        ]

        # Remove any flags specifying architecture. This makes cross compiling
        # impossible, but makes compiling for the current architecture work
        # with both clang and gcc
        self._cxx_compiler.compiler_so = remove_arch(
            self._cxx_compiler.compiler_so
        )
        self._cxx_compiler.linker_so = remove_arch(
            self._cxx_compiler.linker_so
        )

        compiler = check_compiler(self._cxx_compiler.compiler_so[0])
        self.warn_and_remove_openmp_if_clang(compiler, ext)
        self.ignore_unused_arguments_if_clang(compiler, ext)
        _build_ext.build_ext.build_extension(self, ext)
