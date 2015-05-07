import numpy.distutils.command.build_ext as _build_ext

class tbempyBuildExt(_build_ext.build_ext):
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
        def remove_arch(flags):
            new_flags = []
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
        self._cxx_compiler.compiler_so = remove_arch(
            self._cxx_compiler.compiler_so
        )
        self._cxx_compiler.linker_so = remove_arch(
            self._cxx_compiler.linker_so
        )

        _build_ext.build_ext.build_extension(self, ext)
