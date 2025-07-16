import os

import numpy as np
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CustomBuildExt(build_ext):
    """
    Custom build_ext command to handle C++ extensions with specific compiler flags.
    
    This class extends the build_ext command to set extra compile arguments based on the compiler type.
    """
    
    def build_extension(self, ext):
        compiler_type = self.compiler.compiler_type
        
        ext.extra_compile_args = ['-DFMT_UNICODE=0',  # Needed for spdlog
                                  ]
        if compiler_type == 'msvc':
            ext.extra_compile_args += []
        else:
            ext.extra_compile_args += ['-std=c++11']

        super().build_extension(ext)


SRC_DIR = "gropt_dev/src/"

sources = ['gropt_params',
           'gropt_utils',
           'ils',
           'ils_cg',
           'ils_nlcg',
           'ils_bicgstabl',
           'op_main',
           'op_bvalue',
           'op_gradient',
           'op_identity',
           'op_moment',
           'op_safe',
           'op_slew',
           'op_tv',
           'solver',
           'solver_groptsdmm',
           ]

sourcefiles = ['gropt_dev/gropt_wrapper/gropt_wrapper.pyx'] + [os.path.join(SRC_DIR, f'{x}.cpp') for x in sources]

include_dirs = ["gropt_dev/src", "gropt_dev/src/external", np.get_include()]
library_dirs = ["gropt_dev/src"]

include_dirs = [os.path.abspath(x) for x in include_dirs]
library_dirs = [os.path.abspath(x) for x in library_dirs]

ext = Extension("gropt_dev.gropt_wrapper",
                sourcefiles,
                language = "c++",
                include_dirs = include_dirs,
                library_dirs = library_dirs,
                # libraries=libraries,   # Handled in CustomBuildExt now
                # extra_compile_args = extra_compile_args,
                # extra_link_args = extra_link_args,
                # undef_macros=['NDEBUG'], # This will *re-enable* the Eigen assertions, COMMENT THIS FOR RELEASE
            )

setup(
    packages=['gropt_dev', 'gropt_dev.gropt_wrapper'],
    ext_modules=[ext],
    cmdclass={'build_ext': CustomBuildExt},
)
