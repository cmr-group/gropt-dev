import os

import numpy as np
from setuptools import Extension, setup

SRC_DIR = "gropt_dev/src/"

sources = ['gropt_params',
           'ils',
           'ils_cg',
           'op_main',
           'op_bvalue',
           'op_gradient',
           'op_identity',
           'op_moment',
           'op_slew',
           'solver',
           'solver_groptsdmm',
           ]

sourcefiles = ['gropt_dev/gropt_wrapper/gropt_wrapper.pyx'] + [os.path.join(SRC_DIR, f'{x}.cpp') for x in sources]

print(f'{sourcefiles = }')

libraries = []
extra_compile_args = ['-DFMT_UNICODE=0',  # Needed for spdlog
                      '-std=c++11']
extra_link_args = []

include_dirs = ["gropt_dev/src", "gropt_dev/src/external", np.get_include()]
library_dirs = ["gropt_dev/src"]

include_dirs = [os.path.abspath(x) for x in include_dirs]
library_dirs = [os.path.abspath(x) for x in library_dirs]

ext = Extension("gropt_dev.gropt_wrapper",
                sourcefiles,
                language = "c++",
                libraries=libraries,
                include_dirs = include_dirs,
                library_dirs = library_dirs,
                extra_compile_args = extra_compile_args,
                extra_link_args = extra_link_args,
                # undef_macros=['NDEBUG'], # This will *re-enable* the Eigen assertions
            )

setup(
    packages=['gropt_dev', 'gropt_dev.gropt_wrapper'],
    ext_modules=[ext],
)
