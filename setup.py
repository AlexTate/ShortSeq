import os

from setuptools import setup, find_namespace_packages, Extension
from Cython.Build import cythonize

# Package metadata
NAME = 'shortseq'
DESCRIPTION = 'Compact representation of short DNA sequences that uses up to 73% less memory'
URL = 'https://github.com/AlexTate/shortseq'
AUTHOR = 'Alex Tate'
PLATFORM = 'Unix'
REQUIRES_PYTHON = '>3.8, <3.12'
VERSION = '0.0.1'

define_macros = []
short_seq_common_compile_args = [
    '-std=c++20',
    "-O3",
    "-mbmi2",
    "-mpopcnt",
    "-mtune=native",
    # '-march=native',
]

cython_implementations = [
    "shortseq/short_seq.pyx",
    "shortseq/short_seq_var.pyx",
    "shortseq/short_seq_128.pyx",
    "shortseq/short_seq_64.pyx",
    "shortseq/fast_read.pyx",
    "shortseq/counter.pyx",
    "shortseq/util.pyx",
    "shortseq/umi/umi.pyx",
]

if os.environ.get("CYTRACE") == '1':
    define_macros.append(('CYTHON_TRACE_NOGIL', '1'))

extensions = [
    Extension(
        pyx.replace("/", ".").replace(".pyx", ""),
        sources=[pyx],
        extra_compile_args=short_seq_common_compile_args,
        define_macros=define_macros,
        language='c++',
    )
    for pyx in cython_implementations
]

setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    description=DESCRIPTION,
    include_package_data=True,
    packages=find_namespace_packages(),
    python_requires=REQUIRES_PYTHON,
    zip_safe=False,
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': '3'},
        build_dir="build"
    ),
)