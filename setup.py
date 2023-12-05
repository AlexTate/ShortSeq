from setuptools import setup, find_namespace_packages, Extension
from Cython.Build import cythonize

# Package metadata
NAME = 'ShortSeq'
DESCRIPTION = 'Compact representation of short DNA sequences that uses up to 73% less memory'
URL = 'https://github.com/AlexTate/ShortSeq'
AUTHOR = 'Alex Tate'
PLATFORM = 'Unix'
REQUIRES_PYTHON = '>=3.8, <=3.11'
VERSION = '0.0.1'


short_seq_common_compile_args = [
    '-std=c++20',
    "-O3",
    '-march=native']

extensions = [
    Extension("ShortSeq.short_seq",
              sources=['ShortSeq/short_seq.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("ShortSeq.short_seq_var",
              sources=['ShortSeq/short_seq_var.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("ShortSeq.short_seq_128",
              sources=['ShortSeq/short_seq_128.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("ShortSeq.short_seq_64",
              sources=['ShortSeq/short_seq_64.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("ShortSeq.util",
              sources=['ShortSeq/util.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("ShortSeq.fast_read",
              sources=['ShortSeq/fast_read.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
    Extension("ShortSeq.umi.umi",
              sources=['ShortSeq/umi/umi.pyx'],
              extra_compile_args=short_seq_common_compile_args,
              language='c++'),
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