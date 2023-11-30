from setuptools import setup, Extension
from Cython.Build import cythonize

# Package metadata
NAME = 'ShortSeq'
DESCRIPTION = ''
URL = 'https://github.com/AlexTate/ShortSeq'
EMAIL = '0xalextate@gmail.com'
AUTHOR = 'Alex Tate'
PLATFORM = 'Unix'
REQUIRES_PYTHON = '>=3.10.0'
VERSION = '0.1'


short_seq_common_compile_args = [
    '-stdlib=libc++',
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
    Extension("ShortSeq.short_seq_util",
              sources=['ShortSeq/short_seq_util.pyx'],
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
    author_email=EMAIL,
    description=DESCRIPTION,
    include_package_data=True,
    packages=['ShortSeq'],
    zip_safe=False,
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': '3'},
        build_dir="build"
    ),
)