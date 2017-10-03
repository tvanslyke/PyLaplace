from os import environ
from distutils.core import setup, Extension
module1 = Extension('laplace',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    sources = ['laguerre.cpp',
                               'precomputed_laguerre_roots.cpp',
                               'big_numbers.cpp',
                               'PyLaplace.cpp'],
                    extra_compile_args = ['-std=c++17', '-pipe'])

setup(  name = 'laplace',
        version = '1.0',
        description = 'Provides functionality to obtain the Laplace Transform of arbitrary Python functions.',
        author = 'Timothy Van Slyke',
        author_email = 'vanslyke.t@husky.neu.edu',
        ext_modules = [module1])
