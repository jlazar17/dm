from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import cython_gsl

cudainc = '/Developer/NVIDIA/CUDA-5.0/include'
cudalib = '/Developer/NVIDIA/CUDA-5.0/lib'

libdir = [cython_gsl.get_library_dir(),cudalib,'../lib']
incdir = [cython_gsl.get_cython_include_dir(),cudainc,'../inc']

libs = cython_gsl.get_libraries()
libs.extend(['m','stdc++','pthread','neutrinocommon','cublas','cudart'])

print libdir, incdir, libs

ext_modules = [Extension("pyneucpp",
                         ["src/pyneucpp.pyx"],
                         language="c++",
                         include_dirs=incdir,
                         libraries=libs,
                         library_dirs=libdir
                         )]

setup(
    name = 'pyneucpp',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )