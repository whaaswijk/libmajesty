from distutils.core import setup, Extension
from Cython.Build import cythonize
import platform
import pkg_resources
import os
from sys import platform

is_mac = platform == 'darwin'

INCLUDE_DIRS = [pkg_resources.resource_filename('numpy', 'core/include'),
                '../include', 
                '../build/abc/src/abc-project/src',
                '../build/hiredis/src/hiredis-project',
                #os.environ['HIREDIS_HOME'],
                os.environ['BOOST_ROOT']+'/include']
LIBRARY_DIRS = []
LIBRARIES = []
EXTRA_OBJECTS = ['../build/src/libmajesty.a',
               '../build/src/minisat/libMiniSat.a',
               '../build/hiredis/src/hiredis-project/libhiredis.a',
               '../build/abc/src/abc-project/libabc.a',
               os.environ['BOOST_ROOT']+ ('/lib/libboost_filesystem.a' if is_mac else '/lib/libboost_filesystem.so')]
EXTRA_COMPILE_ARGS = ["-std=c++11"]
EXTRA_LINK_ARGS = ["-std=c++11"]

if is_mac:
    EXTRA_COMPILE_ARGS += ['-mmacosx-version-min=10.9']

setup(ext_modules=cythonize(Extension(
    'majespy',
    sources=['majespy.pyx'],
    include_dirs=INCLUDE_DIRS,
    library_dirs=LIBRARY_DIRS,
    libraries=LIBRARIES,
    extra_objects=EXTRA_OBJECTS,
    language='c++',
    extra_compile_args=EXTRA_COMPILE_ARGS,
    extra_link_args=EXTRA_LINK_ARGS
)))
