from distutils.core import setup, Extension
from Cython.Build import cythonize
import platform
import pkg_resources

INCLUDE_DIRS = [pkg_resources.resource_filename('numpy', 'core/include'), '../include']
LIBRARY_DIRS = []
LIBRARIES = []
EXTRA_OBJECTS = ["../build/src/liblibmajesty.a",
               "../build/src/minisat/libMiniSat.a"]
EXTRA_COMPILE_ARGS = ["-std=c++11"]
EXTRA_LINK_ARGS = ["-std=c++11"]

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
