from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext


setup(
    name='GridShiftPP2',
    version='1.0',
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("GridShiftPP2",
                 sources=["gridshiftpp2.pyx"],
                 language="c++",
                 include_dirs=[numpy.get_include()])],
    author='anomyous',
    author_email='anomyous@gmail.com'

)
