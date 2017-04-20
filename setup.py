from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Transform module',
    ext_modules = cythonize("tmpyui/transform.pyx"),
    )
