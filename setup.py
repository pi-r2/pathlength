from setuptools import setup
from Cython.Build import cythonize
import sys
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

"""
directives = {
    'optimize.inline_defnode_calls': True,
    'language_level': sys.version_info[0] # "2" or "3"
}


"""

setup(
    setup_requires=["cysignals"],
    ext_modules = cythonize("pathlength.pyx", annotate=True)
    )


