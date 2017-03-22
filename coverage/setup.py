from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

include_dirs = [
    "../env/lib/python2.7/site-packages",
    "../env/lib/python2.7/site-packages/pysam",
    "../env/lib/python2.7/site-packages/pysam/include/htslib",
    "../env/lib/python2.7/site-packages/pysam/include/htslib/htslib",
    "."
]

library_dirs = [
    "/usr/local/lib"
]

libraries = [
    "hts"
]

cython_directives = {
    "boundscheck": False,
    "nonecheck" : False,
    "cdivision" : True,
    "profile" : False,
    "initializedcheck" : False,
    "wraparound" : True
}

modules = [
    Extension(
        "statistics",
        ["statistics.pyx"],
        include_dirs=include_dirs,
    ),
    Extension(
        "coverage",
        ["calculators.pyx"],
        include_dirs=include_dirs,
        libraries=libraries,
        library_dirs=library_dirs,
    )
]

for module in modules:
    module.cython_directives = cython_directives

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = modules
)
