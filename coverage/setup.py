from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

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
        name="coverage.statistics",
        sources=["statistics.pyx"],
        include_dirs=include_dirs,
    ),
    Extension(
        name="coverage.calculators",
        sources=["calculators.pyx"],
        include_dirs=include_dirs,
        libraries=libraries,
        library_dirs=library_dirs,
    ),
    Extension(
        name="coverage.output",
        sources=["output.pyx"],
        include_dirs=include_dirs,
    ),
    Extension(
        "coverage.reads",
        ["reads.pyx"],
        include_dirs=include_dirs,
    )
]

for module in modules:
    module.cython_directives = cython_directives

setup(
    name="CoverView",
    version=1.0,
    description="Calculation of coverage metrics",
    url='https://github.com/RahmanTeamDevelopment/CoverView',
    author='Marton Munz',
    author_email='munzmarci@gmail.com',
    license='MIT',
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(modules)
)


# from setuptools import setup
#
# setup(
#     name='Cava',
#     version='1.2',
#     description='Annotation of genetic variants',
#     url='https://github.com/RahmanTeamDevelopment/CAVA',
#     author='Marton Munz',
#     author_email='munzmarci@gmail.com',
#     license='MIT',
#     packages=['cava'],
#     scripts=[
#         'bin/cava.py',
#         'bin/dbsnp_prep.py',
#         'bin/ensembl_prep.py'
#     ],
#     zip_safe=False
# )
