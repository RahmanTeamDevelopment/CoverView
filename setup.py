import os
import pysam

from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

compile_flags = [
    '-fgnu89-inline'
]

include_dirs = [
    "coverview_"
]

include_dirs.extend(
    pysam.get_include()
)

cython_directives = {
    "boundscheck": False,
    "nonecheck": False,
    "cdivision": True,
    "profile": False,
    "initializedcheck": False,
    "wraparound": True
}

pysam_library_dirs = list(set(os.path.dirname(x) for x in pysam.get_libraries()))

modules = [
    Extension(
        name="coverview_.statistics",
        sources=["coverview_/statistics.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    ),
    Extension(
        name="coverview_.calculators",
        sources=["coverview_/calculators.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    ),
    Extension(
        name="coverview_.output",
        sources=["coverview_/output.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    ),
    Extension(
        "coverview_.reads",
        ["coverview_/reads.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    )
]

setup(
    name="CoverView",
    version='1.4.3',
    description="A coverage and quality evaluation tool for targeted and whole exome next-generation sequencing",
    url='https://github.com/RahmanTeam/CoverView',
    author='RahmanTeam',
    author_email='rahmanlab@icr.ac.uk',
    license='MIT',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(modules, compiler_directives=cython_directives),
    packages=[
        'coverview_',
        'bamgen',
        'gui_',
        'tgmi',
        'ensembldb'
    ],
    scripts=[
        "bin/CoverView.py",
        "bin/coverview",
        "bin/CoverViewGUI.py",
        "bin/gui",
        "bin/EnsemblDB.py",
        "bin/ensembl_db"
    ],
    zip_safe=False,
    install_requires=[
        "pysam==0.10.0"
    ],
    include_package_data=True
)
