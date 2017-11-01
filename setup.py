import os
import pysam

from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

compile_flags = [
    '-fgnu89-inline'
]

include_dirs = [
    "coverview"
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
        name="coverview.statistics",
        sources=["coverview/statistics.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    ),
    Extension(
        name="coverview.calculators",
        sources=["coverview/calculators.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    ),
    Extension(
        name="coverview.output",
        sources=["coverview/output.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    ),
    Extension(
        "coverview.reads",
        ["coverview/reads.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=compile_flags,
        libraries=['chtslib'],
        library_dirs=pysam_library_dirs,
        runtime_library_dirs=pysam_library_dirs
    )
]

setup(
    name="CoverView",
    version='1.4.1',
    description="A coverage and quality evaluation tool for targeted and whole exome next-generation sequencing",
    url='https://github.com/RahmanTeam/CoverView',
    author='RahmanTeam',
    author_email='rahmanlab@icr.ac.uk',
    license='MIT',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(modules, compiler_directives=cython_directives),
    packages=[
        'coverview',
        'bamgen',
        'testutils',
        'gui'
    ],
    scripts=[
        "bin/CoverView.py",
        "bin/coverview",
        "test/smoke/check_installation_succeeded.bash",
        "bin/CoverViewGUI.py",
        "bin/gui"
    ],
    zip_safe=False,
    install_requires=[
        "pysam==0.10.0"
    ],
    include_package_data=True
)
