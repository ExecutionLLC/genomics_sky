from setuptools import Extension
from setuptools import setup

from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Distutils import Extension


# all project cython-extensions
pyx_modules = (
    ('utils', ['utils.pyx']),
)


# extra options for gcc
extra_compile_args = ['-O3', '-finline-functions']

# create extensions from cython-modules
pyx_extensions = []

for ext_name, ext_source in pyx_modules:
    pyx_extensions.append(
        Extension(
            name=ext_name,
            sources=ext_source,
            extra_compile_args=extra_compile_args
        )
    )


# cython compiler directives, see:
# http://docs.cython.org/src/reference/compilation.html#compiler-directives
compiler_directives = dict(
    boundscheck=False,
    wraparound=False
)


setup(
    name='DraftUtils',
    version='0.2',
    description='Set of various useful utilites',
    author='akonovalov',
    author_email='antonkonovalov1976@gmail.com',
    install_requires=[
        'Cython',
        'numpy',
    ],
    # for generate html-files for cythonize reporting - set annotate=True
    ext_modules=cythonize(
        pyx_extensions,
        compiler_directives=compiler_directives,
        annotate=False
    ),
    cmdclass={
        'build_ext': build_ext
    },
    packages=[
        'utils'
    ]

)
