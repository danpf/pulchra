# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup, find_packages

__version__ = "1.0.0"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
ext_modules = [
    Pybind11Extension(
        "pypulchra._pypulchra",
        ["src/pulchra.cpp", "bindings/python/pulchra_bindings.cpp", "bindings/python/pypulchra_bindings.cpp"],
        include_dirs=["src"],
        define_macros=[("VERSION_INFO", __version__)],
        extra_compile_args=["-static"],
        cxx_std=14,
    )
]


setup(
    name="pypulchra",
    version=__version__,
    author="Danny Farrell",
    author_email="danpfuw@gmail.com",
    description="Python bindings for the pulchra CA -> full bb trace library",
    url="https://github.com/danpf/pulchra",
    long_description="",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11"],
    ext_modules=ext_modules,
    zip_safe=False,
    # extras_require={"test": ["pytest>=6.0"]},
    python_requires=">=3.7",
    packages=["pypulchra"],
    entry_points={
        "console_scripts": [
            "pypulchra = pypulchra:_commandline",
        ],
    },
)
