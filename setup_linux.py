from setuptools import setup, Extension
import glob

with open("README.md", "r") as fp:
    long_desc = fp.read()

source_files = glob.glob("./src/*.c")

setup(
    name='csdpy',
    version='0.1',
    author="Levente Bodnar",
    author_email="bodnalev@gmail.com",
    description="Python wrapper for CSDP",
    long_description_content_type='text/markdown',
    long_description=long_desc,
    license="EPL",
    python_requires='>=3.6',
    ext_modules = [Extension(
        "csdpy", 
        source_files, 
        include_dirs=["./include"], 
        extra_compile_args=[
            "-m64", "-march=native", "-mtune=native", "-O3 -ffast-math", "-fPIC",
            "-fopenmp", "-Wall", "-DBIT64", "-DUSEOPENMP", "-std=c99",
            "-DSETNUMTHREADS", "-DUSESIGTERM", "-DUSEGETTIME", "-I../include", "-g"
        ],
        extra_link_args=[
            "-L/opt/homebrew/opt/libomp/lib",
            "-L../src", "-llapack", "-lblas", "-lm", "-lgomp"
        ]
    )]
)