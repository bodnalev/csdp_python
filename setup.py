from distutils.command.sdist import sdist
from distutils.errors import DistutilsExecError

from setuptools import setup  


class sdist_make(sdist):
	def run(self):
		try:
			self.spawn(['./csdp/Makefile'])
		except DistutilsExecError:
			self.warn('Failed to create the binaries')
		super().run()

with open("README.md", "r") as fp:
	long_desc = fp.read()

setup(
	name='csdpy',
	version='0.1',
	author="Levente Bodnar",
	author_email="bodnalev@gmail.com",
	tests_require=["pytest"],
	description="Python wrapper for CSDP",
	long_description_content_type='text/markdown',
	long_description=long_desc,
	license="EPL",
	py_modules=['csdpy'],
    python_requires='>=3.6', 
	cmdclass={
		'sdist': sdist_make
	}
)