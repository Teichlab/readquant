from glob import glob
from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(
        name='readquant',
	version='1.2.0',
        description='Convenience package for parsing RNA-seq quantification results',
        long_description=readme(),
        packages=find_packages(),
        install_requires=['pandas', 'tqdm'],
        scripts=glob('scripts/*.py'),
        author='Valentine Svensson',
        author_email='valentine@nxn.se',
        license='MIT'
    )
