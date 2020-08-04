from setuptools import setup

setup(
    name='asedlp',
    version='0.1.0',
    packages=[''],
    url='',
    license='MIT',
    author='Zhenhua Zhang',
    author_email='zhenhua.zhang217@gmail.com',
    description='A tool to create regulation sequence matrix and ASE effects from outputs by WASP.',
    install_requires=['scipy', 'tables', 'gffutils', 'numpy', 'pytorch',
                      'biopython']
)
