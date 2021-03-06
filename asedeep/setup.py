"""Setup asedeep package.
"""

from setuptools import setup

setup(
    name='asedeep',
    version='0.1.0',
    packages=[''],
    url='',
    license='MIT',
    author='Zhenhua Zhang',
    author_email='zhenhua.zhang217@gmail.com',
    description='A tool to generate regulation sequence matrix and allele-specific expression effects from allele-specific read counts.',
    install_requires=['scipy', 'vcf', 'pyfaidx', 'pandas', 'gffutils', 'numpy', 'pytorch']
)
