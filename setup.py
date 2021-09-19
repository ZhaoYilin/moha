from setuptools import setup,find_packages

import os
import sys
import platform


setup(name='moha',
    version='0.1',
    description='A QM package',
    url='http://github.com/fhqgfss/penta',
    author='Yilin Zhao',
    author_email='zhaoyilin1991@gmail.com',
    license='MIT',
    packages=find_packages(),
    package_data={'moha': ['io/basis_database/*']})
