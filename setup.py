from distutils.core import setup

from glob import glob,iglob

import os
import sys
import platform


setup(name='moha',
    version='1.0.0',
    description='A QM package',
    url='https://zhaoyilin.github.io/moha',
    author='Yilin Zhao',
    author_email='zhaoyilin1991@gmail.com',
    license='MIT',
    packages=['moha','moha.system','moha.io','moha.system.integral','moha.system.operator',
        'moha.system.hamiltonian','moha.system.wavefunction','moha.io.basis','moha.hf',
        'moha.property','moha.posthf','moha.posthf.ci','moha.posthf.cc','moha.posthf.pt'],
    package_data={'moha': ['io/basis/*']})
