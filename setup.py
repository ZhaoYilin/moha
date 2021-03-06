from distutils.core import setup

from glob import glob

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
    packages=['moha','moha.system','moha.io','moha.system.integral','moha.system.operator','moha.system.hamiltonian','moha.system.wavefunction','moha.io.basis','moha.hf','moha.property','moha.posthf','moha.posthf.ci','moha.posthf.cc','moha.posthf.pt','moha.modelsystem','moha.tn','moha.tools','moha.tools.orbital_order','moha.tools.orbital_entanglement'],
    package_data={'moha': ['io/basis/*']},
    data_files=[
        ('share/moha', glob('data/*.*')),
        ('share/moha/example', glob('data/example/*.*'))])
