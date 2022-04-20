.. _installation:

============
Installation
============

Compiling from source code
==========================
This the the recommended way to install moha. One can find the source code on Github https://github.com/ZhaoYilin/moha.

Prerequisites
-------------
To manual install the package, make sure you have the following dependencies:

* Python >=3.5
* numpy >=1.13.1
* scipy >=0.19.1
* pytest >=6.2.5

These prerequistied modules can be installed with conda or pip:

* Conda

    .. code-block:: bash

        conda install numpy scipy pytest

* Pip

    .. code-block:: bash

        pip install numpy scipy pytest

Compile
-------
The moha package can be downloaded on Github by:

    .. code-block:: bash

        gh repo clone ZhaoYilin/moha

To install the moha package, move to the moha folder and run:

    .. code-block:: bash

        python3 -m pip install -e .

Unit test
---------
To ensure that the installation was successful, you can run an unit test with pytest module. Inside the moha folder and run:

    .. code-block:: bash
        
        pytest 

Installation with pip
=====================

The pip package provides a precompiled moha code (python wheel) which works on almost all Linux systems, and most of Mac OS X systems. you can install moha with:        

    .. code-block:: bash

        pip install moha

If you already have installed moha via pip, you can upgrade it to the new version with:        

    .. code-block:: bash

        pip install --upgrade moha
