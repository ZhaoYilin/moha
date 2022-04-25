.. _installation:

============
Installation
============

Compiling from source code
==========================
This section gives the recommended way to install MoHa. One can find the source code on Github https://github.com/ZhaoYilin/moha.

Prerequisites
-------------
To manually install the package, make sure to have the following dependencies:

* Python >=3.5
* numpy >=1.13.1
* scipy >=0.19.1
* pytest >=6.2.5

To install the prerequisite modules, use the conda or pip command:

* Conda

    .. code-block:: bash

        conda install numpy scipy pytest

* Pip

    .. code-block:: bash

        pip install numpy scipy pytest

Compile
-------
To compile MoHa from the source code, download the package from Github:

    .. code-block:: bash

        gh repo clone ZhaoYilin/moha

Next, move to the :file:`moha` folder and run:

    .. code-block:: bash

        python3 -m pip install -e .

Unit test
---------
To ensure that the installation was successful, run a unit test with the pytest module. Inside the :file:`/moha` folder and run:

    .. code-block:: bash
        
        pytest 

Installation with pip
=====================
The pip package provides a precompiled MoHa code (python wheel) which works on almost all Linux systems and most Mac OS X systems. You can install moha with:

    .. code-block:: bash

        pip install moha

If you already have installed MoHa via pip, you can upgrade it to the new version with:        

    .. code-block:: bash

        pip install --upgrade moha
