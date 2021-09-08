.. moha documentation master file, created by
   sphinx-quickstart on Fri Apr 20 03:14:34 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MoHa |version|
##############
MoHa is the abbreviation of Molecular/Model Hamiltonian. It is a general-purpose electronic structure program written in python 
that emphasizes code simplicity and facilitates new method development. The package provides the standard module including 
custom Hamiltonians, using mean-field and post-mean-field methods with standard Gaussian basis functions，post process for
varieties properties. Besides the modules for a canonical module, MoHa also provides model Hamiltonian and Tensor network
states method.

Contents
========

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: INSTALL:

   install_linux.rst
   install_mac.rst

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: USER DOCUMENTATION:

   user_io.rst
   user_molecular_system.rst
   user_hf.rst
   user_posthf.rst
   user_property.rst

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: DEVELOPER DOCUMENTATION:

   developer_api.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

