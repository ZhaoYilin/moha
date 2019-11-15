Hartree-Fock
############
The SCF procedure is the workhorse of most computational chemistry software packages, as it is generally the first step before doing more advanced calculations (such as MP2, etc). 

Plain Solver
============
Starting from the basics, we'll run a simple calculation for a water molecule with a bond-length of 1.1 Å and a bond angle of 104.0o with an STO-3G basis set. 

* Restricted Hartree-Fock
    .. literalinclude:: ../data/examples/scf/rhf_water_sto3g.py
            :lines: 1-20
            :caption: /data/examples/scf/rhf_water_sto3g.py

* Unrestricted Hartree-Fock
    .. literalinclude:: ../data/examples/scf/uhf_water_sto3g.py
            :lines: 1-20
            :caption: /data/examples/scf/uhf_water_sto3g.py

.. figure:: ./pictures/h2.png
        :scale: 100%
        :align: center

DIIS Solver
============
ToDo
