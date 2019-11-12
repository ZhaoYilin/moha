Tensor Network State
####################
Tensor network states ansatz is used to solve chemistry system with a strong
correlation. Two optimization methods base on tensor network states are
implemented in MoHa, which include variational optimization and time-evolving
block decimation method.

Time-evolving block decimation (TEBD)
======================================
* Consider a Hamiltonian of the form
.. math::
    H = \sum_j h^{[j,j+1]}

* Decompose the Hamiltonian as $H=F+G$
.. math::
    F = \sum_{even} F^{[j]} = \sum_{even j}h^{[j,j+1]}
    G = \sum_{odd} G^{[j]} = \sum_{odd j}h^{[j,j+1]}

* Apply Suzuki-Trotter decomposition, and two chanins of two-site gets
.. math::
    U_{F} = \prod_{even r}exp(-iF^{[r]}\delta t)
    U_{G} = \prod_{odd r}exp(-iG^{[r]}\delta t)

.. figure:: ./pictures/tebd.png
        :scale: 50 %
        :align: center

.. literalinclude:: ../data/examples/tns/tebd_heisenberg.py
        :lines: 1-20
        :caption: /data/examples/tns/tebd_heisenberg.py

Variational Optimization
========================
* Variation principle
.. math::
    \frac{\langle\Psi |H| \Psi\rangle}{\langle\Psi|\Psi\rangle}\geq E_0

* Approach the ground state energy by minimizing
.. math::
    min_{|\Psi\rangle}(\langle\Psi |H| \Psi\rangle - \lambda \langle\Psi|\Psi\rangle)
* Minimize one tensor each time, keeping the other fixed
* sweep over all the tensor several times

.. literalinclude:: ../data/examples/tns/dmrg_heisenberg.py
        :lines: 1-20
        :caption: /data/examples/tns/dmrg_heisenberg.py

Projected Optimization
======================

We can have a trial for solving tensor network state problem by projected approach. 
For simplisity, we start with MPS which is

.. math::
    |MPS\rangle
    =\sum_{\textbf{n}}Tr(U^{\alpha_1}_{m_1}U^{\alpha_2}_{m_1m_2}U^{\alpha_3}_{m_2m_3}\dots)
    |\textbf{n}\rangle

To slove $\hat{H}|MPS\rangle = E|MPS\rangle$, We may set up N equantions to
find out the energy and amplitudes by left projecting a known simple states
$\{\Phi_{\alpha}\}$i,

.. math::
    \langle HF|\hat{H}|MPS\rangle = E\langle HF|MPS\rangle

.. math::
    \langle \Phi_1|\hat{H}|MPS\rangle = E\langle \Phi_1|MPS\rangle

.. math::
    \dots

.. math::
    \langle \Phi_N|\hat{H}|MPS\rangle = E\langle \Phi_N|MPS\rangle

To get the approximate amplitudes, a set of nonlinear equation should be sloved

.. math::
    f_{\alpha}(\textbf{U}) := \langle
    \Phi_{\alpha}|\hat{H}-E|MPS(\textbf{U})\rangle  = 0

Which is usually sloved by quasi-Newton methods,

.. math::
    \textbf{U}^{(\alpha+1)} = \textbf{U}^{(\alpha)} -
    \textbf{F}^{-1}f(\textbf{U}^{(\alpha)}) \ \ \textbf{U}^{0} = 0

The last step of projected MPS is by solving these equantions self consistently, yields approximate amplitudes and an approximate energy. 

TODO
