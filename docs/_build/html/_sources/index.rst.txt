.. hftest documentation master file

hftest — A Hartree-Fock Program in C++23
=========================================

**hftest** is a self-contained Restricted and Unrestricted Hartree-Fock (RHF/UHF)
program written in modern C++23.  It uses `LibTorch <https://pytorch.org/cppdocs/>`_
for all linear algebra, reads molecular geometry and settings from a
`TOML <https://toml.io>`_ input file, and loads basis sets from Basis Set Exchange
(BSE) Gaussian-format files.

.. code-block:: none

   H₂O / STO-3G  ·  RHF

   ╔══════════════════════════════════╗
   ║  Iteration   Energy (Ha)  ΔE        ΔRMS-P  ║
   ╠══════════════════════════════════╣
   ║       1  -74.9420799  —           —        ║
   ║      ...                                   ║
   ║      12  -75.0171735  5.8e-10  3.1e-09  ║
   ╚══════════════════════════════════╝
   Total energy : -75.0171735279 Ha
   HOMO : -0.3971 Ha  LUMO : 0.5758 Ha  gap : 0.9729 Ha

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   usage

.. toctree::
   :maxdepth: 2
   :caption: Theory & Algorithms

   principles
   algorithms

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   api/molecule
   api/basis
   api/boys
   api/integrals
   api/scf
   api/output

Indices
-------

* :ref:`genindex`
* :ref:`search`
