.. -*- rst -*-

Asynchronous Sigma-Delta Modulator Machines
===========================================

Encoding with Asynchronous Sigma-Delta Modulator (ASDM)
-------------------------------------------------------
.. index:: asdm_encode.m

* **ASDM Encoding** (``asdm_encode.m``) |lazar_perfect_2004|_:

  Encodes a bandlimited signal using an ASDM encoder. A signal can be
  encoded by a single ASDM encoder (Single-Input Single-Output
  Encoding), as shown below, or it can be encoded by multiple of such
  ASDM encoders (Single-Input Multiple-Output Encoding).

   .. image:: images/tem-asdm.png
      :scale: 60
      :align: center

Decoding for Signal Encoded with Single-Input Single-Output ASDM Encoder
------------------------------------------------------------------------
.. index:: asdm_decode.m

* **Decoding** (``asdm_decode.m``) |lazar_perfect_2004|_:

  Reconstructs a bandlimited signal encoded by an Single-Input
  Single-Output ASDM encoder using sinc kernels.

   .. image:: images/tdm-sinc.png
      :scale: 60
      :align: center

.. index:: asdm_decode_fast.m

* **Decoding using Fast Approximation Method** (``asdm_decode_fast.m``) |lazar_fast_2005|_:

  Reconstructs a bandlimited signal encoded by a Single-Input
  Single-Output ASDM encoder using a fast approximation method.

   .. image:: images/tdm-fast.png
      :scale: 60
      :align: center

.. index:: asdm_decode_ins.m

* **Decoding using Threshold-Insensitive Method**  (``asdm_decode_ins.m``) |lazar_perfect_2004|_:

  Reconstructs a bandlimited signal encoded by a Single-Input
  Single-Output ASDM encoder. This reconstruction method does not
  require the specification of an integrator threshold.

   .. image:: images/tdm-sinc-ins.png
      :scale: 60
      :align: center

Decoding for Signal Encoded with Single-Input Multiple-Output ASDM Encoder
--------------------------------------------------------------------------
.. index:: asdm_decode_pop.m

* **Population Decoding** (``asdm_decode_pop.m``) |lazar_information_2007|_:

  Reconstructs a bandlimited signal encoded by multiple ASDM encoders
  using sinc kernels.

   .. image:: images/tdm-sinc-miso.png
      :scale: 60
      :align: center

.. index:: asdm_decode_pop_ins.m

* **Population Decoding using Threshold-Insensitive Method** (``asdm_decode_pop_ins.m``) |lazar_perfect_2004|_ |lazar_information_2007|_:

  Reconstructs a bandlimited signal encoded by multiple Asynchronous
  Sigma-Delta Modulators using sinc kernels. This reconstruction
  method does not require the specification of an integrator
  thresholds.

   .. image:: images/tdm-sinc-ins-miso.png
      :scale: 60
      :align: center

.. include:: bibliography.rst
