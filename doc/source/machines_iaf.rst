.. -*- rst -*- 

Integrate-and-Fire Neuron Machines
==================================

Encoding with Integrate-and-Fire (IAF) Neurons
----------------------------------------------
.. index:: iaf_encode.m

* **IAF Encoding** (``iaf_encode.m``) |lazar_perfect_2004|_ |lazar_reconstruction_2009|_:

  Encodes a time-varying signal using an IAF neuron. Leaky and ideal
  neuron models are supported. In addition, IAF neuron with random
  (Gaussian) thresholds is also supported. A signal can be encoded by
  a single IAF encoder (Single-Input Single-Output Encoding), as shown
  in the figure below, or it can be encoded by a population of IAF
  neurons (Single-Input Single-Output Encoding).

   .. image:: images/tem-iaf-rt.png
      :scale: 60
      :align: center

.. index:: iaf_encode_on_off.m

* **ON-OFF IAF Encoding** (``iaf_encode_ideal_on_off.m``) |lazar_consistent_2009|_:

  Encodes a time-varying signal using an ON-OFF IAF neuron pair. Only
  ideal IAF neuron models are supported.

   .. image:: images/tem-iaf-coupled.png
      :scale: 60
      :align: center

Decoding for Signal Encoded with Single-Input Single-Output IAF neuron
----------------------------------------------------------------------
.. index:: iaf_decode.m

* **Decoding** (``iaf_decode.m``) |lazar_time_2004|_:

  Reconstructs a bandlimited signal encoded by an IAF neuron using
  sinc kernels.

   .. image:: images/tdm-sinc.png
      :scale: 60
      :align: center

.. index:: iaf_decode_fast.m

* **Decoding using Fast Approximation Method** (``iaf_decode_fast.m``) |lazar_fast_2005|_:

  Reconstructs a bandlimited signal encoded by an IAF neuron using a
  fast approximation method.

   .. image:: images/tdm-fast.png
      :scale: 60
      :align: center

.. index:: consistent_decoding_LIF.m

* **Decoding using Spline Interpolation** (``consistent_decoding_LIF.m``) |lazar_consistent_2009|_:

  Reconstructs a finite-energy signal encoded by an
  Leaky-Integrate-and-Fire (LIF) neuron. It uses spline interpolation
  algorithm.

   .. image:: images/tdm-spline.png
      :scale: 60
      :align: center

.. index:: LIF_decode_S1.m, LIF_decode_S2.m

* **Decoding using Smoothing Spline** (``LIF_decode_S1.m``, ``LIF_decode_S2.m``) |lazar_reconstruction_2009|_:

  Reconstructs a signal in Sobolev space :math:`S_1` or :math:`S_2`
  encoded by a LIF neuron using smoothing splines. Signals encoded by
  a LIF with random threshold should be decoded using this function.

   .. image:: images/tdm-spline-smoothing.png
      :scale: 60
      :align: center

Decoding for Signal Encoded with Single-Input Multiple-Output IAF neurons
-------------------------------------------------------------------------
.. index:: consistent_decoding_IF_ONOFF.m

* **Decoding for ON-OFF IAF** (``consistent_decoding_IF_ONOFF.m``) |lazar_consistent_2009|_:

  Reconstructs a finite energy signal encoded by ON-OFF IAF neuron
  pair. The reconstruction is performed using spline interpolation
  method.

   .. image:: images/tdm-spline-mimo.png
      :scale: 60
      :align: center

.. index:: iaf_decode_pop.m

* **Population Decoding** (``iaf_decode_pop.m``) |lazar_information_2007|_:

  Reconstructs a bandlimited signal encoded by an ensemble of IAF neurons using sinc kernels.

   .. image:: images/tdm-sinc-miso.png
      :scale: 60
      :align: center

.. index:: LIF_pop_decode_S1.m, LIF_pop_decode_S2.m

* **Population Decoding using Smoothing Splines**
  (``LIF_pop_decode_S1.m``, ``LIF_pop_decode_S2.m``) |lazar_reconstruction_2009|_:

  Reconstructs a signal in Sobolev Space :math:`S_1` or :math:`S_2`
  encoded by a population of LIF neurons. The reconstruction uses
  smoothing spline method in the RKHS. Signals encoded by population
  of LIF with random threshold should be decoded using this function.

   .. image:: images/tdm-spline-smoothing-mimo.png
      :scale: 60
      :align: center

Decoding for Signals Encoded with Multiple-Input Multiple-Output IAF neurons
----------------------------------------------------------------------------
.. index:: consistent_decoding_IF_MIMO.m

* **Population Decoding using Spline Interpolation** (``consistent_decoding_IF_MIMO.m``) |lazar_consistent_2009|_:

  Reconstructs multiple finite energy signals encoded by a population
  of ideal IAF neurons in Multiple-Input Multiple Output setting. The
  reconstruction uses spline interpolation method.

   .. image:: images/tdm-spline-mimo.png
      :scale: 60
      :align: center

.. include:: bibliography.rst
