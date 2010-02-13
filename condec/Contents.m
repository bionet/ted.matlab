% Consistent Decoding Algorithms for Finite Energy Stimuli Encoded with LIF
% Neurons Using Interpolation Splines. Functions and notation follow [1]
%
%   consistent_decoding_LIF       - Reconstruction of signal encoded with a LIF single neuron [Section 2]
%   consistent_decoding_IF_ONOFF  - Reconstruction of signal encoded with an ON-OFF neuron pair [Section 3]
%   consistent_decoding_IF_MIMO   - Reconstruction of a vector valued signal encoded with a population of neurons [Section 4]
%   G_IF                          - Construction of G matrix for the single neuron case (both leaky and ideal) eqs. (17) and (19)
%   G_pop_IF                      - Construction of G matrix for population of neurons
%   Gblock_IF                     - Calculation of G^{ij} for ideal IF neurons eq. (51)
%   p_IF                          - General calculation of p-vector
%   r_IF                          - General calculation of r-vector
%   psi                           - General calculation of reconstruction function psi eqs. (16) and (19) 
%   iaf_encode_ideal_on_off       - Encoding with a symmetric ON-OFF idel IF neuron pair
%   cross_fb                      - Specifying the cross-feedback term
%   q_IF_ONOFF                    - t-transform for the ON-OFF neuron pair

%
% See also condecdemos
%
% References
%
% [1] Aurel A. Lazar and Eftychios A. Pnevmatikakis, Consistent Recovery of
% Sensory Stimuli Encoded with MIMO Neural Circuits, Computational Intelligence 
% and Neuroscience, Volume 2010, February, 2010, Special Issue on Signal 
% Processing for Neural Spike Trains, Article ID 469658, 
% Published Online September 22, 2009.

% Author: Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University