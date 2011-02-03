% Time Encoding and Decoding Toolbox
% Version 0.04
%
% Basic Asynchronous Sigma-Delta Modulator Algorithms
%   asdm_decode                  - ASDM time decoder.
%   asdm_decode_fast             - Fast ASDM time decoder.
%   asdm_decode_ins              - Parameter-insensitive ASDM time decoder.
%   asdm_decode_pop              - Population-based ASDM time decoder.
%   asdm_decode_pop_ins          - Population-based parameter insensitive ASDM time decoder.
%   asdm_encode                  - ASDM time encoder.
%
% Basic Integrate-and-Fire Neuron Algorithms
%   iaf_decode                   - IAF time decoder.
%   iaf_decode_fast              - Fast IAF time decoder.
%   iaf_decode_pop               - Population-based IAF time decoder.
%   iaf_encode                   - IAF time encoder.
%
% Delayed Integrate-and-Fire Neuron Algorithms
%   G_block_delay                - Compute reconstruction matrix.
%   iaf_decode_pop_delay         - Population-based IAF time decoder with delays. 
%
% Gammatone Filter Functions
%   gammatone                    - Create a gammatone filter bank.
%
% Spline Interpolation Functions
%   consistent_decoding_IF_MIMO  - Multiple input multiple IAF neuron time decoder.
%   consistent_decoding_IF_ONOFF - ON-OFF IAF neuron pair time decoder.
%   consistent_decoding_LIF      - IAF time decoder.
%   G_IF                         - Compute reconstruction matrix for single IAF neuron decoder.
%   G_block_IF                   - Compute part of reconstruction matrix for multiple IAF neuron decoder.
%   G_pop_IF                     - Compute reconstruction matrix for multiple IAF neuron decoder.
%   cross_fb                     - Compute cross-feedback for ON-OFF IAF neuron pair.
%   iaf_encode_ideal_on_off      - ON-OFF IAF neuron pair time encoder.
%   p_IF                         - Compute inner product needed by spline interpolation decoders.
%   psi                          - Compute signal reconstruction function needed by spline interpolation decoders.
%   q_IF_ONOFF                   - Compute the t-transform for an ON-OFF IAF neuron pair
%   r_IF                         - Compute inner product needed by spline interpolation decoders.
%
% See also TEDDEMOS

% Author: Lev Givon
% Copyright 2009-2011 Lev Givon
