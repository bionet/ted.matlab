% Time Encoding and Decoding Toolkit
%
% Asynchronous Sigma-Delta Modulator Algorithms
%   asdm_decode                  - ASDM time decoder.
%   asdm_decode_fast             - Fast ASDM time decoder.
%   asdm_decode_ins              - Parameter-insensitive ASDM time decoder.
%   asdm_decode_pop              - Population-based ASDM time decoder.
%   asdm_decode_pop_ins          - Population-based parameter insensitive 
%                                  ASDM time decoder.
%   asdm_encode                  - ASDM time encoder.
%
% Integrate-and-Fire Neuron Algorithms
%   iaf_decode                   - IAF time decoder that uses sinc kernels.
%   iaf_decode_fast              - Fast IAF time decoder.
%   iaf_decode_filt_trig_pop     - Population-based IAF time decoder that uses the 
%                                  trigonometric polynomial approximation.
%   iaf_decode_pop               - Population-based IAF time decoder that uses 
%                                  sinc kernels.
%   iaf_encode                   - IAF time encoder.
%
% Delayed Integrate-and-Fire Neuron Algorithms
%   G_block_delay                - Compute reconstruction matrix.
%   iaf_decode_pop_delay         - Population-based IAF time decoder with delays. 
%
% Spline Interpolation Functions
%   consistent_decoding_IF_MIMO  - Multiple input multiple IAF neuron time decoder.
%   consistent_decoding_IF_ONOFF - ON-OFF IAF neuron pair time decoder.
%   consistent_decoding_LIF      - IAF time decoder.
%   iaf_encode_ideal_on_off      - ON-OFF IAF neuron pair time encoder.
%
% Smoothing Splines Reconstruction Functions
%   LIF_decode_S1                - LIF time decoder using smoothing splines in
%                                  S1 space.
%   LIF_decode_S2                - LIF time decoder using smoothing splines in 
%                                  S2 space.
% See also TEDDEMOS

% Author: Lev Givon
% Copyright 2009-2015 Lev Givon
