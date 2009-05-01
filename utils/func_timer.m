%FUNC_TIMER Time a function invocation.
%   FUNC_TIMER(F,ARG1,ARG2,...) times the execution of function
%   with handle F called with arguments ARG1, ARG2, etc.

function [varargout] = func_timer(f,varargin)

varargout = cell(1,nargout);
t = tic();
[varargout{:}] = f(varargin{:});
t = toc(t);  
fprintf(1,'execution time = %.3f s\n',t);

