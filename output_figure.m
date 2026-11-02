function [fig,cleanup]=output_figure(varargin)
% Create an export figure owned by the caller's function scope.
% Retain cleanup until exporting finishes; user figures are never closed.
fig=figure(varargin{:});
cleanup=onCleanup(@() finish_figure(fig));
end

function finish_figure(fig)
if isgraphics(fig,'figure'), delete(fig); end
end
