function params = Create_Params(U,N2,f,g,hydrostatic,free_surface)
% Creates a structure containing the physical parameters for the problem of
% wave/instability modes on a shelf with a background flow. Can enter a
% single argument, n, to use a preset.
%
% Inputs:
%
% - n: preset;
%       0 - Kelvin waves with no flow
%       1 - ...
%       2 - ...
%
% OR
%
% - U: function, U = U(y,z), describing along shelf (x) flow
% - N2: function, N2 = N2(z), describing background stratification N^2
% - f: Coriolis parameter
% - g: gravitational acceleration
% - hydrostatic: is problem hydrostatic? 1 - yes, 0 - no
% - free_surface: does problem have a (small amplitude) free surface?

if nargin == 1
    if U == 0   % preset 0
        params.U = @(y,z) 0*y;
        params.N2 = @(z) 1e-6*exp(0*z);
        params.f = 1e-4;
        params.g = 10;
        params.hydrostatic = 0;
        params.free_surface = 0;
    end
    if U == 1   % preset 1
        
    end
    if U == 2   % preset 2

    end
end

if nargin == 6
    params.U = U;
    params.N2 = N2;
    params.f = f;
    params.g = g;
    params.hydrostatic = hydrostatic;
    params.free_surface = free_surface;
end

if nargin > 1 && nargin < 6
    error('Insufficient parameters entered.')
end

end