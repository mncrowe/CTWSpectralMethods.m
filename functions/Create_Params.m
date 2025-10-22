function params = Create_Params(U,M2,N2,f,g,hydrostatic,free_surface,nu,Uy,Uz)
% Creates a structure containing the physical parameters for the problem of
% wave/instability modes on a shelf with a background flow. Can enter a
% single argument, n, to use a preset.
%
% Inputs:
%
% - U: function, U = U(y,z), describing along shelf (x) flow
% - N2: function, N2 = N2(y,z), describing vertical stratification N^2
% - M2: function, M2 = M2(y,z), describing horizontal stratification M^2
% - f: Coriolis parameter
% - g: gravitational acceleration
% - hydrostatic: is problem hydrostatic? 1 - yes, 0 - no
% - free_surface: does problem have a (small amplitude) free surface?
% - nu: viscosity in y direction only, use to resolve critical layers
% - Uy: function Uy = Uy(y,z), describing y derivative of U, optional
% - Yz: function Uz = Uz(y,z), describing z derivative of U, optional
%
% Note: if Uy and Uz are not specified, 'Find_Modes.m' will calculate them
% numerically from U.

arguments
    U = @(y,z) 0*y
    M2 = @(y,z) 0*y
    N2 = @(y,z) 1e-6*exp(0*z)
    f = 1e-4
    g = 10
    hydrostatic = 0
    free_surface = 0
    nu = 0
    Uy = 0
    Uz = 0
end

params.U = U;
params.M2 = M2;
params.N2 = N2;
params.f = f;
params.g = g;
params.hydrostatic = hydrostatic;
params.free_surface = free_surface;
params.nu = nu;
params.Uy = Uy;
params.Uz = Uz;

end