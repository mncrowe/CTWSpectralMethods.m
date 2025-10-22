function params = Create_Params_Barotropic(U,f,N2,hydrostatic,nu)
% Creates parameter structure for the Barotropic and Equivalent barotropic
% problems. The parameters 'N2' and 'hydrostatic' are only used for the
% equivalent barotropic problem.

arguments
    U = @(y) erf(2*(y-2))
    f = 1;
    N2 = 1
    hydrostatic = 0
    nu = 0
end

params.U = U;
params.f = f;
params.N2 = N2;
params.hydrostatic = hydrostatic;
params.nu = nu;

end