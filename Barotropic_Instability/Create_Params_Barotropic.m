function params = Create_Params_Barotropic(U,f,N2,hydrostatic)
% Creates parameter structure for the Barotropic and Equivalent barotropic
% problems. The parameters 'N2' and 'hydrostatic' are only used for the
% equivalent barotropic problem.

if nargin < 2; f = 1; end
if nargin < 1; U = @(y) erf(2*(y-2)); end
if nargin < 3; N2 = 1; end
if nargin < 4; hydrostatic = 0; end

params.U = U;
params.f = f;
params.N2 = N2;
params.hydrostatic = hydrostatic;

end