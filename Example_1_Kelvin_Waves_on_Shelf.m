% Considers the case of Kelvin waves on a sloping boundary

close all
clear
addpath("functions")

% define grid and parameters:
Ny = 31;
Nz = 21;
Ly = 4;
H0 = 0.9;
H = @(y) H0 + (1-H0)*tanh(y);

grid = Create_Grid(H,Nz,2,Ly,Ny,4);
params = Create_Params(@(y,z) 0*y,@(z) 1+0*z,1,1,0,0);

% set number of modes and define solver:
n = 6;
Omega_Eval = @(k,n,w0) Find_Modes(grid,params,k,n,w0,0,'lm');

% create wavenumber and frequency arrays:
k = 0.1:0.1:3;
omega = zeros(n,length(k));
w0 = 0.50; %H0*k(1)/pi;

% find initial points on frequency curves:
disp(['Finding modes for k = ' num2str(k(1))])
omega(:,1) = Omega_Eval(0.1,n,w0);

% follow curves through space:
for ik = 2:length(k)
    disp(['Finding modes for k = ' num2str(k(ik))])
    omega(:,ik) = Omega_Eval(k(ik),n,w0);
end

% calculate phase speed and interpolate back to k = 0:
c_p = omega./k;
c_p = [interp1(k,c_p',0,'spline','extrap')' c_p];

% plot c_p curves for each mode:
plot([0 k],c_p,'-','LineWidth',1)
grid on
ylim([0 0.4])
ylabel('c_p'); xlabel('k');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','NumColumns',2)
set(gca,'FontSize',12,'linewidth',0.7);

% surface plots of the pressure signals as a function of y and z:
plot_mode = 1;  % mode number to show
if plot_mode ~= 0
    [omega_01,p_01] = Find_Modes(grid,params,0.1,n,w0,0);
    [omega_1,p_1] = Find_Modes(grid,params,1,n,w0,0);
    [omega_3,p_3] = Find_Modes(grid,params,3,n,w0,0);
    Plot_Mode(p_01(:,:,plot_mode),grid.y,grid.z);
    Plot_Mode(p_1(:,:,plot_mode),grid.y,grid.z);
    Plot_Mode(p_3(:,:,plot_mode),grid.y,grid.z);
end
