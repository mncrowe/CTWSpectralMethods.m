% Topographic waves along a ridge with a linear background flow gradient 

close all
clear
addpath("functions")

% set parameters:

S = -0.2;
f = 1;
type_y = [4 2 2 4]; % Laguerre at ends, Chebyshev in middle

k = 1;
n = 1;
omega_0 = 0.1*1i;   % look for stationary growing modes

plot_sigma = true;
plot_modes = true;
        
% define background flow and buoyancy gradients:

U = @(y,z) S*y;
Uy = @(y,z) S + 0*y;
Uz = @(y,z) 0*y;

M2 = @(y,z) 0*y;
N2 = @(y,z) 1 + 0*y;

% build parameter structure:

params = Create_Params(U,M2,N2,f,1,0,0,0,Uy,Uz);

% plot sigma(H0) curve:

H0 = 0.01:0.01:0.9;

if plot_sigma

    % set parameters:
    
    omega = H0*0;
    Nz = 31;
    Ly = [-5 4 1 1 4];  % 4 regions of widths 4, 1, 1 & 4, starting at y = -5
    Ny = [21 11 11 21];

    % loop through H0 values:
    
    for iH = 1:length(H0)
    
        % define H and build grid:
    
        H = @(y) 1-H0(iH)*sin(pi*(y+1)/2).*(abs(y)<1);
        Hy = @(y) -H0(iH)*pi/2*cos(pi*(y+1)/2).*(abs(y)<1);
    
        grid = Create_Grid(H,Nz,2,Ly,Ny,type_y,Hy);
    
        omega_H = Find_Modes(grid,params,k,n,omega_0,[0 0],'lm');
    
        if isempty(omega_H)
            omega(iH) = 0;
        else
            [~, I] = sort(abs(real(omega_H)),1,"ascend");
            omega_H = omega_H(I);
            if abs(real(omega_H(1))) > 1e-8
                omega(iH) = 0;  % no stationary modes
            else
                omega(iH) = omega_H(1);
            end
        end
    
        disp(['H_0 = ' num2str(H0(iH)) ', omega = ' num2str(omega(iH))])
    
    end
    
    % plot growth rate as a function of H0:
    
    figure; plot(H0, imag(omega), 'k', 'LineWidth', 0.7)
    xlabel('H_0'); ylabel('\sigma'); 
    set(gca,'FontSize',12,'linewidth',0.7,'XGrid','on','YGrid','on');

end

% calculate and plot the mode structure at each growth rate peak:

if plot_modes

    % set parameters:

    Nz = 31;
    Ly = [-6 5 1 1 5];  % 4 regions of widths 5, 1, 1 & 5, starting at y = -6
    Ny = [21 21 21 21];

    % Mode 1:
    
    iH = 17;
    H = @(y) 1-H0(iH)*sin(pi*(y+1)/2).*(abs(y)<1);
    Hy = @(y) -H0(iH)*pi/2*cos(pi*(y+1)/2).*(abs(y)<1);
    
    grid = Create_Grid(H,Nz,2,Ly,Ny,type_y,Hy);
    
    [omega1,p1] = Find_Modes(grid,params,k,1,omega_0,[0 0],'lm');
    
    n1 = (sum(Ny)-length(Ny))/2 + 1;
    phi = atan(real(p1(n1,1)) / imag(p1(n1,1)));
    p1 = exp(1i*phi) .* p1;
    
    disp(num2str(omega1))
    Plot_Mode(real(p1), grid.y, grid.z, H(grid.lambda), 'y', 'z', [-3 3])
    
    % Mode 2:
    
    iH = 47;
    H = @(y) 1-H0(iH)*sin(pi*(y+1)/2).*(abs(y)<1);
    Hy = @(y) -H0(iH)*pi/2*cos(pi*(y+1)/2).*(abs(y)<1);
    
    grid = Create_Grid(H,Nz,2,Ly,Ny,type_y,Hy);
    
    [omega2,p2] = Find_Modes(grid,params,k,1,omega_0,[0 0],'lm');
    
    n1 = (sum(Ny)-length(Ny))/2 + 1;
    phi = atan(real(p2(n1,1)) / imag(p2(n1,1)));
    p2 = exp(1i*phi) .* p2;
    
    disp(num2str(omega2))
    Plot_Mode(real(p2), grid.y, grid.z, H(grid.lambda), 'y', 'z', [-3 3])
    
    % Mode 3:
    
    iH = 76;
    H = @(y) 1-H0(iH)*sin(pi*(y+1)/2).*(abs(y)<1);
    Hy = @(y) -H0(iH)*pi/2*cos(pi*(y+1)/2).*(abs(y)<1);
    
    grid = Create_Grid(H,Nz,2,Ly,Ny,type_y,Hy);
    
    [omega3,p3] = Find_Modes(grid,params,k,1,omega_0,[0 0],'lm');
    
    n1 = (sum(Ny)-length(Ny))/2 + 1;
    phi = atan(real(p3(n1,1)) / imag(p3(n1,1)));
    p3 = exp(1i*phi) .* p3;
    
    disp(num2str(omega3))
    Plot_Mode(real(p3), grid.y, grid.z, H(grid.lambda), 'y', 'z', [-3 3])

end