% =========================================================================
% Project: 3D Transient Heat Conduction in Anisotropic Composites
% Description: ULTIMATE VERSION (Laser Patch Heating & Dynamic Units)
% Features: 
%   1. Temperature-dependent properties (Non-linear)
%   2. Custom Adaptive SOR Solver (Self-tuning acceleration)
%   3. Complex Boundary Conditions (Adiabatic + Localized Laser Patch)
%   4. 3D and 2D Visualization with Dynamic Millimeter (mm) Scaling
% =========================================================================
clear; clc; close all;

%% 1. User Interface (Numeric Choices)
disp('--- 3D Composite Heat Transfer Simulator (Ultimate Version) ---');
disp('Select Material Type:');
disp('1: Carbon-Epoxy (High Anisotropy - Fast X/Y conduction)');
disp('2: Glass-Epoxy (Moderate Anisotropy - Slower conduction)');
mat_choice = input('Enter your choice (1 or 2): ');

disp('Select Numerical Solver:');
disp('1: Standard Gauss-Seidel (Safe, Constant Speed)');
disp('2: Adaptive SOR (AI-like Self-Tuning, High Speed)');
solver_choice = input('Enter your choice (1 or 2): ');

%% 2. Physical and Grid Settings
nx = 15; ny = 15; nz = 15; % Grid dimensions
L = 0.1; % Length of the cubic domain [m]
dx = L/(nx-1); dy = L/(ny-1); dz = L/(nz-1);

t_final = 7200; % Total simulation time: 2 Hours (7200 seconds)
dt = 10;        % Time step increased to 10 seconds for fast computation
num_steps = round(t_final/dt);

% Initialization (Uniform temperature of 300 K)
T_old = ones(nx, ny, nz) * 300; 

% Define Localized Patch Heating (Center of the bottom surface Z=1)
patch_x = round(nx/3):round(2*nx/3);
patch_y = round(ny/3):round(2*ny/3);
T_old(patch_x, patch_y, 1) = 800; % Initial thermal shock at the patch
T_new = T_old;

% Solver settings
tol = 1e-4;
max_iter = 1000;

%% 3. Graphics Preparation
[X, Y, Z] = meshgrid(linspace(0,L,nx), linspace(0,L,ny), linspace(0,L,nz));
fig3D = figure('Name', '3D Transient Heat Transfer', 'Position', [100, 100, 800, 600]);

disp('Starting advanced simulation... Please wait.');

%% 4. Main Transient Loop (Time Marching)
for n = 1:num_steps
    
    % Step A: Update non-linear material properties
    [kx, ky, kz, cp, rho] = local_update_properties(T_old, mat_choice);
    
    % Step B: Assemble FDM coefficients
    a_W = kx / dx^2;  a_E = kx / dx^2;
    a_S = ky / dy^2;  a_N = ky / dy^2;
    a_B = kz / dz^2;  a_T = kz / dz^2;
    
    thermal_inertia = (rho .* cp) / dt;
    a_P = a_W + a_E + a_S + a_N + a_B + a_T + thermal_inertia;
    b = thermal_inertia .* T_old;
    
    % Step C: Execute chosen Solver
    if solver_choice == 1
        T_new = local_gauss_seidel(T_old, a_W, a_E, a_S, a_N, a_B, a_T, a_P, b, max_iter, tol);
    else
        T_new = local_adaptive_sor(T_old, a_W, a_E, a_S, a_N, a_B, a_T, a_P, b, max_iter, tol);
    end
    
    % Step D: Apply Complex Boundary Conditions
    % 1. Adiabatic (Insulated) Boundaries for all 6 faces
    T_new(1, :, :) = T_new(2, :, :);     % Left
    T_new(nx, :, :) = T_new(nx-1, :, :); % Right
    T_new(:, 1, :) = T_new(:, 2, :);     % Rear
    T_new(:, ny, :) = T_new(:, ny-1, :); % Front
    T_new(:, :, 1) = T_new(:, :, 2);     % Bottom
    T_new(:, :, nz) = T_new(:, :, nz-1); % Top
    
    % 2. Constant Localized Laser Heating
    T_new(patch_x, patch_y, 1) = 800; 
    
    % Advance in time
    T_old = T_new;
    
    % Step E: Dynamic 3D Visualization (Update every 2 steps)
    if mod(n, 2) == 0
        figure(fig3D); clf;
        
        % FIX: Permute T_new to align its X/Y dimensions with meshgrid
        T_plot_3D = permute(T_new, [2, 1, 3]);
        
        % Scale coordinates to millimeters (*1000)
        h = slice(X*1000, Y*1000, Z*1000, T_plot_3D, (L/2)*1000, (L/2)*1000, (L/2)*1000);
        
        set(h, 'EdgeColor', 'none');
        shading interp; colormap hot; colorbar; clim([300, 800]);
        view(-45, 30); camlight; lighting gouraud;
        
        xlabel('X-Axis'); ylabel('Y-Axis'); zlabel('Z-Axis');
        xtickformat('%g mm'); ytickformat('%g mm'); ztickformat('%g mm');
        
        title(sprintf('Simulation Time: %.1f seconds', n*dt), 'FontWeight', 'bold');
        drawnow;
    end
end

disp('Simulation Finished Successfully!');

%% 5. Layer-by-Layer 2D Visualization (Post-Processing)
disp('Generating layer-by-layer 2D plots with dynamic units...');
fig2D = figure('Name', 'Layer-by-Layer Distribution', 'Position', [150, 150, 1200, 400]);

layer_indices = [1, round(nz/2), nz];
layer_titles = {'Surface Layer (Z = 0)', 'Middle Layer (Z = L/2)', 'Back Layer (Z = L)'};

for i = 1:3
    k = layer_indices(i); 
    
    % FIX: Squeeze the 3D matrix into 2D, then transpose it to align with X/Y grids
    T_2D = squeeze(T_new(:, :, k))'; 
    
    subplot(1, 3, i);
    contourf(X(:,:,1)*1000, Y(:,:,1)*1000, T_2D, 20, 'LineStyle', 'none'); 
    colormap hot; colorbar;
    clim('auto');
    
    xlabel('Width'); ylabel('Height'); 
    xtickformat('%g mm'); ytickformat('%g mm');
    
    title(layer_titles{i}, 'FontWeight', 'bold');
    axis square;
end

% =========================================================================
% ======================= LOCAL FUNCTIONS =================================
% =========================================================================
function [kx, ky, kz, cp, rho] = local_update_properties(T, mat_choice)
    [nx, ny, nz] = size(T);
    T_ref = 300; 
    
    kx = zeros(nx, ny, nz);
    ky = zeros(nx, ny, nz);
    kz = zeros(nx, ny, nz);
    cp = zeros(nx, ny, nz);
    rho = zeros(nx, ny, nz);
    
    % Laminate Stacking Sequence [0 / 45 / 90 / -45] degrees
    stacking_angles = [0, 45, 90, -45]; 
    num_layers = length(stacking_angles);
    nodes_per_layer = ceil(nz / num_layers);
    
    if mat_choice == 1 % Carbon-Epoxy
        k1_base = 4.5 + 0.005 .* (T - T_ref); 
        k2_base = 0.6 + 0.001 .* (T - T_ref) - 1e-6 .* (T - T_ref).^2;
        cp_base = 900 + 2.5 .* (T - T_ref);
        rho_base = 1550 - 0.05 .* (T - T_ref);
    elseif mat_choice == 2 % Glass-Epoxy
        k1_base = 1.2 + 0.002 .* (T - T_ref);
        k2_base = 0.3 + 0.0005 .* (T - T_ref);
        cp_base = 1100 + 3.0 .* (T - T_ref);
        rho_base = 1900 - 0.08 .* (T - T_ref);
    end
    
    for m = 1:nz
        layer_idx = min(floor((m-1)/nodes_per_layer) + 1, num_layers);
        theta = stacking_angles(layer_idx) * (pi / 180); 
        
        m_cos = cos(theta);
        n_sin = sin(theta);
        
        kx(:, :, m) = k1_base(:, :, m) * m_cos^2 + k2_base(:, :, m) * n_sin^2;
        ky(:, :, m) = k1_base(:, :, m) * n_sin^2 + k2_base(:, :, m) * m_cos^2;
        kz(:, :, m) = k2_base(:, :, m); 
        
        cp(:, :, m) = cp_base(:, :, m);
        rho(:, :, m) = rho_base(:, :, m);
    end
end

function T = local_gauss_seidel(T, a_W, a_E, a_S, a_N, a_B, a_T, a_P, b, max_iter, tol)
    [nx, ny, nz] = size(T);
    iter = 0; max_err = inf;
    while (max_err > tol) && (iter < max_iter)
        max_err = 0;
        for m = 2:nz-1
            for j = 2:ny-1
                for i = 2:nx-1
                    T_old = T(i,j,m);
                    sum_neighbors = a_W(i,j,m)*T(i-1,j,m) + a_E(i,j,m)*T(i+1,j,m) + ...
                                    a_S(i,j,m)*T(i,j-1,m) + a_N(i,j,m)*T(i,j+1,m) + ...
                                    a_B(i,j,m)*T(i,j,m-1) + a_T(i,j,m)*T(i,j,m+1);
                    T_new = (sum_neighbors + b(i,j,m)) / a_P(i,j,m);
                    T(i,j,m) = T_new; 
                    err = abs(T_new - T_old);
                    if err > max_err; max_err = err; end
                end
            end
        end
        iter = iter + 1;
    end
end

function T = local_adaptive_sor(T, a_W, a_E, a_S, a_N, a_B, a_T, a_P, b, max_iter, tol)
    [nx, ny, nz] = size(T);
    iter = 0; max_err = inf;
    
    omega = 1.0; 
    delta_omega = 0.05;
    prev_err = inf; patience = 0;
    
    while (max_err > tol) && (iter < max_iter)
        max_err = 0;
        for m = 2:nz-1
            for j = 2:ny-1
                for i = 2:nx-1
                    T_old = T(i,j,m);
                    sum_neighbors = a_W(i,j,m)*T(i-1,j,m) + a_E(i,j,m)*T(i+1,j,m) + ...
                                    a_S(i,j,m)*T(i,j-1,m) + a_N(i,j,m)*T(i,j+1,m) + ...
                                    a_B(i,j,m)*T(i,j,m-1) + a_T(i,j,m)*T(i,j,m+1);
                    
                    T_GS = (sum_neighbors + b(i,j,m)) / a_P(i,j,m);
                    T_new = (1 - omega) * T_old + omega * T_GS; 
                    
                    T(i,j,m) = T_new; 
                    err = abs(T_new - T_old);
                    if err > max_err; max_err = err; end
                end
            end
        end
        
        if iter > 5
            error_ratio = max_err / prev_err;
            if error_ratio < 0.98 
                patience = patience + 1;
                if patience >= 3
                    omega = min(omega + delta_omega, 1.90); 
                    patience = 0;
                end
            elseif error_ratio >= 1.0 
                omega = max(1.0, omega - 2*delta_omega); 
                delta_omega = delta_omega / 2;
                patience = 0;
            end
        end
        prev_err = max_err;
        iter = iter + 1;
    end
end