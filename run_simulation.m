function results = run_simulation(options)
%RUN_SIMULATION Preliminary 3-D transient composite heat conduction model.
%   results = run_simulation() runs the default carbon-epoxy/SOR example.
%   results = run_simulation(struct('material','glass-epoxy', ...
%       'solver','gauss-seidel','showPlots',false)) overrides defaults.
%   See README.md for assumptions and numerical limitations.
%
%   Coordinates are in meters, time in seconds, and temperature in kelvin.

if nargin == 0
    options = struct();
end
validateattributes(options, {'struct'}, {'scalar'});
config = struct('material','carbon-epoxy', 'solver','adaptive-sor', ...
    'gridSize',[15 15 15], 'length',0.1, 'finalTime',7200, ...
    'timeStep',10, 'initialTemperature',300, 'patchTemperature',800, ...
    'tolerance',1e-4, 'maxIterations',1000, 'showPlots',true, ...
    'plotEvery',2);
names = fieldnames(options);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(config, name)
        error('CompositeHeat:UnknownOption','Unknown option: %s',name);
    end
    config.(name) = options.(name);
end
config.material = validatestring(config.material, {'carbon-epoxy','glass-epoxy'});
config.solver = validatestring(config.solver, {'gauss-seidel','adaptive-sor'});
validateattributes(config.gridSize, {'numeric'}, ...
    {'vector','numel',3,'integer','>=',3,'finite'});
positiveFields = {'length','finalTime','timeStep','initialTemperature', ...
    'patchTemperature','tolerance'};
for index = 1:numel(positiveFields)
    validateattributes(config.(positiveFields{index}), {'numeric'}, ...
        {'scalar','real','finite','positive'});
end
validateattributes(config.maxIterations, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(config.plotEvery, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(config.showPlots, {'logical'}, {'scalar'});
stepRatio = config.finalTime/config.timeStep;
if abs(stepRatio-round(stepRatio)) > 1e-10 || stepRatio < 1
    error('CompositeHeat:TimeGrid','finalTime must be a positive integer multiple of timeStep.');
end

nx = config.gridSize(1); ny = config.gridSize(2); nz = config.gridSize(3);
spacing = config.length./(config.gridSize-1);
steps = round(stepRatio);
T = ones(nx,ny,nz)*config.initialTemperature;
patchX = round(nx/3):round(2*nx/3);
patchY = round(ny/3):round(2*ny/3);
T(patchX,patchY,1) = config.patchTemperature;
iterations = zeros(steps,1);
maxUpdate = zeros(steps,1);
converged = false(steps,1);
finalOmega = ones(steps,1);
solveSeconds = zeros(steps,1);
x = linspace(0,config.length,nx);
y = linspace(0,config.length,ny);
z = linspace(0,config.length,nz);
if config.showPlots
    [X,Y,Z] = meshgrid(x*1000,y*1000,z*1000);
    figure3D = figure('Name','Transient composite heat conduction');
end

simulationTimer = tic;
for step = 1:steps
    [kx,ky,kz,cp,rho] = material_properties(T,config.material);
    if any(~isfinite(kx(:)) | kx(:)<=0) || any(~isfinite(ky(:)) | ky(:)<=0) || ...
            any(~isfinite(kz(:)) | kz(:)<=0) || any(~isfinite(cp(:)) | cp(:)<=0) || ...
            any(~isfinite(rho(:)) | rho(:)<=0)
        error('CompositeHeat:Properties','Property formulas produced nonpositive or nonfinite values.');
    end
    ax = kx/spacing(1)^2;
    ay = ky/spacing(2)^2;
    az = kz/spacing(3)^2;
    inertia = rho.*cp/config.timeStep;
    diagonal = 2*(ax+ay+az)+inertia;
    rhs = inertia.*T;
    solverTimer = tic;
    [T,iterations(step),maxUpdate(step),finalOmega(step)] = ...
        solve_step(T,ax,ay,az,diagonal,rhs,config);
    solveSeconds(step) = toc(solverTimer);
    converged(step) = isfinite(maxUpdate(step)) && maxUpdate(step)<=config.tolerance;
    if any(~isfinite(T(:)))
        error('CompositeHeat:NonfiniteTemperature','Nonfinite temperature at time step %d.',step);
    end

    % Original boundary treatment is retained: copy AFTER the interior solve.
    % This is a lagged approximation, not a coupled zero-flux boundary solve.
    T(1,:,:) = T(2,:,:); T(nx,:,:) = T(nx-1,:,:);
    T(:,1,:) = T(:,2,:); T(:,ny,:) = T(:,ny-1,:);
    T(:,:,1) = T(:,:,2); T(:,:,nz) = T(:,:,nz-1);
    T(patchX,patchY,1) = config.patchTemperature;

    if config.showPlots && isgraphics(figure3D) && ...
            (mod(step,config.plotEvery)==0 || step==steps)
        figure(figure3D); clf(figure3D);
        h = slice(X,Y,Z,permute(T,[2 1 3]), ...
            config.length*500,config.length*500,config.length*500);
        set(h,'EdgeColor','none');
        colormap hot;
        cb = colorbar; cb.Label.String = 'Temperature (K)';
        clim(temperature_limits(config));
        xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
        axis equal tight; view(-45,30);
        title(sprintf('Time: %.0f s | %s',step*config.timeStep,config.material));
        drawnow;
    end
end
elapsed = toc(simulationTimer);
if any(~converged)
    warning('CompositeHeat:IterationLimit', ...
        '%d time steps did not meet the update tolerance. Inspect results.history.',sum(~converged));
end
if config.showPlots
    figure('Name','Temperature at different depths');
    layerIndices = [1 round(nz/2) nz];
    for index = 1:3
        layer = layerIndices(index);
        subplot(1,3,index);
        contourf(x*1000,y*1000,T(:,:,layer)',20,'LineStyle','none');
        colormap hot; clim(temperature_limits(config));
        cb = colorbar; cb.Label.String = 'Temperature (K)';
        xlabel('X (mm)'); ylabel('Y (mm)'); axis equal tight;
        title(sprintf('Z = %.1f mm',z(layer)*1000));
    end
end
results = struct('config',config,'temperature',T,'x',x,'y',y,'z',z, ...
    'elapsedSeconds',elapsed,'solverSeconds',sum(solveSeconds));
results.history = table((1:steps)'*config.timeStep,iterations,maxUpdate, ...
    converged,finalOmega,solveSeconds, ...
    'VariableNames',{'TimeSeconds','Iterations','MaxUpdateKelvin', ...
    'Converged','FinalOmega','SolveSeconds'});
end

function limits = temperature_limits(config)
limits = sort([config.initialTemperature config.patchTemperature]);
if limits(1)==limits(2)
    limits = limits+[-0.5 0.5];
end
end

function [kx,ky,kz,cp,rho] = material_properties(T,material)
% Illustrative formulas retained from the original project, not sourced data.
deltaT = T-300;
if strcmp(material,'carbon-epoxy')
    k1 = 4.5+0.005*deltaT;
    k2 = 0.6+0.001*deltaT-1e-6*deltaT.^2;
    cp = 900+2.5*deltaT;
    rho = 1550-0.05*deltaT;
else
    k1 = 1.2+0.002*deltaT;
    k2 = 0.3+0.0005*deltaT;
    cp = 1100+3.0*deltaT;
    rho = 1900-0.08*deltaT;
end
angles = [0 45 90 -45];
nz = size(T,3);
nodesPerLayer = ceil(nz/numel(angles));
kx = zeros(size(T)); ky = kx; kz = k2;
for layer = 1:nz
    angle = angles(min(floor((layer-1)/nodesPerLayer)+1,numel(angles)))*pi/180;
    kx(:,:,layer) = k1(:,:,layer)*cos(angle)^2+k2(:,:,layer)*sin(angle)^2;
    ky(:,:,layer) = k1(:,:,layer)*sin(angle)^2+k2(:,:,layer)*cos(angle)^2;
end
end

function [T,iteration,maxUpdate,usedOmega] = solve_step(T,ax,ay,az,diagonal,rhs,config)
[nx,ny,nz] = size(T);
omega = 1; increment = 0.05; previousUpdate = inf; patience = 0;
adaptive = strcmp(config.solver,'adaptive-sor');
for iteration = 1:config.maxIterations
    maxUpdate = 0;
    usedOmega = omega;
    for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                oldValue = T(i,j,k);
                neighbors = ax(i,j,k)*(T(i-1,j,k)+T(i+1,j,k)) + ...
                    ay(i,j,k)*(T(i,j-1,k)+T(i,j+1,k)) + ...
                    az(i,j,k)*(T(i,j,k-1)+T(i,j,k+1));
                gsValue = (neighbors+rhs(i,j,k))/diagonal(i,j,k);
                T(i,j,k) = (1-omega)*oldValue+omega*gsValue;
                maxUpdate = max(maxUpdate,abs(T(i,j,k)-oldValue));
            end
        end
    end
    if maxUpdate<=config.tolerance
        break;
    end
    if adaptive && iteration>6
        ratio = maxUpdate/previousUpdate;
        if ratio<0.98
            patience = patience+1;
            if patience>=3
                omega = min(omega+increment,1.9);
                patience = 0;
            end
        elseif ratio>=1
            omega = max(1,omega-2*increment);
            increment = increment/2;
            patience = 0;
        end
    end
    previousUpdate = maxUpdate;
end
end
