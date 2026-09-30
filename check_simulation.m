function check_simulation()
%CHECK_SIMULATION Small behavioral checks, not physical validation.
common = struct('gridSize',[7 7 7], 'finalTime',30, ...
    'showPlots',false,'tolerance',1e-8);
uniform = common;
uniform.patchTemperature = 300;
r = run_simulation(uniform);
assert(max(abs(r.temperature(:)-300))<1e-9, 'Uniform state changed.');
assert(all(r.history.Converged), 'Uniform solve did not converge.');
for material = {'carbon-epoxy','glass-epoxy'}
    common.material = material{1};
    common.solver = 'gauss-seidel';
    gs = run_simulation(common);
    common.solver = 'adaptive-sor';
    sor = run_simulation(common);
    assert(all(gs.history.Converged) && all(sor.history.Converged));
    assert(max(abs(gs.temperature(:)-sor.temperature(:)))<1e-6, ...
        'Solvers disagree for the same discrete problem.');
    assert(all(isfinite(sor.temperature(:))));
    assert(min(sor.temperature(:))>=300-1e-8 && max(sor.temperature(:))<=800+1e-8);
end
bad = common;
bad.finalTime = 31;
try
    run_simulation(bad);
    error('CompositeHeat:MissingError','Invalid time grid was accepted.');
catch exception
    assert(strcmp(exception.identifier,'CompositeHeat:TimeGrid'));
end
fprintf('Behavioral checks passed. These checks do not validate the physical model.\n');
end
