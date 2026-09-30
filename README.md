# Transient Heat Conduction in Layered Composites

A MATLAB learning project connecting concepts from **Mechanics of Composite Materials** and **Numerical Methods**. It explores how localized surface heating spreads through a three-dimensional layered composite with direction-dependent, temperature-dependent properties.

This is a preliminary numerical model. Its outputs have not been validated against an analytical solution, a reference solver, or experimental measurements.

## What the code does

- Offers illustrative carbon-epoxy and glass-epoxy property formulas.
- Represents a nominal `[0 / 45 / 90 / -45]` stacking sequence using directional conductivity components.
- Uses a seven-point finite difference stencil and backward-Euler-style time stepping with properties evaluated at the previous time level.
- Provides Gauss-Seidel and adaptive successive over-relaxation (SOR) solvers.
- Displays temperature slices during the simulation and temperature maps at three depths.
- Returns the final temperature field and per-step iteration, update, timing, and convergence information.

The SOR adjustment is a rule-based numerical heuristic, not an AI or machine-learning method. No performance advantage is claimed without a controlled comparison.

## Requirements and quick start

Use MATLAB R2022a or newer (`clim` is used for plotting). No additional MATLAB toolboxes are required by the code. Execution checks for this release are described in `VALIDATION.md`.

Download or clone this repository, open its folder in MATLAB, and run:

```matlab
results = run_simulation();
```

The default run uses carbon-epoxy and adaptive SOR. For a short, non-graphical run:

```matlab
options = struct('finalTime', 120, 'showPlots', false);
results = run_simulation(options);
disp(results.history);
```

To use the other material and solver:

```matlab
options = struct('material', 'glass-epoxy', 'solver', 'gauss-seidel');
results = run_simulation(options);
```

Save results locally:

```matlab
if ~isfolder('results'), mkdir('results'); end
save(fullfile('results', 'simulation.mat'), 'results');
writetable(results.history, fullfile('results', 'convergence.csv'));
```

The `results/` directory is ignored by Git. Figures are not automatically saved.

## Default configuration

| Option | Default | Meaning |
| --- | --- | --- |
| `material` | `'carbon-epoxy'` | Alternative: `'glass-epoxy'` |
| `solver` | `'adaptive-sor'` | Alternative: `'gauss-seidel'` |
| `gridSize` | `[15 15 15]` | Nodes along x, y, z |
| `length` | `0.1` | Cube side length, m |
| `finalTime` | `7200` | Simulated duration, s |
| `timeStep` | `10` | Time step, s |
| `initialTemperature` | `300` | Initial temperature outside the patch, K |
| `patchTemperature` | `800` | Prescribed patch temperature, K |
| `tolerance` | `1e-4` | Maximum interior temperature update per iteration, K |
| `maxIterations` | `1000` | Iteration cap per time step |
| `showPlots` | `true` | Enable figures |
| `plotEvery` | `2` | Steps between plot updates |

Pass any of these fields in an options struct. `finalTime` must be an integer multiple of `timeStep`. The patch uses node indices `round(n/3):round(2*n/3)` along x and y, as in the original code; its exact extent changes with the grid.

Temperature is stored as `temperature(ix, iy, iz)`; coordinate vectors are returned in meters. Plots use millimeters and kelvin with a common color scale across depths.

## Physical interpretation and limitations

For anisotropic conduction without an internal heat source, the target continuum equation is:

```text
rho * cp * dT/dt = div(K * grad(T))
```

The implementation is a simplified approximation with the following limitations:

1. **Prescribed temperature, not a laser model.** A patch on `z = 0` is held at 800 K by default. Beam power, absorption, beam profile, and heat flux are not modeled.
2. **Incomplete conductivity tensor.** Only `kx`, `ky`, and `kz` are retained. The off-diagonal component produced by rotating an anisotropic ply is omitted. The +45 and -45 degree layers therefore have identical represented conductivities.
3. **Variable-conductivity discretization.** Each node's conductivity is used for both neighboring coefficients in a direction. Face conductivities and conservative interface fluxes are not implemented; the full variable-coefficient divergence operator is not represented.
4. **Lagged boundary treatment.** Intended insulated faces are implemented by copying adjacent temperatures after the interior solve, then restoring the hot patch. The zero-flux condition is not enforced within the iterative system.
5. **Lagged material properties.** Properties are updated once per time step using the preceding temperature field. There is no nonlinear iteration within a time step.
6. **Illustrative properties.** The original polynomial/linear property formulas have no supplied experimental source or validated temperature range. In particular, physical suitability at 800 K is not established; degradation and phase changes are absent. Positive computed properties alone do not establish validity.
7. **Layer assignment.** Layers are assigned in node groups using `ceil(nz/4)`, giving 4, 4, 4, and 3 nodes on the default grid. This is not a precise equal-thickness laminate discretization. Very small grids may not represent all four orientations.
8. **Convergence versus accuracy.** The stopping condition measures temperature changes between iterations, not the equation residual or physical error. Runs reaching the iteration cap are flagged. No mesh/time-step independence, energy-balance validation, or reference-solution comparison has been performed.

## Reading the returned data

`results.history` contains the physical time, iteration count, maximum final interior update, convergence flag, relaxation factor used in the final iteration, and solver time for each step. These diagnostics refer to the interior solve before the boundary copying operation.

`results.solverSeconds` sums solver-only timing. `results.elapsedSeconds` measures the time loop including any in-loop plotting, but excludes final depth plots. For timing comparisons, disable plotting and use identical options; repeat runs to account for MATLAB warm-up and timing variation. A convergence flag does not establish numerical or physical accuracy.

## Next steps

- Implement the full rotated conductivity tensor and conservative interface fluxes.
- Incorporate insulated boundary conditions into the solve.
- Add residual and energy-balance diagnostics.
- Check mesh and time-step independence with a consistently defined physical patch and laminate.
- Compare with a known solution and use sourced material data.

## Repository notes

This function-based version was prepared from the original `Ultimate_Composite_Project.m` supplied for review. An unchanged copy is included in `original/Ultimate_Composite_Project.m` for reference. The recommended entry point is `run_simulation.m`. It retains the original property formulas, stencil, stacking assignment, boundary-update order, and adaptive SOR heuristic. Interactive numeric prompts were replaced with named options; input checks, diagnostics, and consistent plot labeling were added. The original file used the term "laser" for the temperature patch; that terminology is corrected here.

See `PUBLISH.md` for the steps to initialize this folder as a Git repository and publish it on GitHub.

No license has been selected. Add a license of your choice before granting reuse permissions; public GitHub hosting alone does not provide an open-source license.
