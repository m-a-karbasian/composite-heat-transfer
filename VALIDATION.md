# Verification status

Prepared on 2026-09-30.

## Execution status

MATLAB R2025a was found on the preparation machine. Two batch launches were attempted, including a retry with a workspace-local preferences directory. Both failed before executing MATLAB code with:

```text
Fatal Startup Error:
System Error: File system inconsistency
```

Therefore, neither the behavioral checks, full default simulation, graphics checks, nor MATLAB Code Analyzer checks are confirmed to pass. No runtime, speedup, or simulation-result claims are made for this packaged version. The source and documentation were reviewed statically.

## Included checks

From the repository folder in a working MATLAB session, run:

```matlab
check_simulation
```

This checks preservation of a uniform-temperature state, agreement of the two solvers on the same small discrete problem for both materials, finite/bounded temperatures for that case, convergence flags, and rejection of a non-integral time grid. These are implementation checks, not physical validation.

For a full default run without plots:

```matlab
r = run_simulation(struct('showPlots', false));
assert(all(r.history.Converged));
```

For a short graphics check with unequal axis node counts:

```matlab
r = run_simulation(struct('gridSize', [7 8 9], 'finalTime', 20));
```

Physical validation, grid/time-step studies, and numerical comparison against the original executable script remain outstanding. The original script was subsequently supplied from its new location and reviewed alongside the refactored version. An unchanged copy is included in `original/Ultimate_Composite_Project.m`; its SHA-256 was checked against the supplied file at packaging time:

```text
B1F8A37275B55B1E75636ACDC064368C1D0D62EB0243CF222C12334823859F99
```
