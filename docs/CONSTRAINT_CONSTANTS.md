## NanoAmorphizer Constraint Constants

This document records the geometric constraint constants used in NanoAmorphizer v1.0.0 to ensure reproducibility and to preserve their current status for future development.

### Current constants

**DEFAULT_MIN_DISTANCE_SCALE = 0.70**  
Minimum allowed interatomic separation, expressed as a scale of the ASE covalent-radii-based distance criterion.

**DEFAULT_PARENT_MIN_SCALE = 0.85**  
Lower bound for the distance of a parent first-neighbor pair relative to its original distance.

**DEFAULT_PARENT_MAX_SCALE = 1.20**  
Upper bound for the distance of a parent first-neighbor pair relative to its original distance.

**DEFAULT_NEIGHBOR_CUTOFF_SCALE = 1.20**  
Scale factor used with ASE covalent radii to define the geometric first-neighbor network.

**DEFAULT_MAX_ANGLE_DEVIATION = 20.0°**  
Maximum permitted deviation of a parent local angle from its original value.

**DEFAULT_MAX_ATTEMPTS = 50**  
Maximum number of stochastic displacement proposals tested for one shell atom before retaining its original position.

### Status

These values are **provisional geometric defaults**. They currently provide stable and structurally coherent model generation, but they have not yet been established as universal physical limits or optimized parameters.

They should therefore not be interpreted as material-independent bond limits, energetic criteria, or physically validated amorphization thresholds.

### Future improvement

Future NanoAmorphizer versions may evaluate and, where appropriate, refine these defaults using justified or optimized criteria. Possible approaches include:

- crystallographic and materials-science literature;
- distributions of distances and angles derived from the parent structure;
- material- or pair-specific tolerances;
- sensitivity analysis over the constraint parameters;
- comparison with relaxed or experimentally supported structures;
- optimization of the balance between structural disorder and local geometric coherence.

The constants used in NanoAmorphizer v1.0.0 should remain documented to preserve reproducibility as future calibration is introduced.