"""Reduced-order electroosmotic (EO) assistance of liquid transport in a micropillar wick.

Passive wick (Case A) versus the identical wick with an applied in-plane electric field (Case B).
Intrinsic permeability K and the electrical formation factor F are geometric properties and are
the same in both cases; the field adds an independent electroosmotic flux. Equations, assumptions
and validity limits: docs/electroosmosis_model.md.

Modules
    physics          Helmholtz-Smoluchowski mobility, double-layer scales and validity diagnostics
    formation_factor electrical formation factor F of the pillar array (analytical, numerical, COMSOL)
    wick1d           Case A / Case B transport, forced-flow Q(dP), 1D capillary dryout with EO
    inputs           provenance-checked loader for parameters/electroosmosis/eo_inputs.yaml
"""
