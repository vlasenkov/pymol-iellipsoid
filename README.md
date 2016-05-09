# pymol-iellipsoid
Pymol plugin allowing to evaluate and draw inertia ellipsoid for a selection.
Useful to demonstrate mobility of atoms.

Interface:
* `ie_build` build a single ellipsoid for a selection
* `ie_build_all` build ellipsoid for each atom ID
* `ie_build_file` load PDB, align, and apply `ie_build_all` to it
