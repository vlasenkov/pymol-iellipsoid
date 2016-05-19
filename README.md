# pymol-iellipsoid
Pymol plugin allowing to evaluate and draw inertia ellipsoid for a selection.
Useful to demonstrate mobility of atoms.

Interface:
* `ie_build` build a single ellipsoid for a selection
* `ie_build_all` build ellipsoid for each atom ID
* `ie_build_file` load PDB, align, and apply `ie_build_all` to it

Usage example (after installation of the plugin):
```
ie_build_file lipid.pdb, scale=3, col=[0.1, 0.7, 0.9]
```