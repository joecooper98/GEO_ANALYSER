### RMSD

Simple program to calculate RMSD between two .xyz files.

Written in FORTRAN 90 for speed.

I've only included a fairly minimal feature set so far, this can

- Read .xyz files
- Translate the centroid to the origin
- Rotate the geometries to minimise RMSD (Kabsch algorithm using LAPACK SVD)
- Calculate RMSD
- Write .xyz files
- Choose to ignore the hydrogens

Compile using something like
```
gfortran RMSD.f90 -llapack -o RMSD.o
```
And run it like
```
./RMSD.o geom1.xyz geom2.xyz
```
If you want to be able to run it anywhere, just place it in your path, and forego the `./`.

This will print the minimised RMSD for the two geometries.

If you want to only calculate the heavy atoms (i.e. no Hydrogens), then simply change the line
```
   LOGICAL :: HBOOL = .false. ! HBOOL is true if you want to remove Hydrogens
```
to 
```
   LOGICAL :: HBOOL = .true. ! HBOOL is true if you want to remove Hydrogens
```
and recompile (e.g. to `RMSD_NO_H.o`)
