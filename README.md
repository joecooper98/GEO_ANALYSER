# RMSD

Simple program to calculate RMSD between two .xyz files.

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
./RMSD.o comp_geom.xyz ref_geom.xyz
```
If you want to be able to run it anywhere, just place it in your path, and forego the `./`.

This will print the minimised RMSD for `comp_geom.xyz` compared to `ref_geom.xyz`. The geometry `comp_geom.xyz`, can be a multi-geometry .xyz file (i.e. concatenated .xyz files), and the reference geometry will only read the first file.

If you want to only calculate the heavy atoms (i.e. no Hydrogens), then simply change the line
```
   LOGICAL :: HBOOL = .false. ! HBOOL is true if you want to remove Hydrogens
```
to 
```
   LOGICAL :: HBOOL = .true. ! HBOOL is true if you want to remove Hydrogens
```
and recompile (e.g. to `RMSD_NO_H.o`)


## Optimised LAPACK

For better speed, it is recommended you utilise optmised LAPACK libraries. I recommend MKL, which is from INTEL. Utilise their documentation to work out exactly how to compile, but for me it is something like

```
ifort RMSD.f90 -fast -O3 -g ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl  -i8  -I"${MKLROOT}/include" -o RMSD.o 
```
but this depends on your system architecture.
