# GEO_ANALYSER

Simple program to calculate geometric quantities across trajectories.

I've only included a fairly minimal feature set so far, this can

- Read .xyz files
- Calculate
   - Distances between two atoms
   - Angles
   - Dihedrals 
   - RMSD (Kabsch Algorithm)
   - RMSD ignoring the Hydrogens
- Write .xyz files


Compile using something like
```
gfortran GEO_ANALYSER.f90 -llapack -lblas -o GEO_ANALYSER.o
```
or 
```
ifort GEO_ANALYSER.f90 -llapack -lblas -o GEO_ANALYSER.o
```
And run it like
```
./GEO_ANALYSER.o input comp_geom.xyz ref_geom.xyz
```
If you want to be able to run it anywhere, just place it in your path, and forego the `./`.

The input file is just a plain text file with a specific format --- best shown by an example

```
DIST 1 2
ANGL 3 4 8
DIHE 1 2 11 12
RMSD
RNOH
```

This will (in order) compute the


1. The time (given by the second value in the comment line)
2. The distance between atoms 1 and 2 in the `comp_geom.xyz` file
3. The bond angle between atoms 3, 4, and 8 in the `comp_geom.xyz` file
4. The dihedral angle between atoms 1, 2, 11 and 12 in the `comp_geom.xyz` file 
5. The RMSD between `comp_geom.xyz` and `ref_geom.xyz`
6. The RMSD between `comp_geom.xyz` and `ref_geom.xyz`, ignoring the contributions from atoms labelled `H`

If you want to add your own, then it's fairly simple. Give it a 4 letter name, create a subroutine/function that returns a single value, and then add it to the switch statement in the `decider` routine. 

## Trajectories

If the `comp_geom.xyz` file is a `.xyz` trajectory file (i.e. a concatenated series of `.xyz` files), then the program will read each of the geometries in turn, printing the values for each of them in order.


## Optimised LAPACK

For better speed, it is recommended you utilise optmised LAPACK libraries. I recommend MKL, which is from INTEL. Utilise their documentation to work out exactly how to compile, but for me it is something like

```
ifort RMSD.f90 -fast -O3 -g ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl  -i8  -I"${MKLROOT}/include" -o RMSD.o 
```
but this depends on your system architecture.
