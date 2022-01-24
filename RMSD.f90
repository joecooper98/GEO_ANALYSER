!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
! RMSD CACULATION SCRIPT                        !
! Joe Cooper 2022                               !
!                                               !
! MUST BE COMPILED WITH LAPACK LIBRARY, e.g.    !
! $ gfortran RMSD.f90 -llapack -o RMSD.o        !
!                                               !
! USE KABSCH ALGORITHM WITH SVD                 !
! GEOMETRIES ARE .xyz FILES PASSED THROUGH COM- !
! MAND LINE, e.g.                               !
! $ ./RMSD.o geom1.xyz geom2.xyz                !
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program RMSD
   implicit none

   DOUBLE PRECISION, dimension(:, :), ALLOCATABLE :: GEOM1, GEOM2
   DOUBLE PRECISION :: val
   CHARACTER(LEN=50) :: file1, file2
   integer :: f1 = 1, f2 = 2, io1, io2
   CHARACTER(LEN=2), dimension(:), ALLOCATABLE :: atoms, atoms2
   LOGICAL :: HBOOL = .false. ! HBOOL is true if you want to remove Hydrogens

   IF (COMMAND_ARGUMENT_COUNT() .NE. 2) THEN
      WRITE (*, *) 'TWO INPUTS REQUIRED, ABORTING...'
      STOP
   END IF

   !Get command line arguments for filename
   CALL GET_COMMAND_ARGUMENT(1, file1) ! file to be checked
   CALL GET_COMMAND_ARGUMENT(2, file2) ! reference geometry
   
   !Read .xyz files using submodule, save atoms
   if (HBOOL) then
      CALL read_xyz_no_h(GEOM2, file2, f2, atoms2,io2)
   else
      CALL read_xyz(GEOM2, file2, f2, atoms2,io2)
   end if

   do 
   if (HBOOL) then
      CALL read_xyz_no_h(GEOM1, file1, f1, atoms, io1)
   else
      CALL read_xyz(GEOM1, file1, f1, atoms, io1)
   end if
   if (io1 .LT. 0) then
           exit
   end if 

   ! Translate centroid to origin of coordinates
   CALL center(GEOM1)
   CALL center(GEOM2)

   ! Align using Kabsch SVD algorithm
   CALL align(GEOM1, GEOM2)

   ! Calculate RMSD, assign to val
   CALL calc_rmsd(GEOM1, GEOM2, val)
   !print val
   print '(ES12.5)', val
   DEALLOCATE(GEOM1)
   DEALLOCATE(atoms)
   end do
contains

   subroutine center(geom)
      implicit none

      DOUBLE PRECISION, dimension(:, :):: GEOM
      DOUBLE PRECISION, dimension(3) :: CENTROID
      INTEGER :: noatoms, i

      noatoms = SIZE(GEOM, DIM=1)
      CENTROID = SUM(GEOM, DIM=1)/noatoms

      do concurrent (i = 1:noatoms)
         GEOM(i, :) = GEOM(i, :) - CENTROID(:)
      end do
   end subroutine

   subroutine calc_rmsd(G1, G2, val)
      implicit none

      DOUBLE PRECISION, dimension(:, :):: G1, G2
      DOUBLE PRECISION :: val

      if (SIZE(G1, DIM=1) .NE. SIZE(G2, DIM=1)) then
         WRITE (*, *) 'TWO INPUTS HAVE DIFFERENT NUMBERS OF ATOMS, ABORTING...'
         STOP
      end if

      val = sqrt(sum((G1 - G2)**2)/(size(G1, dim=1)))

   end subroutine

   subroutine read_xyz(GEOM, filename, fid, atoms, io)
      implicit none

      DOUBLE PRECISION, dimension(:, :), ALLOCATABLE :: GEOM
      CHARACTER(len=50):: filename
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: atoms
      integer :: noatoms, i, fid, io

      OPEN (fid, file=TRIM(filename))
      READ (fid, *,IOSTAT=io) noatoms 
      !This catches the end of a file to stop the calculation
      if (io .LT. 0) then
              return
      end if 
      READ (fid, *)

      allocate (GEOM(noatoms, 3))
      allocate (atoms(noatoms))

      DO i = 1, noatoms
         READ (fid, *,IOSTAT=io) atoms(i), GEOM(i, :)
      end do

   end subroutine

   subroutine read_xyz_no_h(GEOM, filename, fid, atoms, io)
      implicit none

      DOUBLE PRECISION, dimension(:, :), ALLOCATABLE :: intgeom, GEOM
      CHARACTER(len=50):: filename
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: atoms
      integer :: noatoms, no_non_h, i, fid, counter, io

      OPEN (fid, file=TRIM(filename))
      READ (fid, *,IOSTAT=io) noatoms
      !This catches the end of a file to stop the calculation
      if (io .LT. 0) then
              return
      end if 
      READ (fid, *)

      allocate (intgeom(noatoms, 3))
      allocate (atoms(noatoms))

      no_non_h = 0
      DO i = 1, noatoms
         READ (fid, *, IOSTAT=io) atoms(i), intgeom(i, :)
         if (.not. ((trim(atoms(i)) .EQ. 'H') .OR. (trim(atoms(i)) .EQ. 'h '))) then
            no_non_h = no_non_h + 1
         end if
      end do

      allocate (GEOM(no_non_h, 3))
      counter = 0
      do i = 1, noatoms
         if (.not. ((trim(atoms(i)) .EQ. 'H ') .OR. (trim(atoms(i)) .EQ. 'h '))) then
            counter = counter + 1
            GEOM(counter, :) = intgeom(i, :)
         end if
      end do

      DEALLOCATE (intgeom)
   end subroutine

   subroutine write_xyz(GEOM, filename, fid, atoms)
      implicit none

      DOUBLE PRECISION, dimension(:, :):: GEOM
      CHARACTER(LEN=9):: filename
      CHARACTER(LEN=2), DIMENSION(:):: atoms
      integer :: noatoms, i, fid

      noatoms = size(GEOM, dim=1)

      OPEN (fid, file=filename)
      write (fid, *) noatoms
      write (fid, *) filename

      DO i = 1, noatoms
         write (fid, *) atoms(i), GEOM(i, :)
      end do
   end subroutine

   subroutine align(GEOM1, GEOM2)
      implicit none
      DOUBLE PRECISION, dimension(:, :):: GEOM1, GEOM2
      DOUBLE PRECISION, dimension(3, 3) :: covar, ROTMAT, U, VT, sing, dmat
      DOUBLE PRECISION, dimension(20) :: work
      integer :: info, i

      dmat(:, :) = 0.
      dmat(1, 1) = 1.
      dmat(2, 2) = 1.

      covar = MATMUL(TRANSPOSE(GEOM1), GEOM2)

      CALL DGESVD('A', 'A', 3, 3, covar, 3, sing, U, 3, VT, 3, work, 20, info)

      if (dettt(matmul(TRANSPOSE(VT), TRANSPOSE(U))) .LT. 0) then
         dmat(3, 3) = -1.
      else
         dmat(3, 3) = 1.
      end if

      ROTMAT = matmul(TRANSPOSE(VT), matmul(dmat, TRANSPOSE(U)))

      do concurrent (i = 1:SIZE(GEOM1, DIM=1))
         GEOM1(i, :) = matmul(ROTMAT, GEOM1(i, :))
      end do

   end subroutine

   DOUBLE PRECISION function dettt(MAT)
      implicit none
      DOUBLE PRECISION, dimension(3, 3) :: MAT

      dettt = MAT(1, 1)*(MAT(2, 2)*MAT(3, 3) - MAT(2, 3)*MAT(3, 2))
      dettt = dettt + -MAT(1, 2)*(MAT(2, 1)*MAT(3, 3) - MAT(2, 3)*MAT(3, 1))
      dettt = dettt + MAT(1, 3)*(MAT(2, 1)*MAT(3, 2) - MAT(2, 1)*MAT(3, 2))
   end function
end program
