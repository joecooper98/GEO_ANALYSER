!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Geometry analyser script                        !
! Written by Joe Cooper, 2022                     !
! Takes in .xyz files and performs analyses       !
! Compile with something like                     !
! $ gfortran -lblas -llapack -O3 GEO_ANALYSER.f90 !
! Run with something like                         !
! $ ./GEO_ANALYSER input traj.xyz ref.xyz         !
! Feel free to change to your liking              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program geo_analyser
  implicit none

  DOUBLE PRECISION, dimension(:, :), ALLOCATABLE :: GEOM1, GEOM2
  DOUBLE PRECISION, dimension(:), allocatable  :: val
  double precision :: time
  CHARACTER(LEN=50) :: inpfile, file1, file2
  CHARACTER(LEN=50), allocatable :: str(:)
  integer :: f1 = 1, f2 = 2, finp = 3, io1, io2, no_operations, i
  CHARACTER(LEN=2), dimension(:), ALLOCATABLE :: atoms, atoms2
  LOGICAL :: HBOOL = .false. ! HBOOL is true if you want to remove Hydrogens



  ! Get command line arguments. If it's more than 2, then we use
  ! the final input as a reference geometry for e.g. RMSD
  ! Ideally this will be generalised to as many RMSD geometries as you like
  ! But this is not trivial in fortran - would require a restructuring
  ! to only read the geometries when needed

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    WRITE (*, *) 'AT LEAST TWO INPUTS REQUIRED, ABORTING...'
    STOP
  else if (COMMAND_ARGUMENT_COUNT() .EQ. 2) THEN
    CALL GET_COMMAND_ARGUMENT(1, inpfile) ! file to be checked
    CALL GET_COMMAND_ARGUMENT(2, file1) ! input geometry
    CALL read_xyz(GEOM2, file1, f2, atoms2,io2,time) ! this is just here to simplify later on
    ! it's not at all used.
  else 
    CALL GET_COMMAND_ARGUMENT(1, inpfile) ! file to be checked
    CALL GET_COMMAND_ARGUMENT(2, file1) ! input geometry
    CALL GET_COMMAND_ARGUMENT(3, file2) ! reference geometry
    CALL read_xyz(GEOM2, file2, f2, atoms2,io2,time)
  end if 

  ! read the input file - will return a array of characters deciding the operations
  call read_input(inpfile, finp, str, no_operations)

  ! this will be where the values for each of the operations is held
  allocate(val(no_operations))


  do  ! keep looping until end of file is hit
  ! this will loop over every geometry in a trajectory

  ! read the geometry for this time step
  ! overwrites the geometry from the last time step
  CALL read_xyz(GEOM1, file1, f1, atoms, io1,time)

  if (io1 .LT. 0) then
          exit ! will exit if the reader hits the end of the file
  end if 

  do i = 1,no_operations ! loop over all of the operations
  ! call the main flow controller - this will read the string, work out what to do, and then do it
  ! TO DO improve this by reading through once, and then assigning an integer array to help
  call decider(trim(str(i)), GEOM1, GEOM2, val(i),atoms)
  end do
  ! writes the values to the main out - then pipe to file.
  call printer(val, no_operations)
  end do

  ! deallocate everything at the end.
  if (allocated(GEOM2)) then
          deallocate(GEOM2)
          deallocate(atoms2)
  end if
  if (allocated(val)) then

  deallocate(val)
  end if
  if (allocated(GEOM1)) then
  DEALLOCATE(GEOM1)
  end if
  if (allocated(atoms)) then
  DEALLOCATE(atoms)
  end if
contains

  subroutine decider(str, GEOM1, GEOM2, val, atoms)
    CHARACTER(LEN=*), intent(IN) ::  str
    DOUBLE PRECISION, intent(IN) :: GEOM1(:,:), GEOM2(:,:)
    DOUBLE PRECISION, intent(OUT) :: val
    character(len=2), intent(in) :: atoms(:)
    character(len=5) :: parse
    integer, allocatable :: indices(:)


    call string_split(str,parse, indices)

    select case (StrUpCase(trim(parse)))

    case ('RMSD')
      call RMSD(GEOM1, GEOM2, val)
    case ('RNOH')
      call RNOH(GEOM1, GEOM2, val,atoms)
    case ('DIST')
      val = DIST(GEOM1,indices)
    case ('ANGL')
      val = ANG(GEOM1,indices)
    case ('DIHE')
      val = DIHE(GEOM1,indices)
    case default
      write(*,*) 'Not implemented yet!'
      stop
    end select
    if (allocated(indices)) then
    deallocate(indices)
    end if
  end subroutine

subroutine printer(x,n)
                double precision, intent(in) :: x(:)
                integer,intent(in) :: n

                ! write(*,'(ES16.8,2x)',advance='no') t
                do i = 1,n
                if (i .eq. n) then
                        write(*,'(ES16.8)') x(i)
                else
                        write(*,'(ES16.8,2x)',advance='no') x(i)
                end if 
                end do
                end subroutine

  subroutine RMSD(GEOM1, GEOM2, val)
    DOUBLE PRECISION, intent(IN) :: GEOM1(:,:), GEOM2(:,:)
    DOUBLE PRECISION, intent(OUT) :: val
    ! Translate centroid to origin of coordinates
    CALL center(GEOM1)
    CALL center(GEOM2)
    !
    ! Align using Kabsch SVD algorithm
    CALL align(GEOM1, GEOM2)
    ! Calculate RMSD, assign to val
    CALL calc_rmsd(GEOM1, GEOM2, val)
  end subroutine
  
  subroutine RNOH(GEOM1, GEOM2, val, atoms)
    DOUBLE PRECISION, intent(IN) :: GEOM1(:,:), GEOM2(:,:)
    double precision, allocatable :: intgeom1(:,:), intgeom2(:,:)
    DOUBLE PRECISION, intent(OUT) :: val
    character(len=2), intent(in) :: atoms(:)
    logical :: isH(size(GEOM1,dim=1))
    integer :: i, counter=0
    
    do i = 1, size(atoms)
      if (atoms(i) .ne. 'H ') then
              isH(i) = .true.
      else
              isH(i) = .false.
      end if
      end do

    counter = count(isH)

    ! if (.not. allocated(intgeom1)) then
            !no point allocating what's already allocated
            allocate(intgeom1(counter,3))
            allocate(intgeom2(counter,3))
    ! end if
    
    do i = 1, size(atoms)
      if (isH(i)) then
              intgeom1(i,:) = GEOM1(i,:)
              intgeom2(i,:) = GEOM2(i,:)
      end if
      end do
            
    ! Translate centroid to origin of coordinates
    CALL center(intgeom1)
    CALL center(intgeom2)
    !
    ! Align using Kabsch SVD algorithm
    CALL align(intgeom1, intgeom2)
    ! Calculate RMSD, assign to val
    CALL calc_rmsd(intgeom1, intGEOM2, val)
  end subroutine

  double precision function DIST(GEOM, indices)
    double precision, intent(IN) :: GEOM(:,:)
    integer, dimension(2), intent(IN) :: indices

    DIST = norm(GEOM(indices(1),:)-GEOM(indices(2),:))
  end function

  double precision function ANG(GEOM, indices)
    double precision, intent(IN) :: GEOM(:,:)
    integer, dimension(3), intent(IN) :: indices
    double precision :: vec1(3), vec2(3)

    vec1 = GEOM(indices(1),:) - GEOM(indices(2),:)
    vec2 = GEOM(indices(3),:) - GEOM(indices(2),:)

    ANG = rad2deg(ACOS(DOT_PRODUCT(vec1,vec2)/(norm(vec1)*norm(vec2))))
  end function

  double precision function DIHE(GEOM, indices)
    double precision, external :: ddot
    double precision, intent(IN) :: GEOM(:,:)
    integer, dimension(4), intent(IN) :: indices
    double precision, dimension(3) :: p0,p1,p2,p3,b0,b1,b2,v,w
    double precision :: x,y
    ! based on Praxeolitic formula from Stack overflow

    p0 = GEOM(indices(1),:)
    p1 = GEOM(indices(2),:)
    p2 = GEOM(indices(3),:)
    p3 = GEOM(indices(4),:)

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    !  normalize b1 so that it does not influence magnitude of vector
    !  rejections that come next
    b1 = b1/norm(b1)

    ! # vector rejections
    ! # v = projection of b0 onto plane perpendicular to b1
    ! #   = b0 minus component that aligns with b1
    ! # w = projection of b2 onto plane perpendicular to b1
    ! #   = b2 minus component that aligns with b1
    v = b0 - ddot(3,b0,1,b1,1)*b1
    w = b2 - ddot(3,b2,1,b1,1)*b1

   ! # angle between v and w in a plane is the torsion angle
   ! # v and w may not be normalized but that's fine since tan is y/x
    x = ddot(3,v,1,w,1)
    y = ddot(3,cross_product(b1, v),1,w,1)
    DIHE=rad2deg(atan2(y, x))


  end function

  subroutine read_input(filename, fid, str, no_operations)
    CHARACTER(LEN=*), intent(IN) :: filename
    CHARACTER(LEN=50), dimension(:), allocatable, intent(out) :: str
    integer, intent(in) :: fid
    integer, intent(out) :: no_operations
    integer :: io,  i

    no_operations = 0 
    open(fid, file=trim(filename))
    do
    read(fid, *, iostat=io)
    if (io .lt. 0) then
      exit
    end if
    no_operations = no_operations + 1
    end do 

    allocate(str(no_operations))
    rewind(fid)
    do i = 1,no_operations
    read(fid,'(A)') str(i)
    end do 

  end subroutine


  subroutine string_split(str, output, indices)
    character(len=*), intent(IN) :: str
    character(len=5), intent(out) :: output
    integer, allocatable :: indices(:)
    integer :: noinputs
    character(len=45) :: intstring

    read(str,'(A5,A)') output, intstring

    noinputs = count_spaces(trim(intstring))
    if (noinputs .ne. 0) then
      allocate(indices(noinputs+1))
      read(intstring,*) indices(:)
    else
            allocate(indices(1))
!otherwise it gets mad
    end if 

  end subroutine 

  integer function count_spaces(str) result(n)
          character(len=*):: str
          integer :: i

          n=0
          do i = 1,len_trim(str)
            if (str(i:i) .eq. ' ') then
                    n = n+1
            end if
            end do
            end function




  double precision pure function rad2deg(x) result(y)
          double precision, intent(in) :: x
          double precision, parameter :: conv =57.29577951308

          y = x*conv
          end function

  double precision function norm(X)
    double precision, external ::dnrm2
    double precision, intent(IN) :: x(3)

    norm = dnrm2(3,x,1)
  end function

  function cross_product(x,y)
          double precision, intent(IN) :: x(3), y(3)
          double precision :: ssymvec(3,3)
          double precision :: cross_product(3)

          ssymvec(:,:) = 0.d0
          ssymvec(2,1) = x(3)
          ssymvec(1,2) =-x(3)
          ssymvec(3,1) =-x(2)
          ssymvec(1,3) = x(2)
          ssymvec(2,3) =-x(1)
          ssymvec(3,2) = x(1)

        call dgemv('n',3,3,1.d0,ssymvec,3,y,1,0.,cross_product,1)

        end function

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

  subroutine read_xyz(GEOM, filename, fid, atoms, io, time)
    implicit none

    DOUBLE PRECISION,  ALLOCATABLE, intent(OUT) :: GEOM(:,:)
    double precision :: time
    CHARACTER(len=50):: intstring1, intstring2
    CHARACTER(len=50),intent(IN):: filename
    CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE , intent(OUT):: atoms
    integer :: noatoms, i, fid, io

    OPEN (fid, file=TRIM(filename))
    READ (fid, *,IOSTAT=io) noatoms 
    !This catches the end of a file to stop the calculation
    if (io .LT. 0) then
      return
    end if 
    ! READ (fid,*) intstring1, time
    read(fid,*)

    if (.not. allocated(GEOM)) then
    allocate (GEOM(noatoms, 3))
    allocate (atoms(noatoms))
    end if

    DO i = 1, noatoms
    READ (fid, *,IOSTAT=io) atoms(i), GEOM(i, :)
    end do

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

    do i = 1,SIZE(GEOM1, DIM=1)
    call dgemv('n',3,3,1.d0,ROTMAT,3,GEOM1(i,:),1,0.,GEOM1(i,:),1)
    end do

  end subroutine

  DOUBLE PRECISION pure function dettt(MAT)
    implicit none
    DOUBLE PRECISION, intent(in), dimension(3, 3) :: MAT

    dettt = MAT(1, 1)*(MAT(2, 2)*MAT(3, 3) - MAT(2, 3)*MAT(3, 2))
    dettt = dettt - MAT(1, 2)*(MAT(2, 1)*MAT(3, 3) - MAT(2, 3)*MAT(3, 1))
    dettt = dettt + MAT(1, 3)*(MAT(2, 1)*MAT(3, 2) - MAT(2, 1)*MAT(3, 2))
  end function

  FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String )
    ! -- Argument and result
    CHARACTER( * ), INTENT( IN )     :: Input_String
    CHARACTER( LEN( Input_String ) ) :: Output_String
    ! -- Local variables
    INTEGER :: i, n
    CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    ! -- Copy input string
    Output_String = Input_String
    ! -- Loop over string elements
    DO i = 1, LEN( Output_String )
    ! -- Find location of letter in lower case constant string
    n = INDEX( LOWER_CASE, Output_String( i:i ) )
    ! -- If current substring is a lower case letter, make it upper case
    IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
    END DO 
  end function

end program
