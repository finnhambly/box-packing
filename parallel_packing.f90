!------------------------------------------------------------------------------
!
!   HPC Applications Final Assignment
!
!   (c) Finn Hambly 2020. All rights reserved.
!   File created: 2020-04-18. File last modified: 2020-04-26
!
!   Email: finn.hambly@york.ac.uk
!
!-------------------------------------------------------------------------------
!
! This program simulates ...
!
!
!-------------------------------------------------------------------------------
!   Simple serial implementation of Jason Blevin's
!   FORTRAN Random Number Generator
!
!   https://jblevins.org/log/openmp
!
!-------------------------------------------------------------------------------

module rng
  implicit none

  private
  public :: rng_t, rng_seed, rng_uniform

  ! Dimension of the state
  integer, parameter :: ns = 4

  ! Default seed vector
  integer, parameter, dimension(ns) :: default_seed &
       = (/ 521288629, 362436069, 16163801, 1131199299 /)

  ! A data type for storing the state of the RNG
  type :: rng_t
     integer, dimension(ns) :: state = default_seed
  end type rng_t

contains

  ! Seeds the RNG using a single integer and a default seed vector.
  subroutine rng_seed(self, seed)
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed
    self%state(1) = seed
    self%state(2:ns) = default_seed(2:ns)
  end subroutine rng_seed

  ! Draws a uniform real number on [0,1].
  function rng_uniform(self) result(u)
    type(rng_t), intent(inout) :: self
    real :: u
    integer :: imz

    imz = self%state(1) - self%state(3)

    if (imz < 0) imz = imz + 2147483579

    self%state(1) = self%state(2)
    self%state(2) = self%state(3)
    self%state(3) = imz
    self%state(4) = 69069 * self%state(4) + 1013904243
    imz = imz + self%state(4)
    u = 0.5d0 + 0.23283064d-9 * imz
  end function rng_uniform

end module rng

!-------------------------------------------------------------------------------
!   Main packing factor calculating program
!   (c) Finn Hambly 2020. All rights reserved.
!   File created: 2020-04-18. File last modified: 2020-04-26
!
!   Email: finn.hambly@york.ac.uk
!
!-------------------------------------------------------------------------------
program packing
  use rng
  use mpi
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300), L = 500
  real(kind=dp), parameter :: r = 1.234_dp, pi = 3.14159265358979_dp
  real(kind=dp), allocatable :: box(:,:,:)
  logical, allocatable       :: occupied(:,:)
  real(kind=dp) :: x, y, D
  logical :: pack
  integer :: i, j, k, N_circ, ios, ierr
  type(rng_t):: rng   ! custom type for random number generator
  ! MPI variables
  integer :: size, rank, x_div, y_div, x_max, y_max
  integer :: reqL1,reqL2,reqR1,reqR2,reqB1,reqB2,reqT1,reqT2
  integer :: sreqL1,sreqL2,sreqR1,sreqR2,sreqB1,sreqB2,sreqT1,sreqT2
  integer, dimension(1:MPI_STATUS_SIZE) :: sL1status,sL2status,sR1status, &
    sR2status, sT1status, sT2status, sB1status, sB2status, &
    rL1status,rL2status,rR1status, rR2status, rT1status, rT2status, rB1status,&
    rB2status

  ! Initialise MPI
  call MPI_init(ierr)
  if (ierr/=0) stop 'Error with MPI_init'
  call MPI_comm_rank(MPI_comm_world,rank,ierr)
  if (ierr/=0) stop 'Error with MPI_comm_rank'
  call MPI_comm_size(MPI_comm_world,size,ierr)
  if (ierr/=0) stop 'Error with MPI_comm_size'

  ! Open file for visualising placement
  if ( rank == 0 ) then
    open(unit=11, file='segment.dat', status='replace', iostat=ios)
    if ( ios /= 0 ) stop "Error opening file segment.dat"
  endif

  1 format('Error: number of threads,',i4,', larger than array size',i8)
  if ( size > L**2 ) then
    if ( rank == 0 ) print 1, size, L**2
    call MPI_Finalize(ierr)
    stop
  endif

  ! Decompose problem according to the number of threads
  ! Divide box into squares, or near rectangles if not possible
  x_div = ceiling(sqrt(real(size)))
  if ( mod(size,x_div) == 0 ) then
    y_div = size/x_div ! 2D decomposition
  else ! Resort to 1D decomposition if 2D is not possible
    x_div = size
    y_div = 1
  endif

  ! Allocate number of grid points to each thread, adding leftovers to the
  ! far ends
  if ( mod(rank+1, x_div) /= 0 ) then
    x_max = L/x_div
  else
    x_max = L/x_div + L - x_div*(L/x_div)
  endif

  if ( rank+1 <= size - x_div ) then
    y_max = L/y_div
  else
    y_max = L/y_div + L - y_div*(L/y_div)
  endif

  ! Initialise array on each thread
  2 format('Array length on rank ',i4,' is ', i8, ' by ', i8)
  print 2, rank, x_max, y_max
  ! Allocate arrays for storing coordinates (and occupation states) with a halo
  ! of +3 on each side
  allocate(box(x_max+6,y_max+6,2), stat=ierr)
  if (ierr /= 0) print *, 'Array "box": Allocation request denied'
  allocate(occupied(x_max+6,y_max+6), stat=ierr)
  if (ierr /= 0) print *, 'Array "occupied": Allocation request denied'

  ! set occupation status of box as empty
  box = 0.0_dp
  occupied = .false.

  ! random number seed
  call rng_seed(rng, 362436069+3624*rank)
  N_circ = 0
  do i = 1, 10000000
    ! generate new random coordinates within the box, depending on rank
    if ((mod(rank+1, x_div) == 1) .and. (mod(rank+1, x_div) == 0)) then
      x = (x_max-2*r) * rng_uniform(rng) + r + 4 ! do not place over the edge
    else if (mod(rank+1, x_div) == 1) then ! if on the far left
      x = (x_max-r) * rng_uniform(rng) + r + 4 ! do not place over left edge
    else if (mod(rank+1, x_div) == 0) then ! if on the far right
      x = (x_max-r) * rng_uniform(rng) + 4 ! do not place over right edge
    else
      x = x_max * rng_uniform(rng) + 4 !otherwise, allow overlap with halo
    end if

    if ((rank+1 > size - x_div) .and. (rank+1 <= x_div)) then ! top and bottom:
      y = (y_max-2*r) * rng_uniform(rng) + r + 4 ! do not place over edge
    else if (rank+1 > size - x_div) then ! if at the top
      y = (y_max-r) * rng_uniform(rng) + 4 ! do not place over top edge
    else if (rank+1 <= x_div) then !if at the bottom
      y = (y_max-r) * rng_uniform(rng) + r + 4 ! do not place over bottom edge
    else
      y = y_max * rng_uniform(rng) + 4 ! otherwise, allow overlap with halo
    end if

    ! set boolean test to default value
    pack = .true.

    ! if grid point already occupied, then skip these coordinates
    if (occupied(int(x),int(y))) cycle

    ! Not exhaustive check but skips any that will definitely not fit
    ! Look and, if nearest neighbours are occupied, skip:
    if ( occupied(int(x)+1,int(y)) .or. occupied(int(x)-1,int(y)) .or. &
    occupied(int(x),int(y)+1) .or. occupied(int(x),int(y)-1)) cycle

    ! Update left halos
    if (mod(rank+1, x_div) /= 1) then
      if (int(x) < 7) then! Cordon off this point
        occupied(int(x), int(y)) = .true.
        box(int(x), int(y), 1) = x
        box(int(x), int(y), 2) = y
      endif

      ! Send coordinate data
      call MPI_Issend(box(4:6,1,1),6*(y_max+6), MPI_DOUBLE_PRECISION, &
      & rank-1, 1, MPI_comm_world, sreqL1, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqL1'
      ! Send occupation data
      call MPI_Issend(occupied(4:6,1), 3*(y_max+6), MPI_LOGICAL, &
      & rank-1, 2, MPI_comm_world, sreqL2, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqL2'

      ! Receive coordinate data from the left
      call MPI_Irecv(box(1:3,1,1), 6*(y_max+6), MPI_DOUBLE_PRECISION, &
      & rank-1, 3, MPI_comm_world, reqL1, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqL1'
      ! Receive occupation data from the left
      call MPI_Irecv(occupied(1:3,1), 3*(y_max+6), MPI_LOGICAL, &
      & rank-1, 4, MPI_comm_world, reqL2, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqL2'
    end if

    ! Update right halo
    if (mod(rank+1, x_div) /= 0) then
      ! Cordon off this point
      if (int(x) > x_max) then
        occupied(int(x), int(y)) = .true.
        box(int(x), int(y), 1) = x
        box(int(x), int(y), 2) = y
      endif

      ! Send coordinate data
      call MPI_Issend(box(x_max+1:x_max+3,1,1),6*(y_max+6), MPI_DOUBLE_PRECISION, &
      & rank+1, 3, MPI_comm_world, sreqR1, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqR1'
      ! Send occupation data
      call MPI_Issend(occupied(x_max+1:x_max+3,1), 3*(y_max+6), MPI_LOGICAL, &
      & rank+1, 4, MPI_comm_world, sreqR2, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqR2'

      ! Receive coordinate data from the right
      call MPI_Irecv(box(x_max+4:x_max+6,1,1), 6*(y_max+6), MPI_DOUBLE_PRECISION, &
      & rank+1, 1, MPI_comm_world, reqR1, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqR1'
      ! Receive occupation data from the right
      call MPI_Irecv(occupied(x_max+4:x_max+6,1), 3*(y_max+6), MPI_LOGICAL, &
      & rank+1, 2, MPI_comm_world, reqR2, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqR2'
    end if

    ! Update top halo
    if (rank+1 <= size - x_div) then
      ! Cordon off this point
      if (int(y) > y_max) then
        occupied(int(x), int(y)) = .true.
        box(int(x), int(y), 1) = x
        box(int(x), int(y), 2) = y
      endif

      ! Send coordinate data
      call MPI_Issend(box(:,y_max+1:y_max+3,1), 6*(x_max+6), MPI_DOUBLE_PRECISION, &
      & rank+x_div, 5, MPI_comm_world, sreqT1, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqT1'
      ! Send occupation data
      call MPI_Issend(occupied(:,y_max+1:y_max+3), 3*(x_max+6), MPI_LOGICAL, &
      & rank+x_div, 6, MPI_comm_world, sreqT2, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqT2'

      ! Receive coordinate data from above
      call MPI_Irecv(box(:,y_max+4:y_max+6,1), 6*(x_max+6), MPI_DOUBLE_PRECISION, &
      & rank+x_div, 7, MPI_comm_world, reqT1, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqT1'
      ! Receive occupation data from above
      call MPI_Irecv(occupied(:,y_max+4:y_max+6), 3*(x_max+6), MPI_LOGICAL, &
      & rank+x_div, 8, MPI_comm_world, reqT2, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqT2'
    end if

    ! Update bottom halos
    if (rank+1 > x_div)  then
      ! Cordon off this point
      if (int(y) < 7) then
        occupied(int(x), int(y)) = .true.
        box(int(x), int(y), 1) = x
        box(int(x), int(y), 2) = y
      end if

      ! Send coordinate data
      call MPI_Issend(box(:,4:6,1), 6*(x_max+6), MPI_DOUBLE_PRECISION, &
      & rank-x_div, 7, MPI_comm_world, sreqB1, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqB1'
      ! Send occupation data
      call MPI_Issend(occupied(:,4:6), 3*(x_max+6), MPI_LOGICAL, &
      & rank-x_div, 8, MPI_comm_world, sreqB2, ierr)
      if (ierr/=0) stop 'Error with MPI_Issend sreqB2'

      ! Receive coordinate data from below
      call MPI_Irecv(box(:,1:3,1), 6*(x_max+6), MPI_DOUBLE_PRECISION, &
      & rank-x_div, 5, MPI_comm_world, reqB1, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqB1'
      ! Receive occupation data from below
      call MPI_Irecv(occupied(:,1:3), 3*(x_max+6), MPI_LOGICAL, &
      & rank-x_div, 6, MPI_comm_world, reqB2, ierr)
      if (ierr/=0) stop 'Error with MPI_Irecv reqB2'
    end if

    ! Ensure all data has been received
    ! Left
    if (mod(rank+1, x_div) /= 1) then
      ! Check array 1 has sent
      call MPI_wait(sreqL1,sL1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqL1'

      ! Check array 2 has sent
      call MPI_wait(sreqL2,sL2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqL2'

      ! Ensure array 1 is received
      call MPI_wait(reqL1,rL1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqL1'
      ! Ensure array 2 is received
      call MPI_wait(reqL2,rL2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqL2'

      ! unblock region around test coordinates (when near edge)
      if (int(x) < 7) then
        occupied(int(x), int(y)) = .false.
        box(int(x), int(y), 1) = 0.0_dp
        box(int(x), int(y), 2) = 0.0_dp
      endif
    endif

    ! Right
    if (mod(rank+1, x_div) /= 0) then
      ! Check array 1 has sent
      call MPI_wait(sreqR1,sR1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqR1'
      ! Check array 2 has sent
      call MPI_wait(sreqR2,sR2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqR2'

      ! Ensure array 1 is received
      call MPI_wait(reqR1,rR1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqR1'
      ! Ensure array 2 is received
      call MPI_wait(reqR2,rR2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqR2'

      ! unblock region around test coordinates (when near edge)
      if (int(x) > x_max) then
        occupied(int(x), int(y)) = .false.
        box(int(x), int(y), 1) = 0.0_dp
        box(int(x), int(y), 2) = 0.0_dp
      endif
    endif

    ! Top
    if (rank+1 <= size - x_div) then
      ! Check array 1 has sent
      call MPI_wait(sreqT1,sT1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqT1'
      ! Check array 2 has sent
      call MPI_wait(sreqT2,sT2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqT2'

      ! Ensure array 1 is received
      call MPI_wait(reqT1,rT1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqT1'
      ! Ensure array 2 is received
      call MPI_wait(reqT2,rT2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqT2'

      ! Set arrays back to true values, while these coordinates are tested
      if (int(y) > y_max) then
        occupied(int(x), int(y)) = .false.
        box(int(x), int(y), 1) = 0.0_dp
        box(int(x), int(y), 2) = 0.0_dp
      endif
    endif

    ! Bottom
    if (rank+1 > x_div) then
      ! Check array 1 has sent
      call MPI_wait(sreqB1,sB1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqB1'
      ! Check array 2 has sent
      call MPI_wait(sreqB2,sB2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait sreqB2'

      ! Ensure array 1 is received
      call MPI_wait(reqB1,rB1status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqB1'
      ! Ensure array 2 is received
      call MPI_wait(reqB2,rB2status,ierr)
      if (ierr/=0) stop 'Error with MPI_wait reqB2'

      if (int(y) < 7) then
        occupied(int(x), int(y)) = .false.
        box(int(x), int(y), 1) = 0.0_dp
        box(int(x), int(y), 2) = 0.0_dp
      end if
    endif

    ! If there are circles nearby, see if the new circle will fit next to it
    outer: do j = -3, 3
      do k = -3, 3
        ! skip central cell
        if ( (j == 0) .and. (k == 0) ) cycle
        ! skip cells that are greater than 2*r away
        if ( (abs(j) == 3) .and. (abs(k) == 3) ) cycle
        ! if nearby gridpoint is occupied, test distance between circles
        if ( occupied(int(x)+j, int(y)+k) ) then
          !pythagoros' theroem
          D = (box(int(x)+j, int(y)+k, 1) - x)**2 &
          + (box(int(x)+j, int(y)+k, 2) - y)**2
          ! abandon coords if the distance between the two circles is too small
          if ( D < (2*r)**2 ) then
            pack = .false.
            exit outer
          endif
        end if
      end do
    end do outer

    ! If there is room for the circle, add it to the box
    if ( pack ) then
      ! print *, "Circle added", x, y, "on rank", rank
      N_circ = N_circ + 1
      occupied(int(x), int(y)) = .true.
      box(int(x), int(y), 1) = x
      box(int(x), int(y), 2) = y
      ! check it is working
      if ( rank == 0 .and. (int(x) < 50) .and. (int(y) < 50) ) then
        write(unit=11, fmt=*, iostat=ios) x, y
        if ( ios /= 0 ) stop "Write error in file unit 11"
      end if
      ! if ( i > 100000000) print*, i
    end if

    ! every 10,000 steps, write number of trial circle placements and value of P
      ! if ( modulo(i,10000) == 0 ) then
      !   write(unit=22, fmt=*, iostat=ios) i, pi*N_circ*r**2/L**2
      !   if ( ios /= 0 ) stop "Write error in file unit 22"
      ! endif
  end do

  ! write the final packing fraction and print to terminal
  ! write(unit=22, fmt=*, iostat=ios) i, pi*N_circ*r**2/L**2
  ! print*, pi*N_circ*r**2/L**2

  ! random packings have 55-64% packing fractions
  ! ideal packing of circles is 0.9069

  ! to plot on gnuplot:
  ! gnuplot> set style circle radius 1.234
  ! gnuplot> plot 'segment.dat' with circles

  ! Close data file
  if ( rank == 0 ) then
    close(unit=11, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file unit 11"
  endif

  if (allocated(box)) deallocate(box, stat=ierr)
  if ( ierr /= 0) print *, "box: Deallocation request denied"

  if (allocated(occupied)) deallocate(occupied, stat=ierr)
  if ( ierr /= 0) print *, "occupied: Deallocation request denied"

  call MPI_Finalize(ierr)
end program
