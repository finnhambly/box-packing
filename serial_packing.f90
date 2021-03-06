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
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300), L = 500
  real(kind=dp), parameter :: r = 1.234_dp, pi = 3.14159265358979_dp
  real(kind=dp), dimension(L,L,2) :: box ! 2500 point grid containing
                                        ! coordinates of packed circles
  logical, dimension(L,L) :: occupied
  real(kind=dp) :: x, y, D
  logical :: pack
  integer :: i, j, k, N_circ, left, right, above, below, ios
  type(rng_t):: rng   ! custom type for random number generator

  ! set occupation status of box as empty
  occupied = .false.
  box = 0.0_dp

  ! open file to write segment to
  open(unit=11, file='segment.dat', status='unknown', iostat=ios)
  if ( ios /= 0 ) stop "Error opening file segment.dat"

  open(unit=22, file='packingfrac.dat', status='unknown', iostat=ios)
  if ( ios /= 0 ) stop "Error opening file packingfrac.dat"

  ! random number seed
  call rng_seed(rng, 362436069)
  N_circ = 0
  do i = 1, 10000000
    ! generate new random coordinates within the box
    x = (L-2*r) * rng_uniform(rng) + r
    y = (L-2*r) * rng_uniform(rng) + r
    ! print*, "Testing coordinates", x, y

    ! set boolean test to default value
    pack = .true.

    ! if grid point already occupied, then skip these coordinates
    if (occupied(int(x),int(y))) cycle

    ! if not on the edge of the box... (for simplicity)
    if ( (int(x) /= 1) .and. (int(x) /= L) &
    .and. (int(y) /= 1) .and. (int(y) /= L) ) then
    ! ...look and, if nearest neighbours are occupied, skip
      if ( occupied(int(x)+1,int(y)) .or. occupied(int(x)-1,int(y)) .or. &
      occupied(int(x),int(y)+1) .or. occupied(int(x),int(y)-1)) cycle
    endif

    ! set the surrounding grid points for distance tests
    left = -3
    right = 3
    above = 3
    below = -3
    if (int(x) < 4) then
      left = 1 - int(x)
    end if
    if (int(y) < 4) then
      below = 1 - int(y)
    end if
    if (int(x) > L-3) then
      right = L - int(x)
    end if
    if (int(y) > L-3) then
      above = L - int(y)
    end if

    ! if there are circles nearby, see if the new circle will fit next to it
    outer: do j = left, right
      do k = below, above
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

    ! if there is room for the circle, add it to the box
    if ( pack ) then
      ! print *, "Circle added", x, y
      N_circ = N_circ + 1
      occupied(int(x), int(y)) = .true.
      box(int(x), int(y), 1) = x
      box(int(x), int(y), 2) = y
      if ( (int(x) < 14) .and. (int(y) < 14)) then
        write(unit=11, fmt=*, iostat=ios) x, y
        if ( ios /= 0 ) stop "Write error in file unit 11"
      end if
      ! if ( i > 100000000) print*, i
    end if

    ! every 10,000 steps, write number of trial circle placements and value of P
      if ( modulo(i,10000) == 0 ) then
        write(unit=22, fmt=*, iostat=ios) i, pi*N_circ*r**2/L**2
        if ( ios /= 0 ) stop "Write error in file unit 22"
      endif
  end do

  ! write the final packing fraction and print to terminal
  write(unit=22, fmt=*, iostat=ios) i, pi*N_circ*r**2/L**2
  print*, pi*N_circ*r**2/L**2

  ! random packings have 55-64% packing fractions
  ! ideal packing of circles is 0.9069

  close(unit=11, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 11"

  close(unit=22, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 22"

  ! to plot on gnuplot:
  ! gnuplot> set style circle radius 1.234
  ! gnuplot> plot 'segment.dat' with circles

end program
