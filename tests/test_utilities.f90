module test_utilities
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use utilities
   implicit none
   private

   public :: collect_utilities

contains

   ! Collect all exported unit tests
   subroutine collect_utilities(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("interp2d", test_interp2d), &
         new_unittest("interp1d_even", test_interp2d), &
         new_unittest("qtrap", test_qtrap) &
      !  new_unittest("invalid", test_invalid, should_fail=.true.) &
         ]

   end subroutine collect_utilities

   subroutine test_interp2d(error)
      implicit none
      type(error_type), allocatable, intent(out) :: error

      integer :: nx, ny
      real*8 :: map(5, 5)
      real*8 :: x, y, value, difference

      ! Initialize map values
      nx = 5
      ny = 5
      map = RESHAPE([0.0, 1.0, 2.0, 3.0, 4.0, &
         5.0, 6.0, 7.0, 8.0, 9.0, &
         10.0, 11.0, 12.0, 13.0, 14.0, &
         15.0, 16.0, 17.0, 18.0, 19.0, &
         20.0, 21.0, 22.0, 23.0, 24.0], [nx, ny])

      ! Test interpolation at a set of points
      x = 2.5
      y = 3.5
      call interp2d(map, nx, ny, x, y, value)
      !   difference = (value - 14.d0) / 14.d0
      !   call check(error, value < 0.001)
      call check(error, value, 14.d0, thr=1.d-5, rel=.true.)
      if (allocated(error)) return

      x = 3.8
      y = 1.2
      call interp2d(map, nx, ny, x, y, value)
      call check(error, value, 3.8d0, thr=1.d-3, rel=.true.)
      if (allocated(error)) return

   end subroutine test_interp2d

   subroutine test_interp1d_even(error)
      implicit none
      type(error_type), allocatable, intent(out) :: error

      integer :: n=5
      real*8 :: f(5), x(5), x0, value

      f = [0.0, 3.0, 6.0, 9.0, 12.0]
      x = [0.0, 1.0, 2.0, 3.0, 4.0]

      call interp1d_even(f, x, n, 0.5d0, value)
      call check(error, value, 1.5d0, thr=1.d-6, rel=.true.)
      if (allocated(error)) return

      call interp1d_even(f, x, n, 0.0d0, value)
      call check(error, value, 0d0, thr=1.d-6, rel=.true.)
      if (allocated(error)) return

      call interp1d_even(f, x, n, 4.0d0, value)
      call check(error, value, 12d0, thr=1.d-6, rel=.true.)
      if (allocated(error)) return

      call interp1d_even(f, x, n, 3.0d0, value)
      call check(error, value, 9d0, thr=1.d-6, rel=.true.)
      if (allocated(error)) return

      call interp1d_even(f, x, n, 1.75d0, value)
      call check(error, value, 5.25d0, thr=1.d-6, rel=.true.)
      if (allocated(error)) return

   end subroutine test_interp1d_even

   function integrand1(x)
      real*8 :: x
      real*8 :: integrand1

      integrand1 = 3 * x**2 + x

      return
   end function integrand1

   function integrand2(x)
      real*8 :: x
      real*8 :: integrand2

      integrand2 = 0.d0

      return
   end function integrand2

   function integrand3(x)
      real*8 :: x
      real*8 :: integrand3

      integrand3 = 1/x

      return
   end function integrand3

   function integrand4(x)
      real*8 :: x
      real*8 :: integrand4

      integrand4 = cos(x)

      return
   end function integrand4

   function integrand5(x)
      real*8 :: x
      real*8 :: integrand5

      integrand5 = exp(x)

      return
   end function integrand5

   subroutine test_qtrap(error)
      implicit none
      type(error_type), allocatable, intent(out) :: error

      real*8 :: sum_holder

      sum_holder = 0

      call qtrap(integrand1, 0d0, 10d0, 1.d-12, sum_holder)
      call check(error, sum_holder, 1050.d0, thr=1d-6)
      if (allocated(error)) return

      call qtrap(integrand2, -13.5d0, 10d0, 1.d-12, sum_holder)
      call check(error, sum_holder, 0d0, thr=1d-6)
      if (allocated(error)) return

      call qtrap(integrand3, 1d0, 10d0, 1.d-12, sum_holder)
      call check(error, sum_holder, log(10d0) - log(1d0), thr=1d-6)
      if (allocated(error)) return

      call qtrap(integrand4, 0d0, Pi/2d0, 1.d-6, sum_holder)
      call check(error, sum_holder, 1d0, thr=1d-6)
      if (allocated(error)) return

      call qtrap(integrand5, 0d0, 10d0, 1.d-6, sum_holder)
      call check(error, sum_holder, 1d0 - exp(-10d0), thr=1d-6)
      if (allocated(error)) return

   end subroutine test_qtrap

end module test_utilities
