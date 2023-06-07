MODULE utilities
   USE params
   USE constants

CONTAINS

   SUBROUTINE interp2d(map, nx, ny, x, y, value)

      ! Subroutine to interpolate a value from a map using bilinear
      ! interpolation

      IMPLICIT NONE

      INTEGER                                     :: nx, ny
      REAL*8                                      :: map(nx, ny)
      REAL*8 x, y, value

      REAL*8                                      :: t, u, t1, u1
      INTEGER                                     :: i, j

      i = int(x)
      j = int(y)
      t = x - i
      u = y - j
      t1 = 1 - t
      u1 = 1 - u

      ! The code assumes that the x,y values map to array indicies
      ! Depending on compiler options, this may lead to silent out-of-bounds access.
      ! We therefore add a check for this condition.
      ! Unfortunately I can't figure out how to fit this case into the testing framework.
      if (i == 0 .or. i > nx .or. j == 0 .or. j > ny) then
         write (*, *) 'Warning: Interpolation index out of bounds'
         stop
      end if

      value = t1*u1*map(i, j) + t*u1*map(i + 1, j) + t*u*map(i + 1, j + 1) + t1*u*map(i, j + 1)

   END SUBROUTINE interp2d

!======================================================================

   SUBROUTINE interp1d(f, x, n, x0, f0)

! Subroutine to interpolate a value f0 from a FUNCTION f using linear
! interpolation

      IMPLICIT NONE

      INTEGER                                     :: n
      REAL*8                                      :: f(n), x(n)
      REAL*8                                      :: x0, f0
      REAL*8                                      :: t, u, t1, u1
      INTEGER                                     :: i, j

      ! First find two x values between which x0 lies:

      IF (x(1) >= x0) THEN
         f0 = f(1)
         RETURN
      ELSEIF (x(n) <= x0) THEN
         f0 = f(n)
         RETURN
      ELSE
         CALL locate(x, n, x0, i)
      END IF

      ! Now do the interpolation:

      f0 = f(i) + (f(i + 1) - f(i))*(x0 - x(i))/(x(i + 1) - x(i))

      if (isnan(f0)) then
         if (f(i) .eq. f(i + 1)) then
            f0 = f(i)
         end if
      end if

   END SUBROUTINE interp1d

!-----------------------------------------------------------------------

   subroutine interp1d_even(f, x, n, x0, f0)

      ! Subroutine to interpolate a value f0 from a function f using linear
      ! interpolation, on an evenly spaced grid

      implicit none

      integer n
      double precision f(n), x(n)
      double precision x0, f0
      double precision xmin, xmax, dx, wx
      integer i1

      xmin = x(1)
      xmax = x(n)
      dx = x(2) - x(1)

      ! write(*, *) 'interp1d_even called with x0', x0

      if (xmin > xmax) then
         write (*, *) 'Interpolation bounds reversed'
      end if

      if (x0 < xMin .and. x0 / xmin < 1d-6) then
         write (*, *) 'Interpolation value slightly low. Probably a floating point error, adjusting to limit'
         x0 = xmin
      else if (x0 > xMax .and. x0 / xmax < 1d-6) then
         write (*, *) 'Interpolation value slightly high. Probably a floating point error, adjusting to limit'
         x0 = xmin
      end if

      if (x0 .ge. xmin .and. x0 .le. xmax) then
         i1 = int((x0 - xmin)/dx) + 1
         wx = (x0 - x(i1))/dx
         f0 = (1 - wx)*f(i1) + wx*f(i1 + 1)
      else
         write (*, *) 'x = ', x0, '  is out of table range xmin: ', xmin, ' xmax: ', xmax
         write (*, *) '10**x:', 10**x0, 10**xmin, 10**xmax
         stop
      end if
      return

   end subroutine interp1d_even

!-----------------------------------------------------------------------

   SUBROUTINE locate(xx, n, x, j)

      IMPLICIT NONE

      INTEGER                                     ::n, j
      REAL*8                                      ::xx(n), x
      INTEGER                                     ::jl, ju, jm

      jl = 0
      ju = n + 1

11    IF ((ju - jl) .GT. 1) THEN
         jm = (ju + jl)/2
         IF ((xx(n) .GT. xx(1)) .eqv. (x .GT. xx(jm))) THEN
            jl = jm
         ELSE
            ju = jm
         END IF
         GO TO 11
      END IF

      j = jl

      RETURN
   END SUBROUTINE locate

!=======================================================================

   FUNCTION polar(r, theta, Amp)

      REAL*8                                      ::var_i, var_j, c1
      PARAMETER(var_i=20., var_j=20., c1=1.)
      REAL*8                                      ::r, theta, Amp
      REAL*8                                      ::polar

      polar = c1*Amp*r*cos(theta)*exp(-(r**2*(cos(theta))**2 &
         /(2*var_i)) - (r**2*(sin(theta))**2/(2*var_j)))

      RETURN
   END FUNCTION polar

!=======================================================================

   FUNCTION phlog10(x)

      IMPLICIT NONE

      REAL*8                                      ::x, phlog10

      IF (x .LT. 0.0) THEN
         if (myID == 0) WRITE (*, *) 'ERROR: phlog10: error, tried to take log of ', x
         STOP
      ELSEIF (x .LT. 1d-45) THEN
         phlog10 = -45.0
      ELSE
         phlog10 = log10(x)
      END IF

      RETURN
   END FUNCTION phlog10

!=======================================================================
! Numerical recipes integration routine-the supplied real function
! func is integrated between limits A and B, and the result returned as
! the integral S. Progressively finer discretisations are used until
! either
!
!     | S_{j}-S_{j-1} |
!     ------------------- < eps
!          | S_{j-1} |
!
! or the range has to be split into more than 2^jmax sections. For
! smooth functions the latter shouldnt happen...

   SUBROUTINE qtrap(func, lim1, lim2, eps, S)

      IMPLICIT NONE

      REAL*8                                      ::S, lim1, lim2
      REAL*8                                      ::a1, b1, eps, olds, sign, t1, t2
      INTEGER                                     :: j, jmax
      PARAMETER(jmax=25)

      INTERFACE
         FUNCTION func(zz)
            REAL*8 zz, func
         END FUNCTION func
      END INTERFACE

!-----------------------------------

      a1 = lim1
      b1 = lim2

      ! First catch some stupidities:

      IF (a1 == b1) THEN
         S = 0.d0
         GOTO 30
      ELSEIF (a1 > b1) THEN
         olds = b1
         b1 = a1
         a1 = olds
         sign = -1.d0
      ELSE
         sign = 1.d0
      END IF

      olds = -1.d30
      s = 0.d0

      DO j = 1, jmax
         CALL trapzd(func, a1, b1, S, j)
         IF (j == 1) THEN
            t1 = func(a1)
            t2 = func(b1)
         END IF

         IF (abs(s - olds) < eps*abs(olds) .OR. abs(s - olds) == 0.d0) GOTO 20
         IF (j == jmax) GOTO 10

         olds = s
      END DO

10    if (myID == 0) then
         WRITE (*, *) 'QTRAP error: too many steps...'
         WRITE (*, *) '   S = ', S
         WRITE (*, *) '   oldS = ', oldS
         WRITE (*, *) '   % difference = ', 100.0*(S - oldS)/oldS
         WRITE (*, *) '   limits were ', a1, b1
         WRITE (*, *) '   FUNCTION values at limits were ', t1, t2
         t1 = func(0.5*(a1 + b1))
         WRITE (*, *) '   FUNCTION value at midpoint was ', t1
      end if
      GOTO 30

20    S = S*sign
30    RETURN
   END SUBROUTINE qtrap

!-----------------------------------------------------------------------

   SUBROUTINE trapzd(func, a1, b1, S, n)

      IMPLICIT NONE

      REAL*8                                     :: a1, b1, S
      INTEGER                                     :: n
      INTEGER                                      ::it, j
      REAL*8                                     :: tnm, del, x, sum

      INTERFACE
         FUNCTION func(zz)
            REAL*8                                     :: zz, func
         END FUNCTION func
      END INTERFACE

      IF (n .EQ. 1) THEN
         S = 0.5*(b1 - a1)*(func(a1) + func(b1))
         it = 1
      ELSE
         it = 2**(n - 2)
         tnm = it*1.0
         del = (b1 - a1)/tnm
         x = a1 + 0.5*del
         sum = 0.
         DO 11 j = 1, it
            sum = sum + func(x)
            x = x + del
11       CONTINUE
         s = 0.5*(s + (b1 - a1)*sum/tnm)
      END IF

      RETURN
   END SUBROUTINE trapzd

!=======================================================================

   FUNCTION integrate(y, x, n, xmin, xmax)
!
!  Numerical integration of function y(x): y and x are arrays of
!  length n, and y is integrated between xmin and ymax using the
!  trapezium rule.
!
      IMPLICIT NONE
!
      INTEGER                                      ::n, i, j, k
      REAL*8                                     :: y(n), x(n)
      REAL*8                                      ::sum, integrate
      REAL*8                                     :: xmin, xmax, ymin, ymax
!
!-----------------------------------------------------------------------
!
      sum = 0.0
!
      DO i = 1, n
         IF (xmin .LT. x(i)) GOTO 10
      END DO
!
10    j = i
      ymin = y(j) + (y(j - 1) - y(j))*(x(j) - xmin)/(x(j) - x(j - 1))
      sum = (ymin + y(j))*(x(j) - xmin)/2.0
!
      DO i = j, n
         IF (xmax .LE. x(i)) GOTO 20
      END DO
!
20    k = i - 1
      ymax = y(k) - (y(k) - y(k + 1))*(xmax - x(k))/(x(k + 1) - x(k))
      sum = sum + (ymax + y(k))*(xmax - x(k))/2.0
!
      DO i = j, k - 1
         sum = sum + (y(i + 1) + y(i))*(x(i + 1) - x(i))/2.0
      END DO
!
      integrate = sum
!
      RETURN
   END FUNCTION integrate
!

!======================================================================

   FUNCTION nintegrate(y, x, n, xmin, xmax)
!
!  Numerical integration of function y(x): y and x are arrays of
!  length n, and y is integrated between xmin and ymax using the
!  trapezium rule.
!
      IMPLICIT NONE
!
      INTEGER                                     :: n, i, j, k
      REAL*8                                      :: y(n), x(n)
      REAL*8                                      :: sum, nintegrate
      REAL*8                                      :: xmin, xmax, ymin, ymax
!
!-----------------------------------------------------------------------
!
      sum = 0.0
!
      DO i = 1, n
         IF (xmin .LT. x(i)) GOTO 10
      END DO
!
10    j = i
      ymin = y(j) + (y(j - 1) - y(j))*(x(j) - xmin)/(x(j) - x(j - 1))
      sum = (ymin + y(j))*(x(j) - xmin)/2.0
!
      DO i = j, n
         IF (xmax .le. x(i)) GOTO 20
      END DO
!
20    k = i - 1
      ymax = y(k) - (y(k) - y(k + 1))*(xmax - x(k))/(x(k + 1) - x(k))
      sum = sum + (ymax + y(k))*(xmax - x(k))/2.0
!
      DO i = j, k - 1
         sum = sum + (y(i + 1) + y(i))*(x(i + 1) - x(i))/2.0
      END DO
!
      nintegrate = sum
!
      RETURN
   END FUNCTION nintegrate
!
!======================================================================

   FUNCTION atanh(x)

      REAL*8                                      :: x
      REAL*8                                      :: atanh

      IF (x*x .ge. 1.d0) THEN
         if (myID == 0) WRITE (*, *) 'ERROR: This should not happen!, x=', x
         atanh = 0.0
      ELSE
         atanh = 0.5d0*dlog((1.d0 + x)/(1.d0 - x))
      END IF

      RETURN
   END FUNCTION atanh

!=======================================================================

   FUNCTION FINDR(FUNC, X1, X2, XACC)

      IMPLICIT NONE

      REAL*8                                      :: FINDR, X1, X2, XACC
      REAL*8                                      :: FMID, F, DX, XMID
      INTEGER                                     :: J
      INTEGER, PARAMETER                          :: JMAX = 40

      INTERFACE
         FUNCTION func(zz)
            REAL*8                                      :: zz, func
         END FUNCTION func
      END INTERFACE

      FMID = FUNC(X2)
      F = FUNC(X1)
      IF ((F*FMID) .GT. 0.) THEN
         FINDR = 0.0
         RETURN
      END IF

      IF (F .LT. 0.) THEN
         FINDR = X1
         DX = X2 - X1
      ELSE
         FINDR = X2
         DX = X1 - X2
      END IF
      DO 11 J = 1, JMAX
         DX = DX*.5
         XMID = FINDR + DX
         FMID = FUNC(XMID)
         IF (FMID .LE. 0.) FINDR = XMID
         IF (ABS(DX) .LT. XACC .OR. FMID .EQ. 0.) RETURN
11    CONTINUE
      if (myID == 0) write (*, *) 'ERROR: too many bisections'
      stop
   END FUNCTION FINDR

!=======================================================================
   !Real Incomplete Beta Function
   !Numerical Recipes

   REAL*8 FUNCTION betai(a, b, x)

      IMPLICIT NONE
      REAL*8                                       ::a, b, x
      REAL*8                                       ::bt

      IF (x < 0. .OR. x > 1.) then
         if (myID == 0) write (*, *) 'ERROR: bad argument x in betai'
         stop
      end if

      IF (x == 0. .OR. x == 1.) THEN
         bt = 0.0
      ELSE
         bt = (x**a)*((1.-x)**b)
      END IF
      IF (x < (a + 1.)/(a + b + 2.)) THEN
         betai = bt*betacf(a, b, x)/a
      ELSE
         betai = 1.-bt*betacf(b, a, 1.-x)/b
      END IF
   END FUNCTION betai

!=======================================================================
   !Used by betai: Evaluates continued fraction for incomplete beta function by modified
   !Lentz's method
   !Numerical Recipes

   REAL*8 FUNCTION betacf(a, b, x)

      INTEGER                                      :: MAXIT
      REAL*8                                       :: a, b, x, EPS, FPMIN
      PARAMETER(MAXIT=100, EPS=3.d-7, FPMIN=1.d-30)

      INTEGER                                     :: m, m2
      REAL*8                                      :: aa, c, d, del, h, qab, qam, qap

      qab = a + b         !These q's will be used in factors that occur in the coefficients
      qap = a + 1.
      qam = a - 1.
      c = 1.                 !First step of Lentz\u2019s method.
      d = 1.-qab*x/qap
      IF (abs(d) < FPMIN) d = FPMIN
      d = 1./d
      h = d
      DO m = 1, MAXIT
         m2 = 2*m
         aa = m*(b - m)*x/((qam + m2)*(a + m2))
         d = 1.+aa*d !One step (the even one) of the recurrence.
         IF (abs(d) < FPMIN) d = FPMIN
         c = 1.+aa/c
         IF (abs(c) < FPMIN) c = FPMIN
         d = 1./d
         h = h*d*c
         aa = -(a + m)*(qab + m)*x/((a + m2)*(qap + m2))
         d = 1.+aa*d !Next step of the recurrence (the odd one).
         IF (abs(d) < FPMIN) d = FPMIN
         c = 1.+aa/c
         IF (abs(c) < FPMIN) c = FPMIN
         d = 1./d
         del = d*c
         h = h*del
         IF (abs(del - 1.) < EPS) THEN
            betacf = h
            RETURN
         END IF
      END DO
      if (myID == 0) write (*, *) 'ERROR: a or b too big, or MAXIT too small in betacf'
      betacf = h
      RETURN
   END FUNCTION betacf

!=======================================================================
   !root funding using the Brent method
   !Numerical Recipes

   FUNCTION zbrent(funct, x1, x2, tol)
      INTEGER                                       :: ITMAX
      REAL*8                                        :: zbrent, tol, x1, x2, func, EPS
      PARAMETER(ITMAX=100, EPS=3.d-8)
      INTEGER iter
      REAL*8                                        ::a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm

      INTERFACE
         FUNCTION funct(zz)
            REAL*8 zz, funct
         END FUNCTION funct
      END INTERFACE

      a = x1
      b = x2
      fa = funct(a)
      fb = funct(b)
      IF ((fa > 0. .AND. fb > 0.) .OR. (fa < 0. .AND. fb < 0.)) then
         if (myID == 0) write (*, *) 'root must be bracketed for zbrent'
         stop
      end if
      c = b
      fc = fb
      DO 11 iter = 1, ITMAX
         IF ((fb > 0. .AND. fc > 0.) .OR. (fb < 0. .AND. fc < 0.)) THEN
            c = a
            fc = fa
            d = b - a
            e = d
         END IF
         IF (abs(fc) < abs(fb)) THEN
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
         END IF
         tol1 = 2.*EPS*abs(b) + 0.5*tol
         xm = .5*(c - b)
         IF (abs(xm) <= tol1 .OR. fb == 0.) THEN
            zbrent = b
            RETURN
         END IF
         IF (abs(e) >= tol1 .AND. abs(fa) > abs(fb)) THEN
            s = fb/fa
            IF (a == c) THEN
               p = 2.*xm*s
               q = 1.-s
            ELSE
               q = fa/fc
               r = fb/fc
               p = s*(2.*xm*q*(q - r) - (b - a)*(r - 1.))
               q = (q - 1.)*(r - 1.)*(s - 1.)
            END IF
            IF (p > 0.) q = -q
            p = abs(p)
            IF (2.*p < min(3.*xm*q - abs(tol1*q), abs(e*q))) THEN
               e = d
               d = p/q
            ELSE
               d = xm
               e = d
            END IF
         ELSE
            d = xm
            e = d
         END IF
         a = b
         fa = fb
         IF (abs(d) > tol1) THEN
            b = b + d
         ELSE
            b = b + sign(tol1, xm)
         END IF
         fb = funct(b)
11    CONTINUE
      if (myID == 0) write (*, *) 'zbrent exceeding maximum iterations'
      zbrent = b
      RETURN
   END FUNCTION zbrent

!=======================================================================
   !binary search
   FUNCTION binSearch(n, list, value)

      IMPLICIT NONE

      REAL*8                                        ::  binSearch !largest i for which which has list(i)>=value
      INTEGER                                       ::  n !num of elements in the list
      REAL*8                                        ::  list(n) !arrray that needs to be searched
      REAL*8                                        ::  value !value ot be searched
      INTEGER                                       :: i, j, m1, l1, r1

      !sanity check
      IF (list(1) > value) THEN
         binSearch = 1
         RETURN
      ELSEIF (list(n) < value) THEN
         binSearch = n + 1
         RETURN
      END IF

      m1 = n/2
      l1 = 1
      r1 = n
      DO
         IF (value == list(m1)) THEN
            binSearch = m1
            RETURN
         ELSEIF ((l1 - r1)**2 == 1) THEN
            binSearch = r1
            RETURN
         ELSEIF (value < list(m1)) THEN
            r1 = m1
         ELSEIF (value > list(m1)) THEN
            l1 = m1
         END IF
         m1 = (l1 + r1)/2
      END DO

   END FUNCTION binSearch

!=======================================================================
!        lookup a 2D table, output 2 values
   SUBROUTINE lookUp2D(num1, num2, val1, val2, tab1, tab2, tab3, v)

      IMPLICIT NONE

      INTEGER                                      ::  num1 !total no. of entries in dim 1
      INTEGER                                      ::  num2 !total no. of entries in dim 2
      REAL*8                                       ::  val1 !1st value to look for
      REAL*8                                       ::  val2 !2nd value to look for
      REAL*8                                       :: v !output
      REAL*8                                       ::  tab1(num1) !1st dim
      REAL*8                                       ::  tab2(num2) !2nd dim
      REAL*8                                       ::  tab3(num1, num2) !lookup table

      INTEGER                                      ::  i, j, k
      REAL*8                                       ::  x, y

      !Lookup 1st dim
      i = binSearch(num1, tab1, val1)

      !Lookup 2nd dim
      j = binSearch(num2, tab2, val2)

      IF (tab1(i) == val1 .AND. tab2(j) == val2) THEN
         v = tab3(i, j)
      ELSEIF (tab1(i) == val1 .AND. .not. (j == 0 .OR. j == num2)) THEN
         v = tab3(i, j - 1) + (tab3(i, j) - tab3(i, j - 1))*(val2 - tab2(j - 1))/(tab2(j) - tab2(j - 1))
      ELSEIF (tab2(j) == val2 .AND. .not. (i == 0 .OR. i == num1)) THEN
         v = tab3(i - 1, j) + (tab3(i, j) - tab3(i - 1, j))*(val1 - tab1(i - 1))/(tab1(i) - tab1(i - 1))
         !sanity check
      ELSEIF (i == 0 .OR. j == 0 .OR. i > num1 .OR. j > num2) THEN
         if (myID == 0) WRITE (*, *) "ERROR: problem in lookup2D, value to look for is not within the table bounds"
         stop
      ELSE
         !bilinear interpolation
         x = (i - 1.) + (val1 - tab1(i - 1))/(tab1(i) - tab1(i - 1))
         y = (j - 1.) + (val2 - tab2(j - 1))/(tab2(j) - tab2(j - 1))
         CALL interp2d(tab3, num1, num2, x, y, v)
      END IF

   END SUBROUTINE lookUp2D

!=======================================================================

   SUBROUTINE lookUp1D(num1, val1, tab1, tab3, v)
      ! lookup a 2D table, output 2 values

      IMPLICIT NONE

      INTEGER                                     :: num1 !total no. of entries in dim 1
      REAL*8                                      :: val1 !1st value to look for
      REAL*8                                      :: v !output
      REAL*8                                      ::tab1(num1) !1st dim
      REAL*8                                      ::tab3(num1) !lookup table
      INTEGER                                     ::i, j, k
      REAL*8                                      ::x, y

      !Lookup 1st dim
      i = binSearch(num1, tab1, val1)
      IF (tab1(i) == val1) THEN
         v = tab3(i)
         !sanity check
      ELSEIF (i == 0 .OR. i > num1) THEN
         if (myID == 0) WRITE (*, *) "problem in lookup1D, value to look for is not within the table bounds"
         stop
      ELSE
         !bilinear interpolation
         v = tab3(i - 1) + (tab3(i) - tab3(i - 1))*(val1 - tab1(i - 1))/(tab1(i) - tab1(i - 1))
      END IF

   END SUBROUTINE lookUp1D

!=======================================================================
!        check if the given point p is inside the triangle with vertices a, b & c
   LOGICAL FUNCTION inTriangle(a, b, c, p)

      IMPLICIT NONE

      ! input variables
      REAL*8                                      :: a(2), b(2), c(2) ! the 3 vertices of the triangle
      REAL*8                                      :: p(2) ! the given point
      ! work variables
      REAL*8                                      :: v0(2), v1(2), v2(2), dot00, dot01, dot02, dot11, dot12, invDenom, u, v

      v0 = c - a
      v1 = b - a
      v2 = p - a

      ! compute dot products
      dot00 = dot_product(v0, v0)
      dot01 = dot_product(v0, v1)
      dot02 = dot_product(v0, v2)
      dot11 = dot_product(v1, v1)
      dot12 = dot_product(v1, v2)

      ! Compute barycentric coordinates
      invDenom = 1d0/(dot00*dot11 - dot01*dot01)
      u = (dot11*dot02 - dot01*dot12)*invDenom
      v = (dot00*dot12 - dot01*dot02)*invDenom

      ! Check if point is in triangle
      IF ((u > 0d0) .AND. (v > 0d0) .AND. (u + v < 1d0)) THEN
         inTriangle = .TRUE.
      ELSE
         inTriangle = .FALSE.
      END IF

   END FUNCTION inTriangle

!=======================================================================
   FUNCTION Gammafun(x)

      IMPLICIT NONE
      REAL*8                                      ::   Gammafun, x
      INTEGER                                       ::   i, j, k, nn
      REAL*8                                      ::  w, y, fun
      REAL*8                                       ::  p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13
      PARAMETER(p0=0.999999999999999990d+00, p1=-0.422784335098466784d+00)
      PARAMETER(p2=-0.233093736421782878d+00, p3=0.191091101387638410d+00)
      PARAMETER(p4=-0.024552490005641278d+00, p5=-0.017645244547851414d+00)
      PARAMETER(p6=0.008023273027855346d+00, p7=-0.000804329819255744d+00)
      PARAMETER(p8=-0.000360837876648255d+00, p9=0.000145596568617526d+00)
      PARAMETER(p10=-0.000017545539395205d+00, p11=-0.000002591225267689d+00)
      PARAMETER(p12=0.000001337767384067d+00, p13=-0.000000199542863674d+00)

      nn = nint(x - 2)
      w = x - (nn + 2)
      y = ((((((((((((p13*w + p12)*w + p11)*w + p10)* &
         w + p9)*w + p8)*w + p7)*w + p6)*w + p5)* &
         w + p4)*w + p3)*w + p2)*w + p1)*w + p0
      IF (nn .GT. 0) THEN
         w = x - 1
         DO k = 2, nn
            w = w*(x - k)
         END DO
      ELSE
         w = 1
         DO k = 0, -nn - 1
            y = y*(x + k)
         END DO
      END IF
      fun = w/y
      Gammafun = fun
   END FUNCTION Gammafun
!======================================
   FUNCTION INCOG(A, X)

      IMPLICIT NONE
      REAL*8           ::INCOG
      REAL*8           :: A, X
      REAL*8           :: GIN, GIM, GIP
      INTEGER          :: K
      REAL*8           :: XAM, GA, S, R, T0

      XAM = -X + A*DLOG(X)
      IF (XAM .GT. 700.0 .OR. A .GT. 170.0) THEN
         WRITE (*, *) 'a and/or x too large'
         STOP
      END IF
      IF (X .EQ. 0.0) THEN
         GIN = 0.0
         CALL GAMMA(A, GA)
         GIM = GA
         GIP = 0.0
      ELSE IF (X .LE. 1.0 + A) THEN
         S = 1.0D0/A
         R = S
         DO K = 1, 60
            R = R*X/(A + K)
            S = S + R
            IF (DABS(R/S) .LT. 1.0D-15) GO TO 15
         END DO
15       GIN = DEXP(XAM)*S
         CALL GAMMA(A, GA)
         GIP = GIN/GA
         GIM = GA - GIN
      ELSE IF (X .GT. 1.0 + A) THEN
         T0 = 0.0D0
         DO K = 60, 1, -1
            T0 = (K - A)/(1.0D0 + K/(X + T0))
         END DO
         GIM = DEXP(XAM)/(X + T0)
         CALL GAMMA(A, GA)
         GIN = GA - GIM
         GIP = 1.0D0 - GIM/GA
      END IF
      INCOG = GIN
   END FUNCTION INCOG

!      ===================================================
   SUBROUTINE GAMMA(X, GA)

!     ==================================================
!     Purpose: Compute gamma function â(x)
!     Input :  x  --- Argument of â(x)
!                     ( x is not equal to 0,-1,-2,úúú)
!     Output:  GA --- â(x)
!     ==================================================
      IMPLICIT NONE

      INTEGER          :: M1, K, M
      REAL*8           :: PI, XAM, X, A, GA, Z, R, GR

      REAL*8    :: G(26)
      PI = 3.141592653589793D0
      IF (X .EQ. INT(X)) THEN
         IF (X .GT. 0.0D0) THEN
            GA = 1.0D0
            M1 = X - 1
            DO K = 2, M1
               GA = GA*K
            END DO
         ELSE
            GA = 1.0D+300
         END IF
      ELSE
         IF (DABS(X) .GT. 1.0D0) THEN
            Z = DABS(X)
            M = INT(Z)
            R = 1.0D0
            DO K = 1, M
               R = R*(Z - K)
            END DO
            Z = Z - M
         ELSE
            Z = X
         END IF
         DATA G/1.0D0, 0.5772156649015329D0, &
            -0.6558780715202538D0, -0.420026350340952D-1, &
            0.1665386113822915D0, -.421977345555443D-1, &
            -.96219715278770D-2, .72189432466630D-2, &
            -.11651675918591D-2, -.2152416741149D-3, &
            .1280502823882D-3, -.201348547807D-4, &
            -.12504934821D-5, .11330272320D-5, &
            -.2056338417D-6, .61160950D-8, &
            .50020075D-8, -.11812746D-8, &
            .1043427D-9, .77823D-11, &
            -.36968D-11, .51D-12, &
            -.206D-13, -.54D-14, .14D-14, .1D-15/

         GR = G(26)
         DO K = 25, 1, -1
            GR = GR*Z + G(K)
         END DO
         GA = 1.0D0/(GR*Z)
         IF (DABS(X) .GT. 1.0D0) THEN
            GA = GA*R
            IF (X .LT. 0.0D0) GA = -PI/(X*GA*DSIN(PI*X))
         END IF
      END IF
      RETURN
   END SUBROUTINE GAMMA
!======================================
END MODULE utilities
