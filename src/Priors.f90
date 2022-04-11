MODULE  priors1

CONTAINS

!=======================================================================
! Prior distribution functions: r is a unIForm deviate from the unit
! Cube, Prior returns the corresponding deviate drawn from the desired
! distribution (labelled k).
!
! 
!
       FUNCTION Prior(k,r,x1,x2,x3)

            IMPLICIT NONE
	    
	        INTEGER                         :: k
	        REAL*8                           ::   r,x1,x2,x3,Prior
   

 	           IF (k==0) THEN
 	                  Prior=DeltaFunctionPrior(r,x1,x2)
  	          ELSEIF (k==1) THEN
 	                  Prior=UniformPrior(r,x1,x2)
  	          ELSEIF (k==2) THEN
 	                  Prior=LogPrior(r,x1,x2)
  	          ELSEIF (k==3) THEN
 	                  Prior=GaussianPrior(r,x1,x2)
                  ELSEIF (k==4) THEN
  	                  Prior=LogNormalPrior(r,x1,x2)
                  ELSEIF (k==5) THEN
 	                  Prior=SinPrior(r,x1,x2)
                  ELSEIF (k==6) THEN
 	                  Prior=CauchyPrior(r,x1,x2)
                  ELSEIF (k==7) THEN
  	                  WRITE(*,*)"can't use the spectral index prior for any other parameter"
 	                  STOP
                  ELSEIF (k==8) THEN
 	                  WRITE(*,*)"can't use the mass FUNCTION prior for any parameter other than M200 & z"
 	                  STOP
                  ELSEIF (k==10) THEN
  	                 Prior=ExpPrior(r,x1,x2,x3)
                  ELSEIF (k==11) THEN
 	                  Prior=PowerPrior(r,x1,x2,x3)
                  ELSEIF (k==12) THEN
 	                  Prior=TruncatedGaussianPrior(r,x1,x2,x3)
                  ENDIF

            RETURN
       END FUNCTION Prior

!=======================================================================
! Uniform[0:1]  ->  Delta[x1]

       FUNCTION DeltaFunctionPrior(r,x1,x2)

            IMPLICIT NONE

                REAL*8                          ::  r,x1,x2,DeltaFunctionPrior

                   DeltaFunctionPrior=x1

            RETURN
       END FUNCTION DeltaFunctionPrior

!=======================================================================
! Uniform[0:1]  ->  Uniform[x1:x2]

       FUNCTION UniformPrior(r,x1,x2)

            IMPLICIT NONE

                REAL*8                         ::  r,x1,x2,UniformPrior

                   UniformPrior=x1+r*(x2-x1)

            RETURN
      END FUNCTION UniformPrior

!=======================================================================
! Uniform[0:1]  ->  LogUniform[x1:x2]

       FUNCTION LogPrior(r,x1,x2)

            IMPLICIT NONE

                REAL*8                         ::  r,x1,x2,LogPrior
                REAL*8                         ::  lx1,lx2

                   IF (r < 0.0d0) THEN
                       LogPrior=-1.0d32
                   ELSE
                       lx1=dlog10(x1)
                       lx2=dlog10(x2)
                       LogPrior=10.d0**(lx1+r*(lx2-lx1))
                   ENDIF

            RETURN
       END FUNCTION LogPrior

!=======================================================================
! Uniform[0:1]  ->  Sin[x1:x2]  (angles in degrees):

       FUNCTION SinPrior(r,x1,x2)

            IMPLICIT NONE

                REAL*8                                     :: r,x1,x2,SinPrior
                REAL*8                                     ::  cx1,cx2
                REAL*8 ,PARAMETER                          :: deg2rad=0.017453292
      

                   cx1=cos(x1*deg2rad)
                   cx2=cos(x2*deg2rad)
                   SinPrior=1.d0*acos(cx1+r*(cx2-cx1))

            RETURN
       END FUNCTION SinPrior

!=======================================================================
! Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]

       FUNCTION GaussianPrior(r,mu,sigma)

            IMPLICIT NONE

                REAL*8                                     ::  r,mu,sigma,GaussianPrior
                REAL*8 ,PARAMETER                          ::  SqrtTwo=1.414213562d0
      

                   IF (r < 1.0d-16.or.(1.0d0-r) < 1.0d-16) THEN
                       GaussianPrior=-1.0d32
                   ELSE
                       GaussianPrior=mu+sigma*SqrtTwo*dierfc(2.d0*(1.d0-r))
                   ENDIF

            RETURN
       END FUNCTION GaussianPrior

!=======================================================================
! x = Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2] with (x > limit IF limit < mu) 
! or (x < limit IF limit > mu)

       FUNCTION TruncatedGaussianPrior(r,mu,sigma,limit)

            IMPLICIT NONE

                REAL*8                         ::  r,mu,sigma,limit,TruncatedGaussianPrior
                REAL*8, PARAMETER              ::  SqrtTwo=1.414213562d0
      

                   IF (r < 1.0d-16.or.(1.0d0-r) < 1.0d-16) THEN
                       TruncatedGaussianPrior=-1.0d32
                   ELSE
                       TruncatedGaussianPrior=mu+sigma*SqrtTwo*dierfc(2.d0*(1.d0-r))
                   IF( ( limit > mu .and. TruncatedGaussianPrior > limit ) .or. ( limit < mu .and. &
                       TruncatedGaussianPrior < limit ) ) TruncatedGaussianPrior=-1.0d32
                   ENDIF

            RETURN
                   END FUNCTION TruncatedGaussianPrior

!=======================================================================
! Uniform[0:1]  ->  Cauchy[mean=x0,FWHM=2*gamma]

       FUNCTION CauchyPrior(r,x0,gamma)

            IMPLICIT NONE

                REAL*8                          :: r,x0,gamma,CauchyPrior
                REAL*8 , PARAMETER              :: Pi=3.141592654
     
                    CauchyPrior=x0+gamma*tan(Pi*(r-0.5))

            RETURN
       END FUNCTION CauchyPrior

!=======================================================================
! Uniform[0:1]  ->  LogNormal[mode=a,width parameter=sigma]

       FUNCTION LogNormalPrior(r,a,sigma)

            IMPLICIT NONE

                   REAL*8                          :: r,a,sigma,LogNormalPrior
                   REAL*8                          :: bracket
                   REAL*8 , PARAMETER              :: SqrtTwo=1.414213562d0 
     

                   bracket=sigma*sigma+sigma*SqrtTwo*dierfc(2.d0*r)
                   LogNormalPrior=a*dexp(bracket)

            RETURN
       END FUNCTION LogNormalPrior

!=======================================================================
! Uniform[0:1]  ->  lambda * exp(-lambda * x)[min=a, max=b]

       FUNCTION ExpPrior(r,a,b,lambda)

            IMPLICIT NONE

                REAL*8                          :: ExpPrior,r,a,b,lambda,temp
      
                   IF( a > b ) THEN
      	               temp = a
      	               a = b
      	               b = temp
                   ELSEIF( a == b ) THEN
      	               ExpPrior = a
      	               RETURN
                   ENDIF

                   IF( a < 0d0 ) THEN
	               WRITE(*,*)'One of the limits for exponential prior is less than 0. Aborting'
      	               STOP
                   ENDIF

                   ExpPrior = a - log( 1d0 - r + r * exp( -lambda * (b - a) ) ) / lambda

       END FUNCTION ExpPrior

!=======================================================================
! Uniform[0:1]  ->  x^{-alpha}[min=a, max=b]

       FUNCTION PowerPrior(r,a,b,alpha)

          IMPLICIT NONE

 	        REAL*8                          :: PowerPrior,r,a,b,alpha,p,temp
      
                   IF( a > b ) THEN
      	               temp = a
      	               a = b
      	               b = temp
                   ELSEIF( a == b ) THEN
      	               PowerPrior = a
      	               RETURN
                   ENDIF

                   IF( a < 0d0 ) THEN
	              WRITE(*,*)'ERROR: One of the limits for power prior is less than 0. Aborting'
      	              STOP
                   ENDIF

                   p = 1d0 - alpha
                   PowerPrior = ( ( 1d0 - r ) * ( a ** p ) + r * ( b ** p ) ) ** ( 1d0 / p )

       END FUNCTION PowerPrior

!=======================================================================       
! Inverse of complimentary error FUNCTION in double precision
!
       FUNCTION dierfc(y)
            IMPLICIT REAL*8 (a-h, o-z)
             REAL*8 , PARAMETER              ::infinity=5.0d0
    
                   PARAMETER (qa=9.16461398268964d-01, &
                              qb=2.31729200323405d-01, &
                              qc=4.88826640273108d-01, &
                              qd=1.24610454613712d-01, &
                              q0=4.99999303439796d-01, &
                              q1=1.16065025341614d-01, &
                              q2=1.50689047360223d-01, &
                              q3=2.69999308670029d-01, &
                              q4=-7.28846765585675d-02)
                   PARAMETER (pa=3.97886080735226000d+00, &
                              pb=1.20782237635245222d-01, &
                              p0=2.44044510593190935d-01, &
                              p1=4.34397492331430115d-01, &
                              p2=6.86265948274097816d-01, &
                              p3=9.56464974744799006d-01, &
                              p4=1.16374581931560831d+00, &
                              p5=1.21448730779995237d+00, &
                              p6=1.05375024970847138d+00, &
                              p7=7.13657635868730364d-01, &
                              p8=3.16847638520135944d-01, &
                              p9=1.47297938331485121d-02, &
                              p10=-1.05872177941595488d-01, &
                              p11=-7.43424357241784861d-02)
                   PARAMETER (p12=2.20995927012179067d-03, &
                              p13=3.46494207789099922d-02, &
                              p14=1.42961988697898018d-02, &
                              p15=-1.18598117047771104d-02, &
                              p16=-1.12749169332504870d-02, &
                              p17=3.39721910367775861d-03, &
                              p18=6.85649426074558612d-03, &
                              p19=-7.71708358954120939d-04, &
                              p20=-3.51287146129100025d-03, &
                              p21=1.05739299623423047d-04, &
                              p22=1.12648096188977922d-03)
                   IF (y==0.0) THEN
                       x=infinity
                       GOTO 999
                   ENDIF  
                   z=y
                   IF (y > 1) z=2-y
                       w=qa-log(z)
                       u=sqrt(w)
                       s=(qc+log(u))/w
                       t=1/(u+qb)
                       x=u*(1-s*(0.5d0+s*qd))-((((q4*t+q3)*t+q2)*t+q1)*t+q0)*t
                       t=pa/(pa+x)
                       u=t-0.5d0
                       s=(((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12
                       s=((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2) &
                          *u+p1)*u+p0)*t-z*exp(x*x-pb)
                       x=x+s*(1+x*s)
                   IF (y > 1) x=-x
 999               dierfc=x
            RETURN
       END FUNCTION dierfc
!=======================================================================
END MODULE priors1
