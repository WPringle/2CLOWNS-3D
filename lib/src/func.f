c********************************************************
      FUNCTION RHO0T_MG(Ten,cc,rhos)
c********************************************************
      implicit real*8 (a-h,o-z)
         rho0 =
     1       4.0704d-5  *Ten**3
     1      -7.7617d-3  *Ten**2
     1      +5.5301d-2  *Ten
     1      +999.91d0

CC ppm cc
c     RHO0T_MG=rho0+(rhos-rho0)*cc*1.0d-6
CC ppm cc

CC mg/L cc
      RHO0T_MG=rho0+(rhos-rho0)*cc*1.0d-3/rhos
CC mg/L cc
      return
      end


c********************************************************
      FUNCTION RHO0T2(Ten,cc,rhos)
c********************************************************
c
      implicit none
      doubleprecision::   RHO0T2
      doubleprecision,intent(in)::   Ten,cc,rhos     
      doubleprecision::   rho0,r0rs
 

         rho0 =
     1      (4.0704d-8  *Ten**3
     1      -7.7617d-6  *Ten**2
     1      +5.5301d-5  *Ten
     1      +0.99991d0)
     1    *1000.d0
            r0rs=(rhos/rho0)-1.0d0
c           r0rs=1.0d0-(rho0/rhos)
            RHO0T2=rho0*(1.0d0+r0rs*cc/1.0d6)
      return
      end

c********************************************************
      FUNCTION RHO0T(Ten)
c********************************************************
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8  RHO0T
      real*8 data(0:50)/0.99984d0,0.99990d0,0.99994d0,
     1     0.99996d0,0.99997d0,0.99996d0,
     1     0.99994d0,0.99990d0,0.99985d0,0.99978d0,0.99970d0,0.99961d0,
     1     0.99949d0,0.99938d0,0.99924d0,0.99910d0,0.99894d0,0.99877d0,
     1     0.99860d0,0.99841d0,0.99820d0,0.99799d0,0.99777d0,0.99754d0,
     1     0.99730d0,0.99704d0,0.99678d0,0.99651d0,0.99623d0,0.99594d0,
     1     0.99565d0,0.99534d0,0.99503d0,0.99470d0,0.99437d0,0.99403d0,
     1     0.99368d0,0.99333d0,0.99297d0,0.99259d0,0.99222d0,0.99183d0,
     1     0.99144d0,0.99104d0,0.99063d0,0.99021d0,0.98979d0,0.98936d0,
     1     0.98893d0,0.98849d0,0.98804d0/
c
c
      if((Ten .LE. 0.0d0) .or. (Ten .GE. 50.0)) then
       rho0T=0.0d0
      else
c
      iT1=nint(Ten)
      rho0T=(data(iT1)+(data(iT1+1)-data(iT1))*(Ten-dfloat(iT1))
     1      )*1.0d3      
      end if

      return
      end

        doubleprecision function rho0T2S(tt,cc) !‰–•ª”Z“x
        doubleprecision::tt,cc
        rho0T2S=1028.14d0-0.0735d0*tt-0.00469*tt**2+(0.802d0-0.002d0*tt)*(cc-35.0d0)
        end function