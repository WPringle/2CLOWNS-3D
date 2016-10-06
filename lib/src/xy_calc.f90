real*8 function wdiv(upper,lower)
implicit none
real*8,intent(in)::upper,lower
!Divide only when the lower part of the fraction is non-zero
if (lower.eq.0.0d0) then
    wdiv=0.0d0
else
    wdiv=upper/lower
endif
end function wdiv    
!    
real*8 function LI_1(fm,fp,dxm,dxp)
implicit none
real*8,intent(in) :: fm,fp,dxm,dxp
! Function LI_1 : Linear interpolation in 1D
! fm/fp :: value upwind and downwind
! dxm/dxp :: spatial separation upwind and downwind
LI_1=fm+(fp-fm)*dxm/(dxm+dxp)
endfunction LI_1
!
!********************************************************
!*          Lagrange interpolation subroutine           *
!* ---------------------------------------------------- *
!* n is the level of the interpolation ( Ex. n=2 is     *
!* quadratic ). v is the total number of table values.  *
!* X(i), Y(i) are the coordinate table values, Y(i)     *
!* being the dependant variable. The X(i) may be arbi-  *
!* trarily spaced. xx is the interpolation point which  *
!* is assumed to be in the interval  with at least one  *
!* table value to the left, and n to the right. If this *
!* is violated, n will be set to zero. It is assumed    *
!* that the table values are in ascending X(i) order.   *
!*                                                      *
!*     Reference: BASIC Scientific Subroutines, Vol. II *
!*     By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
!*                                                      *
!*                 F90 Version by J.-P. Moreau, Paris.  *
!*                      (www.jpmoreau.fr)               *
!********************************************************
real*8 function Interpol_Lagrange(iv,n,xx,X,Y)  
  !Labels: 100,200,300
  integer i,iv,j,k,n
  real*8  X(0:iv), Y(0:iv)
  real*8  xx,yy
  real*8  XL(0:9)
  ! Check to see if interpolation point is correct
  if (xx < X(0)) goto 100 
  if (xx <= X(iv)) goto 200
  ! An error has been encountered
100 return
  ! Find the relevant table interval
200 i=0
300 i = i + 1
  if (xx > X(i)) goto 300
  i = i - 1
  ! Begin interpolation
  do j = 0, n
    XL(j)=1.d0
  end do
  Interpol_Lagrange=0.d0
  do k = 0, n
    do j = 0, n
      if (j.eq.k) goto 400
      XL(k)=XL(k)*(xx-X(j+i))/(X(i+k)-X(j+i))
400 end do
    Interpol_Lagrange=Interpol_Lagrange+XL(k)*Y(i+k)
  end do
  return
endfunction Interpol_Lagrange    
!    
real*8 function lat_to_m(latp,latm)
implicit none
real*8,intent(in) :: latp,latm
real*8 :: alat,rlat,m
real*8,parameter :: PI_8 = 4d0*atan(1.0d0)
! lat_to_m = latitude difference in meters
! latp = latitude at the north cell boundary
! latm = latitude at the south cell boundary
! Reference: American Practical Navigator, Vol II, 1975 Edition, p 5 
alat = 0.5d0*(latp+latm)
rlat = alat*PI_8/180d0
m = 111132.09d0*rlat - 566.05d0*cos(2d0*rlat) + 1.2d0*cos(4d0*rlat)
lat_to_m = (latp-latm)*m
endfunction lat_to_m
!    
real*8 function lon_to_m(dlon,alat)
implicit none
real*8,intent(in) :: dlon,alat
real*8 :: rlat,p
real*8,parameter :: PI_8 = 4d0*atan(1.0d0)
! lon_to_m = longitude difference in meters
! dlon = longitude difference in degrees
! alat = average latitude between the two fixes
! Reference: American Practical Navigator, Vol II, 1975 Edition, p 5 
rlat = alat*PI_8/180d0
p = 111415.13d0 * cos(rlat) - 94.55d0 * cos(3d0*rlat)
lon_to_m = dlon*p
endfunction lon_to_m     
!
real*8 function round(num,dec)
implicit none
real*8,intent(in) :: num
integer,intent(in) :: dec
real*8 :: A
!Round a number to the desired number of decimal places
! num :: the number
! dec :: number of decimal places for rounding
A = 10d0**dec
round = dfloat(nint(num*A))/A   ! must be dfloat
!write(6,*) num,a,num*a,dfloat(nint(num*A)),round
!
endfunction round
!Filtering routines below:
real*8 function Shapiro_filt(val,point,nl,nr)  
implicit none
integer,intent(in) :: point,nl,nr
real*8,intent(in) :: val(-nl:nr)
integer :: nm,e
real*8,parameter,dimension(4,6) :: &
c=reshape((/4.0d0,16.0d0,64.0d0,256.0d0&
,2.0d0,10.0d0,44.0d0,186.0d0,1.0d0,4.0d0,15.0d0,56.0d0&
,0.0d0,-1.0d0,-6.0d0,-28.0d0,0.0d0,0.0d0,1.0d0,8.0d0&
,0.0d0,0.0d0,0.0d0,-1.0d0/),shape(c))
!Shapiro filter to desired number of points (3,5,7 or 9)
if (nl.ne.nr) then
    write(6,*) 'Error in Shapiro filter'
    return
endif
Shapiro_filt=0.0d0
do nm=-nl,nr 
    e=abs(nm)+2
    Shapiro_filt=Shapiro_filt+val(nm)*c(point/2,e)
enddo
Shapiro_filt=Shapiro_filt/c(point/2,1)
return
endfunction Shapiro_filt    
!    
!real*8 function SavGol_filt(val,point,m,nl,nr)  
!implicit none
!integer,intent(in) :: point,m,nr,nl
!real*8,intent(in) :: val(-nl:nr)
!integer :: nm,e
!real*8 :: c(point)
!!
!!Savitzky-Golay filter to desired number of points
!!5-point uses the second order polynomial (m=2)
!!9-point uses the fourth order polynomial (m=4)
!call savgol(c,point,nl,nr,0,m,2*(m+1))
!SavGol_filt=0.0d0
!do nm=-nl,nr 
!    e=mod(point-nm,point)+1
!    SavGol_filt=SavGol_filt+val(nm)*c(e)
!enddo
!return
!endfunction SavGol_filt
!!
!subroutine savgol(c,np,nl,nr,ld,m,Lwork)
!integer,intent(in) :: ld,m,nl,np,nr,Lwork
!real*8,intent(out) :: c(np)
!integer,parameter :: MMAX=6
!!USES dgels
!!Returns in c(1:np), in wrap-around order (N.B.!) consistent with the argument respns
!!in SavGol_filt, a set of Savitzky-Golay lter coefficients. nl is the number of leftward
!!(past) data points used, while nr is the number of rightward (future) data points, making
!!the total number of data points used nl+nr+1. ld is the order of the derivative desired
!!(e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also
!!equal to the highest conserved moment; usual values are m = 2 or m = 4.
!!
!!Original By: Numerical Recipes in Fortran 77: The art of scientific computing Copywright (C) 1986-1992
!!(Has been modified for Intel Fortran and for the matrix solving routine by William Pringle 2013)
!integer :: imj,ipj,k,kk,mm,info
!real*8 :: fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1),work(Lwork)
!if (np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m) then
!    write(6,*) 'bad args in savgol'
!    return
!endif
!do ipj=0,2*m                        !Set up the normal equations of the desired least-
!    sum=0.0d0                       !squares fit.
!    if(ipj.eq.0) sum=1.0d0
!    do k=1,nr
!        sum=sum+dfloat(k)**ipj
!    enddo
!    do k=1,nl
!        sum=sum+dfloat(-k)**ipj
!    enddo
!    mm=min(ipj,2*m-ipj)
!    do imj=-mm,mm,2
!        a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
!    enddo
!enddo
!b=0.0d0;b(ld+1)=1.0d0                !Solve the system of equations using LAPACK dgels
!call dgels('N',MMAX+1,m+1,1,a(1:MMAX+1,1:m+1),MMAX+1,b,MMAX+1,work,Lwork,info)
!            if (info.ne.0) then
!                write(6,*) 'Matrix calc. in savgol unsuccessful, info=',info
!                return
!            endif                     
!c=0.0d0                              !Zero the output array (it may be bigger than number of coefficients).
!do k=-nl,nr                          !Each Savitzky-Golay coefficient is the dot product                                    
!sum=b(1)                             !of powers of an integer with the inverse matrix row
!fac=1.0d0
!    do mm=1,m
!        fac=fac*k
!        sum=sum+b(mm+1)*fac
!    enddo
!kk=mod(np-k,np)+1                    !Store in wrap-around order.
!c(kk)=sum
!enddo
!endsubroutine savgol    