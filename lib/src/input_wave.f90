doubleprecision function solitary_wave(h_shoki,h_wave,t_adjust,t,g)
implicit none
doubleprecision,intent(in)::h_shoki,h_wave,t_adjust,t,g
doubleprecision::c,kpara,eta
!  eta(x,t)=H(sech(k(x-ct))**2
!  H:: h_wave is the maximum wave height(which occurs at x=0 when t=0)
!  k:: kpara
!  c:: wave speed
!  t_adjust:: adjustment for input 
c=sqrt(g*(h_shoki+h_wave)) 
kpara=sqrt(3.0d0*h_wave/4.0d0/h_shoki**3)
eta=h_wave/cosh(c*kpara*(t-t_adjust))**2 
!eta=h_wave/(exp(-c*kpara*(t-t_adjust))+exp(c*kpara*(t-t_adjust)))**2*4.0d0  
solitary_wave=h_shoki+eta
end function solitary_wave
!===================================================
! Functions introduced by William Pringle 2014 below
!===================================================
real*8 function calc_sw_per(h_shoki,h_wave,g)
implicit none
real*8,intent(in) :: h_shoki,h_wave,g
real*8 :: c,k
! Calculates the period of a solitary wave
c = sqrt(g*(h_shoki+h_wave))
k = sqrt(3.0d0*h_wave/(4.0d0*h_shoki**3))
calc_sw_per = 10.0d0 / k / c ! From COULWAVE
end function
!
subroutine sw_eta_vel(eta,vel,period,H,d,t,x,xx,g,kss,kee,z)
implicit none
integer,intent(in) :: kss,kee
real*8,intent(in) :: H,period,t,x,xx,g,d,z(kss:kee)
real*8,intent(out) :: eta, vel(min(kee-kss,2),kss:kee-1)
real*8  :: t_a,kn,c,M,N,Htemp,zm,zz
integer :: k
real*8,parameter :: TT=(2d0/3d0),eps=1d-7
!  eta(x,t) = H(sech(k(x-ct))**2  !Solitary wave theory
!  vel(x,t) = eta*sqrt(g/d)       !Shallow water Linear assumption
!                                 !or from McCowan's theory for 3D   
!  H:: Is the maximum wave height (which occurs at x=0 when t=0)
!  t_a :: adjusted t by half the wave period
!  x:: x adjustment for input
!  k:: k is the wave number
!  c:: c is the wave celerity
kn = sqrt( 3.0d0*H / ( 4.0d0 * d**3 ) )
c  = sqrt( g * ( H + d) )
t_a = t - 0.5d0 * period
eta = H * (cosh(kn * ( x - t_a * c ) ) )**-2
if (kss.eq.kee-1) then
    vel(1,:) = (H * (cosh(kn * ( xx - t_a * c ) ) )**-2) * sqrt( g / d )
else
    M = sqrt( 3.0d0 * H / d ) ! M is guessed as non-dimensional wave number 
    !Now lets iterate to get the real values
    do
        N = TT * ( sin(M*(1.0d0+TT*H/d)) )**2
        !Calculate H =>
        Htemp = (d*N/M) * tan(0.5d0*M*(1.0d0+H/d)) 
        if (abs(Htemp-H).lt.eps) then
            exit !M and N are satisfied
        else
            !Update M
            M = M / sqrt( Htemp / H )
        endif
    enddo
    !Now lets find u and w
    do k = kss,kee-1
        zm = z(k)
        ZZ = 0.5d0 * (z(k) + z(k+1))
        if (z(k) > eta) then
            vel(:,k) = 0.0d0
        else
            vel(1,k) = N * c * ( 1.0d0 + cos(M*(zz+d)/d) * cosh(M*(xx - t_a*c)/d))& 
                     / ( cos(M*(zz+d)/d) + cosh(M*(xx - t_a*c)/d) )**2
            vel(2,k) = N * c * (sin(M*(zm+d)/d) * sinh(M*(x - t_a*c)/d) )         &
                     / ( cos(M*(zm+d)/d) + cosh(M*(x - t_a*c)/d) )**2
        endif
    enddo
endif
endsubroutine sw_eta_vel
!           
subroutine rw_eta_vel(eta,velocity,period,a_wave,h_shoki,k,t,x,g)
implicit none
real*8,intent(in) :: period,a_wave,t,x,k,g,h_shoki
real*8,intent(out) :: eta, velocity
real*8,parameter :: TWOPI = 8d0*atan(1.0d0)
!  eta(x,t) = Hcos(kx-omega*t) : Linear wave theory
!  velocity(x,t) =  eta*sqrt(g/h_shoki) : Shallow water Linear wave theory
!  a:: a_wave is the wave amplitude (which occurs at x=0 when t=0)
!  k:: wave number
!  omega:: angular wave period
!  x:: x adjustment for input
!  period:: wave period
eta = a_wave * cos( k*x - TWOPI*t/period + 0.25d0*TWOPI )
velocity = eta * sqrt( g / h_shoki ) 
endsubroutine rw_eta_vel
!
real*8 function calc_erf_vel(period,t)
implicit none
real*8,intent(in) :: period,t
real*8 :: X
real*8,parameter :: PI_8 = 4d0*atan(1.0d0)
! Calculates the velocity at time t from the error function wavemaker
X = 2.0d0 * PI_8 * t / period
calc_erf_vel = 4.0d0 * sqrt(PI_8) * exp(- X * X) / period
end function
!
real*8 function calc_ns_per(h_wave,size,dt) 
implicit none
integer,intent(in) :: size
real*8,intent(in) :: dt,h_wave(size)
integer :: i,ier,lensav,lenwrk,NH
real*8 :: Emax,freq,r(size)
real*8,allocatable :: wsave(:),work(:),mag(:)
!==============================================
!Calculates the main wave period of the dataset
!==============================================
!
!Perform a FFT of the wave data: (see FFTpack5)
r=h_wave
lensav= size + int ( log ( real ( size, kind = 8 ) ) / log ( 2.0d0 ) ) + 4
lenwrk= size
allocate(wsave(lensav),work(lenwrk))
call dfft1i ( size, wsave, lensav, ier )
call dfft1f ( size, 1, r, size, wsave, lensav, work, lenwrk, ier )
!Now find the largest spectral energy
if (mod(size,2) == 0) then
    NH=size/2-1;  
else
    NH=(size-1)/2;   
endif
allocate(mag(NH))
do i=1,NH
    mag(i)=sqrt(r(2*i)**2+r(2*i+1)**2)
enddo
Emax=maxval(mag)
do i=1,NH
    if (mag(i).eq.Emax) then
        !The frequency at the largest spectral energy
        freq=dfloat(i)/dt/NH
        !Thus the period becomes:
        calc_ns_per=1d0/freq   
    endif
enddo
end function
!
real*8 function calc_k(period,h,g)  
implicit none
real*8,intent(in) :: period,h,g
real*8 :: w,og,test
real*8,parameter :: TWOPI = 8d0*atan(1.0d0),eps=1d-10
!Calculating the wave number k (calc_k) from linear dispersion relation
!Works best for waves that are in shallow or intermediate depths
w=TWOPI/period;og=w**2/g
!First guess using the shallow water limit:
calc_k=sqrt(og/h);goto 124
!For reiteration
123 calc_k=calc_k*sqrt(og/test)
!Now test
124 test=calc_k*tanh(calc_k*h)
if (abs(og-test).lt.eps) return !k is found within the limits of eps
goto 123
end function   
!
real*8 function set_wl(t,dt,size,hwave,order) 
implicit none
integer,intent(in) :: size,order
real*8,intent(in) :: t,dt
real*8,intent(in) :: hwave(size)
integer :: i,j,jm,jp,jn
real*8 :: Interpol_Lagrange,twavem,hwavem,hwavep,x(0:14)
!
!Calculating water level from prescribed wave data
i=ceiling(t/dt)
twavem=dt*(i-1)
if (twavem.eq.t) then
    set_wl=hwave(i)
else
jp=min(7+max(0,8-i),size-i);jm=max(0,jp-14)
jn=0
    do j=jm,jp
        x(jn)=twavem+float(j)*dt
        jn=jn+1
    enddo
    set_wl=Interpol_Lagrange(jn-1,order,t,x(0:jn-1),hwave(i+jm:i+jp))
endif  
!
endfunction set_wl    
!
real*8 function calc_tadjust(wave,depth,cflag)
implicit none
real*8,intent(in) :: wave,depth
character*2,intent(in) :: cflag
real*8 :: k,c,g=9.81d0
real*8,parameter :: PI_8 = 4d0*atan(1.0d0)
!Calculates the adjusted t based on the selected wave theory
if (cflag.eq.'SW') then
    !Solitary Wave
k=(3.0d0*g*wave)/(4.0d0*depth**3)
c=sqrt(g*(wave+depth))
calc_tadjust=2.0d0*PI_8/(k*c)
endif
!
end function calc_tadjust    
!    
real*8 function greensfunc(eta0,h0,hp)
implicit none
real*8,intent(in)::eta0,h0,hp
!Gives the modified wave free surface based on Green's Function
!Eta :: Reference wave free surface
!h0 :: Reference water depth
!hp :: Initial water depth at the location of interest
if (h0.eq.0.0d0) then
    greensfunc=eta0; return
endif
greensfunc=eta0*(hp/h0)**(-0.25d0)
end function
!
real*8 function hn_sanders_obc(h_wave,h_shoki,eta,mom,g)    
implicit none
real*8,intent(in) :: h_shoki,h_wave,eta,mom,g
real*8 :: Dwave,Db,wdiv,w_depth,w_ini
real*8,parameter :: eps = 1.0d-12, COMFLOW = 0.35d0
w_depth = eta + h_shoki
if (h_shoki.lt.0.0d0) then
    !For initially dry grid points
    Dwave = h_wave + sqrt(2d0)  * COMFLOW * w_depth
    w_ini = sqrt(2d0) * COMFLOW * w_depth
else
    Dwave = h_wave + h_shoki; w_ini = h_shoki
endif
Db = ( 4.0d0 * sqrt(g * Dwave  ) - 2.0d0 * sqrt(g * w_ini) +        &
       2.0d0 * sqrt(g * w_depth) - mom / w_depth ) **2 / ( 16d0 * g ) 
if (abs(Db-h_shoki).lt.eps) then
    hn_sanders_obc = 0.0d0
else    
    hn_sanders_obc = Db-h_shoki
endif
end function
!
real*8 function obc(vb,vtp,vtm,vim,vjp,vjm,vbjp,vbjm,f,zkb)
implicit none
real*8,intent(in) :: vb,vtp,vtm,vim,vjp,vjm,vbjp,vbjm,zkb
logical,intent(in) :: f
real*8 :: rx,ry,dvdt,dvdx,dvdy,eps=1.0d-20,g=9.8d0
!Calculates the parameter of interest at the boundary
!(e.g. water level) using an Open boundary formula
!Based on the equations proposed by Marchesiello et al. (2001)
!
!vb :: Current/previous value at the boundary (may be prescribed by a forcing condition)
!vtp :: Current value normally adjacent to vb
!vtm :: Previous value (in time) normally adjacent to vb
!vim :: Current value normally adjacent to vtp
!vjp/vjm :: Current values tangentially adjacent to vtp (in each direction)
!vbjp/vbjm :: Current/previous values tangentially adjacent to vb (in each direction)
!f :: logical that indicates whether a forcing condition is used or not
!
!Calculate the differences
dvdt=vtp-vtm
dvdx=vtp-vim
dvdy=vjp-vjm
    if (abs(dvdt).lt.eps) dvdt=0.0d0
    if (abs(dvdx).lt.eps) dvdx=0.0d0
    if (abs(dvdy).lt.eps) dvdy=0.0d0
    !Converting to upwind method values
    if (dvdt*dvdy.lt.0.0d0) then
        dvdy=vjp-vtm
    elseif (dvdt*dvdy.gt.0.0d0) then
        dvdy=vtm-vjm
    endif
!Calculate normal wave speed
if (dvdx.eq.0.0d0) then
    rx=0.0d0
else
    rx=-dvdt*dvdx/(dvdx**2+dvdy**2)
endif
    !When the wave is travelling into the calc. domain
    if (rx.lt.0.0d0) then
        if (f) obc=vb 
        if (.not.f) obc=vtp
        return
    endif
!Calculate tangential wave speed
if (dvdy.eq.0.0d0) then
    ry=0.0d0
else    
    ry=-dvdt*dvdy/(dvdx**2+dvdy**2)
endif
    !Calculate the new boundary value
    if (ry.gt.0.0d0) then
        obc=(1.0d0/(1.0d0+rx))*(vb+rx*vtp+ry*(vbjm-vb))
    elseif (ry.lt.0.0d0) then
        obc=(1.0d0/(1.0d0+rx))*(vb+rx*vtp+ry*(vb-vbjp))
    else
        obc=(1.0d0/(1.0d0+rx))*(vb+rx*vtp)
    endif
    if (abs(obc-vb).lt.eps) obc=vb
    if (obc.lt.zkb) obc=zkb
!
end function    