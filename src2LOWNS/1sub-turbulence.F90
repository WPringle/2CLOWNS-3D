!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-turbulence.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Updates the turbulence through the k-epsilon method
subroutine cal_turbulence
      ! Calculate the production term from new velocities
      call pro_k !;write(6,*) "end product"
      
      ! Update k and e through their respective transport equations
      call cal_ke_ex !;write(6,*) "end cal_ke"
    
      ! Update the turbulent viscosity from k and e
      call kado !;write(6,*) "kado"

endsubroutine cal_turbulence    

subroutine pro_k
    use interface_list,only: bibun, vbnd, vbnd_1, wdiv
    use variables, only: inns, inne, ifln, n, t, nu, js, je, ieee, mu, rho0,  &
                         ZERO, HALF, TWO, ONE, ONEHALF, THIRD, SQRSIX,        &
                         pro_lim, vt_op
    use arrays, only: a, nff, un, dx, fb, mn, fln, in, jn, kn,                &
                      pro, prop, d, r, SS, AsUs
    implicit none    
    integer :: i,k,j,nn,im,ip,jm,jp,km,kp,nnf
    real*8 :: uxp,uc,uxm
    real*8 :: vyp,vc,vym
    real*8 :: wzp,wc,wzm
    real*8 :: nu1
    real*8 :: ukp,ukm,ujp,ujm
    real*8 :: vip,vim,vkp,vkm
    real*8 :: wip,wim,wjp,wjm
    real*8 :: dzc,dzp,dzm  
    real*8 :: u2,v2,w2
    real*8 :: W, SS3, Phi, As, U_star ! For realizable k-eps Oct 5 2016
    real*8 :: dudy,dudz,dvdx,dvdz,dwdx,dwdy,QQ2 = ZERO
    real*8 :: dyc1,dym1,dyp1,dzc1,dzm1,dzp1,SS2
    real*8,parameter :: C_lim = 10d0

    ieee = 0
 !  akrateS=0.8d0
#ifdef DAKU
    mu=nu*rho0 
#else
    nu1 = nu
#endif
!$omp parallel do private(nn,i,k,j,im,ip,jm,jp,km,kp,dzc,dzp,dzm&
!$omp ,vyp,vym,uxp,uxm,wzm,wzp,uc,vc,wc,ukm,ukp,dudz,ujm,ujp&
!$omp ,DUDY,dvdz,dzc1,dzm1,dzp1,vkm,vkp,vim,vip,dvdx,wim,wip&
!$omp ,dwdx,wjp,wjm,dyc1,dym1,dyp1,DWDY,u2,v2,w2,SS2,qq2    &
!$omp ,As,U_star,SS3,W,Phi)
PRO_LOOP: DO nnf = 1,ifln
    !Map to the real cell number
    nn = fln(nnf)
#ifdef BUGWIN
    if (nn.eq.459692) then
        continue
    endif
#endif
    !Get integers
    i = in(nn) ; j=jn(nn) ; k = kn(nn)
    !Get adjacent cell integers
    im = i-1 ; ip = i+1 ; jm = j-1 ; jp = j+1 ; km = k-1 ; kp = k+1
    !If an object cell or air cell production term is zero 
    if (nff(nn)%f.eq.-1.or.nff(nn)%f.eq.0.or.fb(nn).eq.ZERO) then
        pro(nn) = ZERO ; cycle
    endif
    !CC                    
    !CC  REYNOLDS STRESS   
    !CC                    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !William, 08 Aug 2014  
    !Added new function vbnd_1 above to tidy up code
    !FREESURF preprocessor is read in the function
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !Determine velocities at the cell boundary
    vyp = vbnd_1(un(1,mn(i,k,jp)),un(1,nn),a(1,mn(i,k,jp)),1)
    vym = vbnd_1(un(1,nn),un(1,mn(i,k,jp)),a(1,nn),1)
    uxp = vbnd_1(un(0,mn(ip,k,j)),un(0,nn),a(0,mn(ip,k,j)),0)
    uxm = vbnd_1(un(0,nn),un(0,mn(ip,k,j)),a(0,nn),0)
    wzp = vbnd_1(un(2,mn(i,kp,j)),un(2,nn),a(2,mn(i,kp,j)),2)
    wzm = vbnd_1(un(2,nn),un(2,mn(i,kp,j)),a(2,nn),2)
    !Calculate velocities at cell center
    uc = ( uxm + uxp ) * HALF
    vc = ( vym + vyp ) * HALF
    wc = ( wzp + wzm ) * HALF
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !William, 08 Aug 2014  
    !Altered to include better performance for sloping staircase chikei
    !Put under FREESURF preprocessor and read in vbnd function
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%% Calculation for X direction %%%%%%%%%%%%%%%%%%%%%%%%
    dzc = dx(2,k)
    if (mn(i,kp,j).eq.0) then
        dzp = dzc
    else
        dzp = dx(2,kp)
    endif
    if (mn(i,km,j).eq.0) then
        dzm = dzc
    else
        dzm = dx(2,km)
    endif
    !Determine velocities at neighbouring cells
    ukm = vbnd(nff(mn(i,km,j))%f,nff(mn(i,km,j))%b,uc,un(0,mn(i,km,j)),un(0,mn(ip,km,j)), &
               dzc,dzm,a(2,nn),a(0,mn(ip,km,j)),a(0,mn(i,km,j)))
    ukp = vbnd(nff(mn(i,kp,j))%f,nff(mn(i,kp,j))%b,uc,un(0,mn(i,kp,j)),un(0,mn(ip,kp,j)), &
             dzc,dzp,a(2,mn(i,kp,j)),a(0,mn(ip,kp,j)),a(0,mn(i,kp,j)))
    ujm = vbnd(nff(mn(i,k,jm))%f,nff(mn(i,k,jm))%b,uc,un(0,mn(i,k,jm)),un(0,mn(ip,k,jm)), &
             dx(1,j),dx(1,jm),a(1,nn),a(0,mn(ip,k,jm)),a(0,mn(i,k,jm)))
    ujp = vbnd(nff(mn(i,k,jp))%f,nff(mn(i,k,jp))%b,uc,un(0,mn(i,k,jp)),un(0,mn(ip,k,jp)), &
             dx(1,j),dx(1,jp),a(1,mn(i,k,jp)),a(0,mn(ip,k,jp)),a(0,mn(i,k,jp)))
    !Calculate the first-order derivatives
    DUDZ = BIBUN(ukp,uc,ukm,dzp,dzc,dzm) 
    DUDY = BIBUN(ujp,uc,ujm,dx(1,jp),dx(1,j),dx(1,jm)) 
    !
    !%%%%%%%%%%% Calculation for Y direction %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nff(mn(i,k,j))%f.eq.2.and.abs(nff(mn(i,k,j))%b).eq.2&
        .and.dx(1,j).gt.dx(2,k)*1.5d0) then
        DVDZ = ZERO
    else
        dzc1 = dx(2,k)
        if(nff(mn(i,km,j))%f.eq.-1) then
            dzm1 = dzc1
        else
            dzm1 = dx(2,km)
        endif
        if(nff(mn(i,kp,j))%f.eq.-1) then
            dzp1 = dzc1
        else
            dzp1 = dx(2,kp)
        endif
        !William, 08 Aug 2014  
        vkm=vbnd(nff(mn(i,km,j))%f,nff(mn(i,km,j))%b,vc,un(1,mn(i,km,j)),un(1,mn(i,km,jp)),   &
                 dzc1,dzm1,a(2,nn),a(1,mn(i,km,jp)),a(1,mn(i,km,j)))       
        vkp=vbnd(nff(mn(i,kp,j))%f,nff(mn(i,kp,j))%b,vc,un(1,mn(i,kp,j)),un(1,mn(i,kp,jp)),   &
                 dzc1,dzp1,a(2,mn(i,kp,j)),a(1,mn(i,kp,jp)),a(1,mn(i,kp,j)))
        !Calculate the first-order derivative
        DVDZ=BIBUN(vkp,vc,vkm,dzp1,dzc1,dzm1)
    endif
    !Determine velocities at neighbouring cells
    vim=vbnd(nff(mn(im,k,j))%f,nff(mn(im,k,j))%b,vc,un(1,mn(im,k,j)),un(1,mn(im,k,jp)),       &
             dx(0,i),dx(0,im),a(0,nn),a(1,mn(im,k,jp)),a(1,mn(im,k,j)))  
    vip=vbnd(nff(mn(ip,k,j))%f,nff(mn(ip,k,j))%b,vc,un(1,mn(ip,k,j)),un(1,mn(ip,k,jp)),       &
             dx(0,i),dx(0,ip),a(0,mn(ip,k,j)),a(1,mn(ip,k,jp)),a(1,mn(ip,k,j)))
    !Calculate the first-order derivative
    DVDX=BIBUN(vip,vc,vim,dx(0,ip),dx(0,i),dx(0,im))
    !
    !%%%%%%%%%%%% Calculation for Z direction %%%%%%%%%%%%%%%%%%%%%%
    !Determine velocities at neighbouring cells
    wim=vbnd(nff(mn(im,k,j))%f,nff(mn(im,k,j))%b,wc,un(2,mn(im,k,j)),un(2,mn(im,kp,j)),       &
             dx(0,i),dx(0,im),a(0,nn),a(2,mn(im,kp,j)),a(2,mn(im,k,j)))  
    wip=vbnd(nff(mn(ip,k,j))%f,nff(mn(ip,k,j))%b,wc,un(2,mn(ip,k,j)),un(2,mn(ip,kp,j)),       &
             dx(0,i),dx(0,ip),a(0,mn(ip,k,j)),a(2,mn(ip,kp,j)),a(2,mn(ip,k,j)))
    !Calculate the first-order derivative
    if (nff(mn(im,k,j))%f.eq.-1.and.nff(mn(im,kp,j))%f.eq.1.and.nff(mn(im,km,j))%f.eq.1) then
        DWDX=((wc+wip)*HALF-(wc+wim)*HALF)/dx(0,i)
    else
        DWDX=BIBUN(wip,wc,wim,dx(0,ip),dx(0,i),dx(0,im))         
    end if 
    !Determine velocities at neighbouring cells     
    wjm=vbnd(nff(mn(i,k,jm))%f,nff(mn(i,k,jm))%b,wc,un(2,mn(i,k,jm)),un(2,mn(i,kp,jm)),       &
             dx(1,j),dx(1,jm),a(1,nn),a(2,mn(i,kp,jm)),a(2,mn(i,k,jm)))
    wjp=vbnd(nff(mn(i,k,jp))%f,nff(mn(i,k,jp))%b,wc,un(2,mn(i,k,jp)),un(2,mn(i,kp,jp)),       &
             dx(1,j),dx(1,jp),a(1,mn(i,k,jp)),a(2,mn(i,kp,jp)),a(2,mn(i,k,jp)))  
    !
    dyc1=dx(1,j)
    if(nff(mn(i,k,jm))%f.eq.-1) then
        dym1=dyc1
    else
        dym1=dx(1,jm)
    endif
    if(nff(mn(i,k,jp))%f.eq.-1) then
        dyp1=dyc1
    else
        dyp1=dx(1,jp)
    endif
    !Calculate the first-order derivative
    DWDY = BIBUN(wjp,wc,wjm,dyp1,dyc1,dym1)
    !Calculate the remaining first-order derivatives
    u2 = ( Uxp - Uxm ) / dx(0,i)
    v2 = ( vyp - Vym ) / dx(1,j)       
    w2 = ( wzp - Wzm ) / dzc
    !Evaulate the modulus of the mean rate-of-strain tensor, S
    SS2 = (DUDZ + DWDX )**2 + (DUDY + DVDX )**2 + ( DVDZ + DWDY )**2                       &
        + TWO * ( u2**2 + v2**2 + w2**2 )
    !
    if (vt_op.eq.2.or.vt_op.eq.3) then
        !For Smagorinsky model or realizable k-eps
        SS(nn) = sqrt(SS2)
        if (vt_op.eq.2) cycle ! Obtained strain-rate, cycle
        if (vt_op.eq.3) then 
            !For realizable k-eps only
            SS3      = HALF * HALF * ( (DUDZ + DWDX )**2 + (DUDY + DVDX )**2 +             &
                       ( DVDZ + DWDY )**2 ) + u2**2 + v2**2 + w2**2
            W        = TWO * SQRT(TWO) * wdiv( SS3 , SS(nn)**3 ) 
            Phi      = THIRD * acos(max(-ONE,min(SQRSIX*W,ONE)))
            As       = SQRSIX * cos(phi)
            QQ2      = (DUDZ - DWDX )**2 + ( DUDY - DVDX )**2 + ( DVDZ - DWDY )**2 
            U_star   = sqrt( HALF * SS2 + HALF * QQ2 )
            AsUs(nn) = As * U_star
        endif
    endif
    !
    if (pro_lim.eq.1.or.pro_lim.eq.3) then 
    !<<<<<<<<<<<<<<<<< This is the Kato-Launder modification <<<<<<<<<<<<<<<<<<
        !Evaluate QQ2 (the vorticity)
        QQ2 = (DUDZ - DWDX )**2 + ( DUDY - DVDX )**2 + ( DVDZ - DWDY )**2
        !Evaluate the production term 
        PRO(nn) = D(nn) * sqrt(SS2) * sqrt(QQ2)
    else
    !<<<<<<<<<<<<<<<< This is the standard k-eps approach <<<<<<<<<<<<<<<<<<<<<
        !Evaluate the production term 
        PRO(nn) = SS2 * D(nn)
    endif
    !Limit increase in production term
    if( (n.ge.2.or.t.gt.ZERO).and.prop(nn).ne.ZERO) then
        if (abs(pro(nn)-prop(nn)).gt.prop(nn)) then    
           pro(nn) = prop(nn) + dsign(ONE,pro(nn)-prop(nn)) * prop(nn)
        endif    
    endif       
    if (pro_lim.eq.2.or.pro_lim.eq.3) then
    ! Set a limit to production term to be less than ten times the dissipation
    ! see: Menter et al. (2003). Ten Years of Industrial Experience with the 
    !      SST Turbulence Model. Turbulence Heat and Mass Transfer 4. 
        pro(nn) = min( C_lim * r(nn)%e, pro(nn) ) 
    endif
enddo PRO_LOOP
!$omp end parallel do
end subroutine pro_k   
   
function vbnd(nfm,nfbm,wc,wmm,wmp,dxc,dxm,az,axp,axm)
use variables, only: ZERO, HALF, ONE
use interface_list, only: ubnds, ubnds2
implicit none
real*8 :: vbnd
real*8,intent(in)  :: wc,wmp,wmm,dxc,dxm,az,axm,axp
integer,intent(in) :: nfm,nfbm
real*8 :: wf = ONE
if(nfm.eq.-1) then
    if (nfbm.eq.-1) then
        vbnd = wc
    elseif (nfbm.eq.2) then 
        vbnd = HALF * ( wmm + wmp )
    elseif (nfbm.eq.10) then 
        vbnd = ubndS2(wf,wc,dxc,dxm)
    else
        vbnd = ubndS(wf,wc,dxc,dxm)
    endif
!William, 13 June 2014  
!Altered to include better performance for sloping staircase chikei
!#ifdef FREESURF
elseif (axp.eq.ZERO) then
        vbnd = wmm  
elseif (axm.eq.ZERO) then
        vbnd = wmp  
!#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elseif (az.eq.ZERO) then
    vbnd = ubndS(wf,wc,dxc,dxm)
elseif (nfm.eq.0) then   !!!!!!　地形が階段状で、水面の乱流生成が過大になる場合
    vbnd = wc              !!!!!!  地形が階段状で、水面の乱流生成が過大になる場合
else
    vbnd = HALF * ( wmm + wmp )
endif
endfunction vbnd

!William, 04 July 2014  
!Added new function vbnd_1 to tidy up code in PRODUCT subroutine
function vbnd_1(v,va,ax,iflag)
use variables, only: ZERO
implicit none
real*8 :: vbnd_1
real*8,intent(in)  :: v, va, ax
integer,intent(in) :: iflag
!
if (iflag.eq.2) then
    vbnd_1 = v * ax
else
    vbnd_1 = v
endif
!#ifdef FREESURF
if(ax.eq.ZERO) vbnd_1 = va
!#endif
end function vbnd_1
   
subroutine KADO
    use variables, only: inns, inne, iob, akmin, cm, emin, ZERO, ONE, &
                         D_rlim, vt_op, nu, THIRD, c_smag, A0
    use arrays, only : in, jn, kn, r1, d, nff, mn, dx, SS, AsUs
    implicit none
    integer :: i, k, j, nn
    real*8  :: Rt, f, delta, C_nu
#ifdef FURYOKU
    real*8  :: bfp,rhop,rhom,zp,zm,ak1bye1
!$omp parallel do private(i,k,j,bfp,rhop,rhom,zp,zm,ak1bye1)
#else
!$omp parallel do private(i,k,j,Rt,f,delta,C_nu)
#endif
DO nn = inns,inne
    i = in(nn) ; j = jn(nn) ; k = kn(nn)
#ifdef BUGWIN
    if (nn.eq.302285) then
        continue
    endif
#endif
    if (nff(nn)%f.eq.0.and.nff(mn(i,k-1,j))%f.ne.2) cycle
    if (nff(nn)%f.eq.iob.and.nff(nn)%b.ge.3) cycle

    if (r1(nn)%e.eq.ZERO) then
        r1(nn)%k = akmin
        r1(nn)%e = emin
        d(nn) = cm * akmin * akmin / emin
    else
#ifdef FURYOKU
        if(nf(i,k+1,j).ge.1) then
        rhop=rho(mn(i,k+1,j));zp=(z(k+1)+z(k+2))*HALF
        else
        rhop=rho(nn);zp=(z(k+1)+z(k))*HALF
        endif
        if(nf(i,k-1,j).ge.1) then
        rhom=rho(mn(i,k-1,j));zm=(z(k)+z(k-1))*HALF
        else
        rhom=rho(nn);zm=(z(k+1)+z(k))*HALF
        endif

        if(zp-zm.gt.ZERO.and.rho(nn).ne.ZERO.and.r1(nn)%e.ne.ZERO) then  !.and.rhom-rhop.gt.8.0d0) then
        ak1bye1=r1(nn)%k/r1(nn)%e
        Bfp=-g/rho(nn)*min( ZERO,rhop-rhom)/(zp-zm)*ak1bye1**2   !/r1(nn)%e**2*ak1(nn)**2
        continue
        else
        bfp=ZERO
        endif
        if(1.0d0/(1.0+bfp).lt.HALF) then
        continue
        endif
        frate(nn)=1.0d0/(1.0+0.1d0*bfp)
        D(nn)=CM*r1(nn)%k**2/r1(nn)%e*1.0d0/(1.0+stpara*bfp)
#else
        if (vt_op.eq.0) then
            !Use the standard linear turbulent eddy viscosity relation
            D(nn) = CM * r1(nn)%k * r1(nn)%k / r1(nn)%e
        elseif (vt_op.eq.1) then
            !Use the adjustment for low-turbulence regions, see:
            !(Bradford, 2000. Numerical Simulation of Surf Zone Dynamics,
            ! Journal of Waterway, Port, Coastal, and Ocean Engineering)
            Rt = r1(nn)%k * r1(nn)%k / r1(nn)%e / nu
            if (Rt.ge.1d5) then
                f = ONE
            else
                f = exp(-3.4d0/(ONE + Rt/50d0)**2)
            endif
            D(nn) = CM * f * r1(nn)%k * r1(nn)%k / r1(nn)%e
        elseif (vt_op.eq.2) then
            ! Smagorinsky model which only depends on the strain rate
            ! and the size of the cell
            delta = ( dx(0,i) * dx(1,j) * dx(2,k) )**(THIRD)
            D(nn) = SS(nn) * (c_smag * delta)**2 
        elseif (vt_op.eq.3) then
            ! Realizable k-epsilon model
            ! Shih et al. (1993). A realizable reynolds stress
            ! algebraic equation model
            C_nu   = min( CM, ONE / ( A0 + AsUs(nn) * r1(nn)%k / r1(nn)%e ) )
            D(nn)  = C_nu * r1(nn)%k * r1(nn)%k / r1(nn)%e
        else
            !Use nonlinear algebraic Reynolds stress model, see:
            !(Lin & Liu, 1998. Turbulence transport, vorticity dynamics, 
            !and solute mixing under plunging breaking waves in surf zone,
            !Journal of Geophysical Research)
            stop 'Under construction: Cannot use vt_op > 3'
        endif
        ! Limit turbulent viscosity to ratio of molecular viscosity
        if (D(nn).gt.D_rlim*nu) D(nn) = D_rlim * nu
#endif 
    endif
enddo
!$omp end parallel do
endsubroutine kado   
    
function ustar(uc,dy,ks_r)
    use variables, only: nu, HALF, ZERO, karman, LITTLE, FIVE
    implicit none
    real*8 :: ustar 
    real*8,intent(in) :: uc, dy, ks_r
    real*8 :: ksstar, ksstar1, ustar1, eps, E, yplus, ustar_v
    real*8 :: a = 8.74d0, E_s = 7.8d0, E_r = 0.36d0, R_c = 26.6d0,            &
              SMOOTH = 5d0, ROUGH = 70d0
    !Guess ustar from the linear log law (in viscous sub layer)
    ustar = sqrt(abs(uc)*nu/dy)
    ! Check for validility...
    yplus = dy * ustar / nu
    if (yplus < FIVE) return !Definitely inside the viscous sub-layer. 
    !Use log law for larger values of y+
    eps = abs(uc) * LITTLE
    !Log - law => u_c = (ustar/karman)*ln(dy/yk) 
    ! Check hydraulic smoothness
    ustar_v = ustar
    ksstar = ks_r * ustar / nu
    do ! Iterate until ksstar and ustar converge
        if (ksstar.le.SMOOTH) then
            ! Use log law considering hydraulically smooth surface
            do
                ustar1 = ustar
                ustar = abs(uc) * karman / log(dy * ustar1 * E_s / nu)
                if (abs(ustar1 - ustar)/ustar.lt.eps) exit
            enddo
        elseif (ksstar.ge.ROUGH) then
            ! Use log law considering hydraulically rough surface
            ustar = abs(uc) * karman / log(dy * R_c / ks_r )
        else
            ! Transitional regime (assume E varies linearly from E_s to E_r)
            E = E_s + (ksstar - SMOOTH) / (ROUGH-SMOOTH) * (E_r - E_s)
            do 
                ustar1 = ustar
                ustar = abs(uc) * karman / log(dy * ustar1 * E / nu)
                if (abs(ustar1 - ustar)/ustar.lt.eps) exit
            enddo
        endif
        ksstar1 = ks_r * ustar / nu
        ! Exit if ksstar has converged
        if (ksstar.eq.ZERO) exit
        if ( abs(ksstar1 - ksstar)/ksstar.lt.eps.or.                          & 
            (ksstar.le.SMOOTH.and.ksstar1.le.SMOOTH).or.                      &
            (ksstar.ge.ROUGH.and.ksstar1.ge.ROUGH) ) exit 
        ! Loop if ksstar has not converged
        ksstar = ksstar1         
    enddo
    ! Check y+
    !yplus = dy * ustar / nu
    !if (yplus > 500) then
    !    write(6,*) 'yplus =',yplus,'ksstar',ksstar 
    !endif
    !Check maximum from linear or log-law
    ustar = max(ustar,ustar_v)
    !================== OLD power law =========================================
    ! Guess ustar using Prandtl's 1/7 law 
    ! (for small Reynolds numbers using Blasius's resistance formular argument)
    ! See: Schlichting, Herrmann, 1979. Boundary layer Theory, 7th ed. pg 601 
    !ustar = (abs(uc)/a)**(7./8.)*(nu/dy)**(1./8.)
    !R = abs(uc) * dy / nu  
    !7d4) return ! Limit as suggested by 水理学 (1967) pg. 66
    !if (R.lt.50d0) return ! Includes when uc is zero 
    !(avoiding divide by zero in log-law)
end function 

function rbnd(rc,vvv,dxb,ks_in)
    use type_list,only: rancomp
    use interface_list,only: ustar, wdiv
    use variables, only: akrate, erate, akini, eini, cm, dt,                  &
                         karman, HALF, TWO, ks_r, ZERO, nu, ONE 
    implicit none
    type(rancomp) :: rbnd
    type(rancomp),intent(in) :: rc
    real*8,intent(in) :: vvv, dxb, ks_in
    real*8 :: ust, dy, vand
    !
    if (ks_in.ge.ZERO) then
        dy = HALF * dxb
        ! Wall function,: ks_r = 0: smooth wall, ks_r > 0: rough wall
        ! Get the friction velocity from log-law
        ust = ustar(vvv,dy,ks_in)
        ! Now we can find k and e from the assumption 
        ! that Pro_k = rho*e in log-law region:
        rbnd%k = ust**2 / sqrt(cm)     
        ! decreasing Van driest function
        vand   = ONE - exp(-ust*dy/nu/26d0) 
        rbnd%e = wdiv( ust**3 / karman / dy , vand)
        !rbnd%e = ust**3 / karman / dy
        rbnd%d = ust * karman * dy * vand !(Since d = cm*k^2/e)    
    else
        !No wall function (no gradient or some percentage gradient)
        rbnd%k = max( rc%k * akrate , akini )  
        rbnd%e = max( rc%e * erate ,  eini  )   
        rbnd%d = cm * rbnd%k**2 / rbnd%e  
    endif
    rbnd%x = dxb
endfunction
!
! Second-order TVD William April 8 2015
subroutine TVD_ke(fx,ae,aw,ue,uw,fE,fp,fW,fWW,dxx)
    use type_list, only: rancomp, ranke
    use variables, only: HALF, TWO
    implicit none
    real*8,intent(in) :: ae, aw, ue, uw, dxx
    type(rancomp),intent(in) :: fE, fp, fW, fWW
    type(ranke),intent(out)  :: fx
    real*8                   :: Phie, Phiw, MPe, MPw 
    type(ranke)              :: Le, Lw
    ! New Second-order TVD scheme
    ! Get limiter values
    call limiter_ke(Le,fp,fW,fE)
    call limiter_ke(Lw,fW,fWW,fp)
    !Get multiples for non-uniform grid
    MPe = (HALF*dxx) / (fE%x-fp%x)
    MPw = (HALF*dxx) / (fp%x-fW%x)
    !%%%%%%%%%%%!
    !     k     !
    !%%%%%%%%%%%!
    ! Get east and west values
    Phie = fp%k + MPe * Le%k * (fE%k - fp%k)
    Phiw = fW%k + MPw * Lw%k * (fp%k - fW%k)
    ! Now minus from each other and divide by dx
    fx%k = (Phie * ae * ue - Phiw * aw * uw) / dxx
    !%%%%%%%%%%%!
    !     e     !
    !%%%%%%%%%%%!
    ! Get east and west values
    Phie = fp%e + MPe * Le%e * (fE%e - fp%e)
    Phiw = fW%e + MPw * Lw%e * (fp%e - fW%e)
    ! Now minus from each other and divide by dx
    fx%e = (Phie * ae * ue - Phiw * aw * uw) / dxx
end subroutine
    
subroutine limiter_ke(L,Vp,Vw,Ve) 
    use variables, only: ONE
    use type_list, only: rancomp, ranke
    use interface_list, only: wdiv
    type(rancomp),intent(in) :: Vp,Vw,Ve
    type(ranke),intent(out)  :: L
    real*8 :: r
    !%%%%%%%%%%%!
    !     k     !
    !%%%%%%%%%%%!
    ! Calculate r taking into account possibility of different dx
    r = wdiv( (Vp%k - Vw%k) * (Ve%x - Vp%x) , (Ve%k - Vp%k) * (Vp%x - Vw%x) )
    ! Van Leer limiter
    L%k = wdiv( r + abs(r) , ONE + r )
    !%%%%%%%%%%%!
    !     e     !
    !%%%%%%%%%%%!
    ! Calculate r taking into account possibility of different dx
    r = wdiv( (Vp%e - Vw%e) * (Ve%x - Vp%x) , (Ve%e - Vp%e) * (Vp%x - Vw%x) )
    ! Van Leer limiter
    L%e = wdiv( r + abs(r) , ONE + r )
end subroutine    
!    
function center_diff_r(anu1,rc,rxp,rxm,ryp,rym,rzp,rzm,dt)
use type_list, only: rancomp, velocity, ranke
use interface_list, only: center_diff2
use variables, only: Se, Sk, Ser, ONE, HALF, ZERO, vt_op
implicit none
type(ranke) :: center_diff_r
real*8,intent(in) :: anu1, dt
real*8 :: kx, ky, kz, ex, ey, ez, Se1
type(rancomp),intent(in) :: rc,rxp,rxm,ryp,rym,rzp,rzm
type(velocity) :: Vp, Vm, rkc, rec

if (vt_op.eq.3) then
   Se1 = Ser
else
   Se1 = Se
endif

!Addition of aperture ratios included -> Oct 14 2015

rkc%v=rc%k;rkc%x=rc%x;
rec%v=rc%e;rec%x=rc%x;

vp%x=rxp%x;vm%x=rxm%x
#ifdef BUGWIN
if(vp%x.eq.rc%x.or.vm%x.eq.rc%x.or.dt.eq.ZERO) then
    write(6,*) 'wewerqwe0',vp%x,rc%x,vm%x,'dt=',dt
    continue
endif
#endif
vp%v=rxp%k;vm%v=rxm%k
kx=center_diff2(Vp,rkc,Vm,rxp%a*(anu1+(rc%d+rxp%d)/SK*HALF),  &
                          rxm%a*(anu1+(rc%d+rxm%d)/SK*HALF),dt)
vp%v=rxp%e;vm%v=rxm%e
ex=center_diff2(Vp,rec,Vm,rxp%a*(anu1+(rc%d+rxp%d)/Se1*HALF), &
                          rxm%a*(anu1+(rc%d+rxm%d)/Se1*HALF),dt)

vp%x=ryp%x;vm%x=rym%x
#ifdef BUGWIN
if(vp%x.eq.rc%x.or.vm%x.eq.rc%x) then
    write(6,*) 'wewerqwe1',vp%x,rc%x,vm%x
    continue
endif
#endif
vp%v=ryp%k;vm%v=rym%k
ky=center_diff2(Vp,rkc,Vm,ryp%a*(anu1+(rc%d+ryp%d)/SK*HALF),  &
                          rym%a*(anu1+(rc%d+rym%d)/SK*HALF),dt)
vp%v=ryp%e;vm%v=rym%e
ey=center_diff2(Vp,rec,Vm,ryp%a*(anu1+(rc%d+ryp%d)/Se1*HALF), &
                          rym%a*(anu1+(rc%d+rym%d)/Se1*HALF),dt)

vp%x=rzp%x;vm%x=rzm%x
#ifdef BUGWIN
if(vp%x.eq.rc%x.or.vm%x.eq.rc%x) then
    write(6,*) 'wewerqwe1',vp%x,rc%x,vm%x
    continue
endif
#endif
vp%v=rzp%k;vm%v=rzm%k
kz=center_diff2(Vp,rkc,Vm,rzp%a*(anu1+(rc%d+rzp%d)/SK*HALF),  &
                          rzm%a*(anu1+(rc%d+rzm%d)/SK*HALF),dt)
vp%v=rzp%e;vm%v=rzm%e
ez=center_diff2(Vp,rec,Vm,rzp%a*(anu1+(rc%d+rzp%d)/Se1*HALF), &
                          rzm%a*(anu1+(rc%d+rzm%d)/Se1*HALF),dt)

center_diff_r%k=kx+ky+kz            
center_diff_r%e=ex+ey+ez

endfunction

subroutine ke_chousei(nn1)
use arrays, only: r, r1
use variables, only: akini, eini, ieee, r_lim, HALF, cm, TENTH
implicit none
integer,intent(in) :: nn1
!
if (r1(nn1)%e.lt.eini.or.r1(nn1)%k.lt.akini.or.                               &
    r1(nn1)%e.gt.r_lim%e.or.r1(nn1)%k.gt.r_lim%k) then
    ieee = ieee + 1
    if (r1(nn1)%k.gt.r_lim%k.or.r(nn1)%k.gt.r_lim%k) then
        r1(nn1)%k = r_lim%k          
    elseif (r1(nn1)%k.lt.akini) then
        r1(nn1)%k = akini
    endif 
    if (r1(nn1)%e.gt.r_lim%e.or.r(nn1)%e.gt.r_lim%e) then
        r1(nn1)%e = r_lim%e     
    elseif (r1(nn1)%e.lt.eini) then
        r1(nn1)%e = eini
    endif 
endif   
end subroutine
    
subroutine turbulence_surface_boundary
    use arrays, only: nff, mn, in ,jn, kn, f, fln, r1, d, pro
    use variables, only: akrates, akini, eini, cm, ifln, ZERO
    implicit none
    real*8  :: ww, akww, eww, dww
    integer :: nn, nnf, i, k, j
    ! Finds previous air cells that have now become fluid/surface ones
    ! and seeds an initial estimate of the turbulence quantity
    ! based on the surrounding values for the next time step
!$omp parallel do default(none)                                                &
!$omp shared(ifln,fln,in,jn,kn,nff,f,r1,mn,d,pro,akini,eini)                   &
!$omp private(i,k,j,nn,ww,akww,eww,dww)
    FLUID_LOOP: do nnf = 1,ifln
      nn = fln(nnf)
      i = in(nn) ; j = jn(nn) ; k = kn(nn)
      if (nff(nn)%fp.eq.0) then
          pro(nn) = ZERO
          ww = f(mn(i-1,k,j)) + f(mn(i,k-1,j)) + f(mn(i,k,j-1))                &
             + f(mn(i+1,k,j)) + f(mn(i,k+1,j)) + f(mn(i,k,j+1))
          if (ww.le.ZERO) then
            if (ww.eq.ZERO) then
                r1(nn)%k = akini ; r1(nn)%e = eini
            else
                write(6,*) 'ww=',ww
                r1(nn)%k = max(akww/ww,akini) ; r1(nn)%e = max(eww/ww,eini)
            endif    
            d(nn) = cm * akini * akini / eini
          else
              akww = r1(mn(i-1,k,j))%k * f(mn(i-1,k,j))&
                   + r1(mn(i,k,j-1))%k * f(mn(i,k,j-1))&
                   + r1(mn(i,k-1,j))%k * f(mn(i,k-1,j))&
                   + r1(mn(i+1,k,j))%k * f(mn(i+1,k,j))&
                   + r1(mn(i,k,j+1))%k * f(mn(i,k,j+1))&
                   + r1(mn(i,k+1,j))%k * f(mn(i,k+1,j))
           
               eww = r1(mn(i-1,k,j))%e * f(mn(i-1,k,j))&
                   + r1(mn(i,k,j-1))%e * f(mn(i,k,j-1))&
                   + r1(mn(i,k-1,j))%e * f(mn(i,k-1,j))&
                   + r1(mn(i+1,k,j))%e * f(mn(i+1,k,j))&
                   + r1(mn(i,k,j+1))%e * f(mn(i,k,j+1))&
                   + r1(mn(i,k+1,j))%e * f(mn(i,k+1,j))
        
               dww = d(mn(i-1,k,j)) * f(mn(i-1,k,j))&
                   + d(mn(i,k,j-1)) * f(mn(i,k,j-1))&
                   + d(mn(i,k-1,j)) * f(mn(i,k-1,j))&
                   + d(mn(i+1,k,j)) * f(mn(i+1,k,j))&
                   + d(mn(i,k,j+1)) * f(mn(i,k,j+1))&
                   + d(mn(i,k+1,j)) * f(mn(i,k+1,j))
            
               r1(nn)%k = max( akww / ww, akini)
               r1(nn)%e = max( eww / ww,  eini )
               d(nn) = dww / ww
           endif
      endif
    enddo FLUID_LOOP
!$omp end parallel do
end subroutine

subroutine set_turbulence
    use variables, only: akrates, akrate, eini, cm, akini, akmin, emin, nu,   &
                         ike, inns, inne, g, ZERO, HALF, t, r_lim, ONE, inflow
    use arrays, only: r, r1, d, prop, pro
    implicit none
    integer :: nn

    akrates = ONE !0.8d0

    write(6,*) 'define AKZERO akrate=', akrate
    write(6,*) 'define akrates=', akrates
    if (t.eq.ZERO.or.ike.eq.0) then
        r(:)%k = akini ; r(:)%e = eini ; prop(:) = ZERO
    endif
    r1(inns:inne)%k = r(inns:inne)%k
    r1(inns:inne)%e = r(inns:inne)%e
    pro(inns:inne)  = prop(inns:inne)
    do nn = inns,inne
        if (r(nn)%e.ne.ZERO) then
            d(nn) = cm * r(nn)%k**2 / r(nn)%e
        else
            d(nn) = cm * akmin**2 / emin
        endif
    enddo
end subroutine set_turbulence
!    
