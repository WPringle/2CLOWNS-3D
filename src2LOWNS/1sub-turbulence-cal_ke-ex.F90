subroutine cal_ke_ex
    use type_list, only: rancomp, ranke
    use interface_list, only: rbnd, center_diff_r, wdiv
    use variables, only: ifln, nu, akrates, dt, ce2i, ce2r, t, ieee, ce1, mu, &
                         vt_op, ZERO, ONE, HALF, ks_r, emin, akmin, TWO, TINY,&
                         SIXTH, FIVE, c1r
    use arrays, only: fln, fb, in, jn, kn, mn, un, nff, d, f, SS,             &
                      dx, pro, prop, r, r1, a, rhon, rho, lf, x
    implicit none
    real*8 :: akox1, eox1
    integer :: nn, im, ip, nnf, i, j, k, nfbnw, nnc
    logical :: blayer
    integer :: ipkj, ipkja, ippkj, immkj, l, ll, lll, lg(0:2)
    real*8 :: uc(0:2,-1:1), vv(0:2), dxc, nu1, ce2, r6, fm4, eta, C1
    real*8 :: akcf, ecf, rhocf, bfp, bfm, bfp1, bfm1
    type(rancomp) :: rc, rc0, rx(0:2,-1:1), wf(0:2,-1:1) 
    type(ranke)   :: fx(0:2), DIFF

    if (vt_op.eq.2) return !Smagorinsky 
#ifndef DAKU
    nu1 = nu
#endif

!$omp parallel do default(none)                                               &
!$omp shared(ifln,fln,f,fb,in,jn,kn,un,mn,r1,nff,r,d,dx,x,a,lf,ks_r           &
#ifndef DAKU
!$omp ,nu1                                                                    &
#endif
!$omp ,vt_op,akrates,dt,pro,prop,SS)                                          &
!$omp private(i,k,j,lg,ll,l,lll,nn,im,ip,nnc,uc,vv,rc,rc0,ce2,R6,FM4,dxc,wf   &
!$omp        ,blayer,rx,nfbnw,ipkj,ippkj,immkj,ipkja,fx,diff,akox1,eox1,C1,eta&
#ifdef DAKU
!$omp ,nu1&
#endif
#ifdef FURYOKU
!$omp ,rhocf,akcf,ecf,bfm,bfp,bfm1,bfp1&
#endif
!$omp ) !schedule(dynamic)
DO  nnf = 1,ifln
    nn = fln(nnf)
#ifdef BUGWIN
    if (nn.eq.365020) then
        continue
    endif
#endif
    if (f(nn)*fb(nn).eq.ZERO) cycle    
    ! Initialise k, e, x, nu, ce2
    rc%k  = r(nn)%k ; rc%e  = r(nn)%e ; rc%d  = d(nn) ; rc%x  = ZERO
    rc0%k = r(nn)%k ; rc0%e = r(nn)%e ; rc0%d = d(nn) ; rc0%x = ZERO
#ifdef DAKU
    nu1 = mu / rhon(nn)
#endif
    blayer = .false.

    if (vt_op.ne.1.or.rc%e.lt.TINY) then
        ! No adjustment for low-turbulence region
        if (vt_op.eq.3) then
            ce2 = ce2r  ! Realizable
        else
            ce2 = ce2i  ! Standard
        endif
    else
        ! Use the adjustment for low-turbulence regions
        R6 = -(SIXTH*rc%k**2/rc%e/nu1)**2
        if ( R6.LT.LOG(TINY) ) then
         FM4 = ONE
        else
         FM4 = ONE - TWO / 9.0d0 * EXP(R6)
        endif
        ce2 = ce2i * FM4
    endif

    i = in(nn); j = jn(nn); k = kn(nn)
    lg = [ i , j , k ]
    
    ! Get cell boundary and center velocities
    do l = 0,2
        uc(l,-1) = un(l,nn) ; uc(l,1) = un(l,mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))
        uc(l,0)  = ( uc(l,-1) + uc(l,1) ) * HALF
    enddo
    vv(0) = sqrt( sum(uc(1:2,0)**2) )
    vv(1) = sqrt( sum(uc(0:2:2,0)**2) )
    vv(2) = sqrt( sum(uc(0:1,0)**2) )
    wf%k = ZERO ; wf%e = ZERO 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!           Start loop over dimensions                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!    
    do l = 0,2
        dxc = dx(l,lg(l))
        ! Loop backwards and forwards
        do ll = -1,1,2
            ip = lg(l) + ll
            ipkja = mn(i+max(0,ll)*lf(l,0),                                   &
                       k+max(0,ll)*lf(l,2),                                   &
                       j+max(0,ll)*lf(l,1))
            ipkj  = mn(i+ll*lf(l,0),k+ll*lf(l,2),j+ll*lf(l,1))
            ! Add aperture ratio to rx Oct 14th 2015
            rx(l,ll)%a = a(l,ipkja)
            rx(l,ll)%x = ( dxc + dx(l,ip) ) * HALF
            if (nff(ipkj)%f.ge.1) then
                rx(l,ll)%k = r(ipkj)%k ; rx(l,ll)%e = r(ipkj)%e
                rx(l,ll)%d = d(ipkj)
            elseif (a(l,ipkja).eq.ZERO) then
                if (nff(ipkj)%b.eq.-1) then
                    rx(l,ll)%k = rc%k ; rx(l,ll)%e = rc%e
                else
                    !First get bc with no wall function
                    rx(l,ll)   = rbnd(rc,vv(l),dxc,-1d4)
                    rx(l,ll)%a = a(l,ipkja)
                    if (ks_r(nff(ipkj)%b).ge.ZERO) then
                        blayer = .true.
                        ! Now get bc with wall function
                        wf(l,ll) = rbnd(rc,vv(l),fb(nn)*dxc,ks_r(nff(ipkj)%b))
                    endif
                endif
            elseif (nff(ipkj)%f.eq.0) then
                ! Apply the no gradient boundary condition across 
                ! the free surface for k and epsilon
                nfbnw = nff(nn)%b ; lll = abs(nfbnw) - 1
                if (nfbnw.eq.-(l+1)*ll) then    
                    rx(l,ll)%k = rc%k * akrateS ; rx(l,ll)%e = rc%e
                    rx(l,ll)%d = rc%d * akrateS**2
                elseif (l.ne.lll) then
                    nnc = mn(i+ll*lf(l,0)+lf(lll,0),                          &
                             k+ll*lf(l,2)+lf(lll,0),                          &
                             j+ll*lf(l,1)+lf(lll,0))
                    if (nff(nnc)%b.eq.nfbnw) then
                        rx(l,ll)%k = r(nnc)%k * akrateS
                        rx(l,ll)%e = r(nnc)%e
                        rx(l,ll)%d = d(nnc) * akrateS**2            
                    else
                        rx(l,ll)%k = ( r(nnc)%k + rc%k ) * akrateS * HALF
                        rx(l,ll)%e = ( r(nnc)%e + rc%e ) * HALF
                        rx(l,ll)%d = ( d(nnc)   + rc%d ) * akrateS**2 * HALF 
                    endif
                endif
            elseif (nff(ipkj)%f.eq.-1.and.nff(ipkj)%b.gt.0) then
                rx(l,ll)%k = r(ipkj)%k ; rx(l,ll)%e = r(ipkj)%e
                rx(l,ll)%d = d(ipkj)
            else 
                continue
            endif
            rx(l,ll)%x = real(ll) * rx(l,ll)%x
        enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!           Estimate the advection terms of k and e                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        ! Add Second-order TVD - William April 8 2015
        ippkj = mn(i+2*lf(l,0),k+2*lf(l,2),j+2*lf(l,1))
        immkj = mn(i-2*lf(l,0),k-2*lf(l,2),j-2*lf(l,1))
        if (uc(l,1).gt.ZERO.and.uc(l,-1).gt.ZERO.and.nff(immkj)%f.ge.1) then
            rx(l,0)%k = r(immkj)%k ; rx(l,0)%e = r(immkj)%e
            rx(l,0)%x = HALF * ( x(l,lg(l)-1) + x(l,lg(l)-2)                  &
                                 - x(l,lg(l)) - x(l,lg(l)+1) )
            call TVD_ke(fx(l),rx(l,1)%a,rx(l,-1)%a,uc(l,1),uc(l,-1),          &
                        rx(l,1),rc,rx(l,-1),rx(l,0),dxc)
        elseif (uc(l,1).lt.ZERO.and.uc(l,-1).lt.ZERO.and.nff(ippkj)%f.ge.1) then
            rx(l,0)%k = r(ippkj)%k ; rx(l,0)%e = r(ippkj)%e
            rx(l,0)%x = HALF * ( x(l,lg(l)+1) + x(l,lg(l)+2)                  &
                                 - x(l,lg(l)) - x(l,lg(l)-1) )
            call TVD_ke(fx(l),rx(l,-1)%a,rx(l,1)%a,uc(l,-1),uc(l,1),          &
                        rx(l,-1),rc,rx(l,1),rx(l,0),-dxc)
        else    
            ! First-order upwind
            fx(l)%k = (rx(l,1)%a * ( max(uc(l,1),ZERO) * rc%k                 &
                                  + min(uc(l,1),ZERO) * rx(l,1)%k )           &
                    - rx(l,-1)%a * ( max(uc(l,-1),ZERO)* rx(l,-1)%k           &
                                  + min(uc(l,-1),ZERO)* rc%k) ) / dxc
            fx(l)%e = (rx(l,1)%a * ( max(uc(l,1),ZERO) * rc%e                 &
                                  + min(uc(l,1),ZERO) * rx(l,1)%e )           &
                    - rx(l,-1)%a * ( max(uc(l,-1),ZERO)* rx(l,-1)%e           &
                                  + min(uc(l,-1),ZERO)* rc%e) ) / dxc
        endif
    enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!           Estimate the viscous stress terms of k and e                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#ifdef FURYOKU
    rhocf=(rho(nn)+rho(mn(i,k-1,j)))*HALF
    if (rhocf.gt.600.d0) then
        akcf=(ak(nn)+ak(mn(i,k-1,j)))*HALF
        ecf=(e(nn)+e(mn(i,k-1,j)))*HALF
        Bfm=-g/rhocf*min( ZERO,rho(nn)-rho(mn(i,k-1,j)) ) &
           /(dx(2,k)+dx(2,k-1))/(ecf/akcf)**2
    else
        bfm=ZERO
    endif
        rhocf=(rho(nn)+rho(mn(i,k+1,j)))*HALF
    if(rhocf.gt.600.d0) then
        akcf=(ak(nn)+ak(mn(i,k+1,j)))*HALF
        ecf=(e(nn)+e(mn(i,k+1,j)))*HALF
        Bfp=-g/rhocf*min( ZERO,rho(mn(i,k+1,j))-rho(nn) ) &
            /(dx(2,k)+dx(2,k+1))/(ecf/akcf)**2
    else
        Bfp=ZERO
    endif
    bfm1 = ONE/(ONE+stpara*bfm);bfp1=ONE/(ONE+stpara*bfp);
    diff = center_diff_r_f(nu1,rc,rx(0,1),rx(0,-1),rx(1,1),rx(1,-1),rx(2,1),  &
                           rx(2,-1),bfp1,bfm1)
#else
    diff = center_diff_r(nu1,rc,rx(0,1),rx(0,-1),rx(1,1),rx(1,-1),rx(2,1),    &
                         rx(2,-1),dt) 
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!           Evaluate new k and e based on the transport equations             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Note addition of aperture ratios and void ratios in DIFF - Oct 14 2015
    ! Find the new turbulent dissipation rate
    if (rc%k.gt.ZERO) then
        if (vt_op.eq.3) then
           ! The realizable k-eps model equation for eps
           eta  = SS(nn) * rc%k / rc%e
           C1   = max( c1r , eta / ( eta + FIVE ) )
           eox1 = ( DIFF%e - sum(fx%e) ) / fb(nn) + C1 * SS(nn) * rc%e
           ! Semi-implicit method
           r1(nn)%e = ( rc0%e + DT * eox1 ) / ( ONE + DT * ce2 * rc%e /      &
                                              ( rc%k + sqrt(nu1 * rc%e) ) ) 
        else
           ! The standard k-eps model equations for eps
           eox1 = ( DIFF%e - sum(fx%e) ) / fb(nn) + ce1 * PRO(nn) * rc%e / rc%k 
           ! Find the new turbulent dissipation rate through semi implicit scheme
           r1(nn)%e = ( rc0%e + DT * eox1 ) / ( ONE + DT * ce2 * rc%e / rc%k ) 
        endif
    else
        r1(nn)%e = rc0%e
    endif
    ! Find the new turbulent kinetic energy
    ! Consider new form for PRO and E to be 
    ! half of the previous and current time step
    ! (Lemos, 1992. A simple numerical technique for turbulent flows 
    !  with free surfaces. Int. J. for Num. Methods in Fluids)
    akox1 = ( DIFF%k - sum(fx%k) ) / fb(nn)                                   &
          + HALF * ( PRO(nn) + PROP(nn) - r1(nn)%e - rc%e )   
          !akox1 =  ( DIFF%k - sum(fx%k) ) / fb(nn) + PRO(nn) - rc%e   
             !* f(nn) * fb(nn)min(1.0d0,PRO(nn))
    r1(nn)%k = rc0%k + DT * akox1
    ! Consider wall function 
    if (blayer) then
        ! Wall function only for full fluid cells
        !if (nff(nn)%f.eq.1) then
        r1(nn)%k = maxval(wf%k) ; r1(nn)%e = maxval(wf%e)
        !endif
    endif
    call ke_chousei(nn)
enddo
!$omp end parallel do
end subroutine cal_ke_ex

