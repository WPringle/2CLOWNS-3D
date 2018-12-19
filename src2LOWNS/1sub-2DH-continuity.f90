!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                       SUBROUTINE: LEAPFROG_CONT                         &&!
!&&        THIS SUBROUTINE CALCULATES THE FREE SURFACE LEVEL AND WATER      &&!
!&&         DEPTH FROM THE CONTINUITY EQN WITH 2ND ORDER ACCURACY           &&!
!&&                                                        - PRINGLE (2016) &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! 
subroutine leapfrog_cont(L,LN,DT_CF)
    use type_list, only: Layer
    use variables, only: LNUM, ZERO, HALF, VSMALL, nstart, inflow, t
    implicit none
    !----------------- Input/Output variables ---------------------------------
    real*8,intent(in)         :: DT_CF
    integer,intent(in)        :: LN
    type(Layer),intent(inout) :: L
    !----------------- Temporary variables used in the routine ----------------
    integer :: nnw, nn, i, j, ii, nnpx, nnpy, nnmx, loop, iter
    real*8  :: Qflux, Nsouth, Nnorth, c1
!
    loop = 0
141 iter = 0 ; loop = loop + 1
    if (loop.gt.10) then
        write(6,*) 'Iteration of depth correction not converging'
        do nn = 1,L%inne
            i = L%in(1,nn);j=L%in(2,nn)
            if (i.lt.L%is.or.i.ge.L%ie.or.j.lt.L%js.or.j.ge.L%je) cycle
            if (L%ETAn(nn).lt.L%ZK(nn)) then
                write(6,*) 'LN=',LN,'i=',i,'j=',j,                            &
                           'ETA=',L%ETAn(nn),'ZK=',L%ZK(nn)
            endif
        enddo
        return
    endif
!========== Calculate continuity in 1D, 2D ====================================
!
!$omp parallel private(i,j,ii,nn,nnpx,nnpy,nnmx,Qflux,Nsouth,Nnorth,c1)
!$omp do
!
CONT_LOOP: do nnw = 1,L%inwet
    !
    i = L%inw(1,nnw) ; j = L%inw(2,nnw) ; nn = L%mn(i,j)
    if (i.lt.L%is.or.i.ge.L%ie.or.j.lt.L%js.or.j.ge.L%je) cycle  
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%    If flux adjusted for negative depths..       %
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (loop.gt.1) then
        if (L%ETAn(nn).ge.L%ZK(nn)) then
            nnpx = L%mn(i+1,j) ; nnmx = L%mn(i-1,j)
            if (L%ETAn(nnmx).ge.L%ZK(nnmx).and.                               &
                L%ETAn(nnpx).ge.L%ZK(nnpx)) then
                if (L%DIM.eq.1) then
                    cycle !Cycle where flux is not adjusted
                elseif (L%DIM.eq.2) then
                    nnpx = L%mn(i,j+1) ; nnmx = L%mn(i,j-1)
                    if (L%ETAn(nnmx).ge.L%ZK(nnmx).and.                       &
                        L%ETAn(nnpx).ge.L%ZK(nnpx)) cycle
                          !Cycle where flux is not adjusted
                endif
            endif
        endif
    endif
!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%    Calculate the flux gradients                 %
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qflux = ZERO
    do ii = 1,L%DIM !Loop over the number of dimensions
        if (ii.eq.1) nnpx = L%mn(i+1,j)
        if (ii.eq.2) nnpx = L%mn(i,j+1)
        !Calculating one component of the flux derivative
        Qflux = Qflux + L%RX(ii,nn) * ( L%Qn(ii,nnpx) - L%Qn(ii,nn) )
    enddo
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% We can now calculate the new free surface level %
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !   
    ! Inside inner layer we need to extrapolate the 
    ! free-surface ahead using continuity by some time-fraction
    if (L%coupling_dim.eq.2.and.(nstart.ne.0.or.LN.ne.LNUM).and.              &
        i.ge.L%is2.and.i.lt.L%ie2.and.j.ge.L%js2.and.j.lt.L%je2) then
        L%ETAn(nn) = L%ETA(nn) - DT_CF * L%dt * Qflux
    else
        L%ETAn(nn) = L%ETA(nn) - L%dt * Qflux
    endif
!
!------------- Limit the value of the water surface level ---------------------
    if (abs(L%ETAn(nn)).lt.VSMALL)          L%ETAn(nn) = ZERO
    if (abs(L%ETAn(nn)-L%ZK(nn)).lt.VSMALL) L%ETAn(nn) = L%ZK(nn)
    if (L%ETAn(nn).gt.L%Z(L%KE))            L%ETAn(nn) = L%Z(L%KE)
!
enddo CONT_LOOP
!$omp end do
!
!Check for negative depths and adjust outgoing momentum fluxes by lambda ------
if (L%CD.eq.1) then
!$omp do
    CD_LOOP: do nnw = 1,L%inwet
        i = L%inw(1,nnw) ; j = L%inw(2,nnw) ; nn = L%mn(i,j)
        if (i.lt.L%is.or.i.ge.L%ie.or.j.lt.L%js.or.j.ge.L%je) cycle
        ! Dont worry about free surface within next layer
        if (L%coupling_dim.eq.2.and.(nstart.ne.0.or.LN.ne.LNUM).and.          &
            i.ge.L%is2.and.i.lt.L%ie2.and.j.ge.L%js2.and.j.lt.L%je2) cycle
        ! Cycle at free-surfaces greater than the ground level
        if (L%ETAn(nn).ge.L%ZK(nn)) cycle
        nnpx = L%mn(i+1,j)
        if (L%DIM.eq.1) then
            Nsouth = ZERO ; Nnorth = ZERO
        elseif (L%DIM.eq.2) then
            nnpy   = L%mn(i,j+1)
            Nsouth = L%Qn(2,nn) ; Nnorth = L%Qn(2,nnpy)
        endif
        call adjust_flux(Nsouth,Nnorth,L%Qn(1,nn),L%Qn(1,nnpx),               &
                         L%ETAn(nn)-L%ZK(nn),L%ETA(nn)-L%ZK(nn),L%DX(1,nn),   &
                         L%DX(L%DIM,nn),L%Dmin,L%dt,loop) 
        if (L%DIM.eq.2) then
            L%Qn(2,nn) = Nsouth; L%Qn(2,nnpy) = Nnorth
        endif    
        if (iter.ne.1) iter=1
    enddo CD_LOOP
!$omp end do
endif
!$omp end parallel
if (iter.eq.1) goto 141
if (LN.gt.1) call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%ETAn,1,1)
!
endsubroutine leapfrog_cont
!
subroutine adjust_flux(Ns,Nn,Mw,Me,Dc,Dcm,DX,DY,Dmin,dt,loop) 
    use variables, only: SMALL, ZERO, ONE
    use interface_list, only: wdiv
    implicit none
    !----------------- Input/Output variables ---------------------------------
    integer,intent(in) :: loop
    real*8,intent(in) :: Dcm,DX,DY,Dmin,dt !(Opposite for negative dt)
    real*8,intent(inout) :: Nn,Ns,Me,Mw
    !----------------- Temporary variables used in the routine ----------------
    integer :: iter
    real*8 :: Qin, Qout, lambda, Dcn, Dc
    !=========================================================================!
    !  Adjusts the outflow fluxes so that the new depth remains non-negative  !
    !=========================================================================!
    !
    Qin  = max(Ns,ZERO)      / DY + max(Mw,ZERO)      / DX +                  &
           abs(min(Nn,ZERO)) / DY + abs(min(Me,ZERO)) / DX
    Qout = abs(min(Ns,ZERO)) / DY + abs(min(Mw,ZERO)) / DX +                  &
           max(Nn,ZERO)      / DY + max(Me,ZERO)      / DX
!
    if (dt.gt.ZERO) then !Normal calculation
        lambda = wdiv( Dcm + Qin*dt , Qout*dt )
    elseif (dt.lt.ZERO) then !Switch Qin to Qout when dt is neg.
        lambda = wdiv( Dcm + Qout*dt , Qin*dt )
    endif
    iter = 0
    do ! Infinite loop reiterating lambda until depth is zero or positive
        if (Ns*dt.lt.ZERO) Ns = Ns * lambda
        if (Mw*dt.lt.ZERO) Mw = Mw * lambda
        if (Nn*dt.gt.ZERO) Nn = Nn * lambda
        if (Me*dt.gt.ZERO) Me = Me * lambda
        ! Return if depth is >= 0
        Dcn = Dcm - dt * ( Nn - Ns ) / DY - dt * ( Me - Mw ) / DX
        if ( Dcn >= ZERO ) return
        ! Relax requirement
        lambda = max(ZERO,min(ONE - Dcn / Dc,ONE - SMALL))
        iter = iter + 1
        Dc = Dcn
        if (iter.gt.10000) then
            write(6,*) lambda,Dcn,Dc,Dcm,Ns,Nn,Mw,Me
            stop 'Many iterations in adjust_flux (>10000)'
        endif
    enddo
!
endsubroutine adjust_flux
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                            SUBROUTINE: CALC_DEPTH                       &&!
!&&         THIS SUBROUTINE CALCULATES THE WATER DEPTHS AT THE CELL         &&!
!&&         BOUNDARIES USING LINEAR INTERPOLATION FOR NON-UNIFORM           &&!
!&&                      CELL SIZES AND TIME STEPPING.                      &&!
!&&                                                        - PRINGLE (2016) &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! 
subroutine calc_depth(L,LN)
    use type_list, only: Layer
    use arrays, only: ix
    implicit none
    !----------------- Input/Output variables ---------------------------------
    integer,intent(in)        :: LN
    type(Layer),intent(inout) :: L
    !----------------- Temporary variables used in the routine ----------------
    integer :: nn, nnw, i, j, ii
    !
    !==========================================================================
    !----------------- Calculate cell boundary depths -------------------------
!$omp parallel do private(nn,i,j,ii)
DEPTH_LOOP: do nnw = 1,L%inwet 
    !
    i = L%inw(1,nnw) ; j = L%inw(2,nnw) ; nn = L%mn(i,j)
    if (i.lt.L%is.or.i.gt.L%ie.or.j.lt.L%js.or.j.gt.L%je) cycle
!
    do ii = 1,L%DIM !Loop over the number of dimensions
        ! Cycle at domain edges 
        if ((ii.eq.1.and.j.eq.L%je).or.(ii.eq.2.and.i.eq.L%ie)) cycle
        ! Get the necessary depths for the cell boundary
        call get_depths(L,LN,nn,i,j,ii)
        ! Get depth for adjacent boundary or for a shoreline cell boundary 
        if (L%mn(i+ix(ii,1),j+ix(ii,2)).eq.0) cycle
        if ((L%NF(i+ix(ii,1),L%ks,j+ix(ii,2)).eq.-1.and.                      &
             L%NFB(i+ix(ii,1),L%ks,j+ix(ii,2)).ne.0).or.                      &
            (L%NF(i,L%ks,j).eq.0.and.L%NFB(i,L%ks,j).ne.0)) then
            call get_depths(L,LN,L%mn(i+ix(ii,1),j+ix(ii,2)),                 &
                                      i+ix(ii,1),j+ix(ii,2),ii)
        endif
    enddo
!
enddo DEPTH_LOOP
!$omp end parallel do
call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%D,L%DIM,2)
call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%Dn,L%DIM,2)
!
endsubroutine calc_depth  
!    
subroutine get_depths(L,LN,nn2,i2,j2,ii2)
    use type_list, only: Layer
    use variables, only: ZERO, LNUM, nstart, HALF
    implicit none
    !----------------- Input/Output variables ---------------------------------
    integer,intent(in) :: LN, nn2, i2, j2, ii2
    type(Layer),intent(inout) :: L
    !----------------- Temporary variables used in the routine ----------------
    integer :: nnmx,nnpx,nnmmx
    real*8 :: find_depth
    !------------ Determine the depths at the cell boundaries ----------------- 
    if (ii2.eq.1) then
        nnmx = L%mn(i2-1,j2) ; nnpx = L%mn(i2+1,j2) ; nnmmx = L%mn(i2-2,j2) 
    elseif (ii2.eq.2) then
        nnmx = L%mn(i2,j2-1) ; nnpx = L%mn(i2,j2+1) ; nnmmx = L%mn(i2,j2-2) 
    endif
    if (L%SWE.eq.1) then
        !The depth at the n+1/2 time step
        L%Dn(ii2,nn2) = find_depth(L%ETAn(nnmx),L%ETAn(nn2),L%H(ii2,nn2),     &
                                   L%ZK(nnmx),L%ZK(nn2),L%DX(ii2,nnmx),       &
                                   L%DX(ii2,nn2),L%Qn(ii2,nn2),L%Dmin)
    elseif (L%SWE.eq.0.and.L%H(ii2,nn2).gt.ZERO) then 
        !Calculate depth only at initial stage
        L%Dn(ii2,nn2) = L%H(ii2,nn2) ; L%D(ii2,nn2) = L%H(ii2,nn2) ; return
    endif   
    ! Within next layer depths for nonlinear terms
    ! are directly obtained from that layer
    if (L%coupling_dim.eq.2.and.(nstart.ne.0.or.LN.ne.LNUM).and.              &
        i2.ge.L%is2.and.i2.lt.L%ie2.and.j2.ge.L%js2.and.j2.lt.L%je2) then 
        ! Cycle inside of the next layer domain
        if (ii2.eq.1.and.i2.gt.L%is2) return
        if (ii2.eq.2.and.j2.gt.L%js2) return
        !if (ii2.eq.1.and.i2.gt.L%is2+1.and.i2.lt.L%ie2-1) return
        !if (ii2.eq.2.and.j2.gt.L%js2+1.and.j2.lt.L%je2-1) return
    endif
    ! The depth at the n time step interpolating in time and space
    ! - for nonlinear & bed friction terms, and for calculating velocity
    L%D(ii2,nn2) = find_depth(HALF * (L%ETAn(nnmx) + L%ETA(nnmx)),            &
                              HALF * (L%ETAn(nn2) + L%ETA(nn2)),              & 
                              L%H(ii2,nn2),L%ZK(nnmx),L%ZK(nn2),L%DX(ii2,nnmx)&
                             ,L%DX(ii2,nn2),L%Qn(ii2,nn2),L%Dmin)
endsubroutine get_depths  
!
function find_depth(etam,etap,h,zkm,zkp,dxm,dxp,q,Dmin)
    use variables, only: ZERO, ONE, HALF
    use interface_list, only: LI_1
    implicit none 
    real*8,intent(in)  :: etam, etap, h, zkm, zkp, dxm, dxp, Dmin, q
    real*8 :: find_depth
    !=========================================================================!
    !  Finds the water depth at the cell boundary based on various inputs     !
    !=========================================================================!
    !
    if (-h.gt.zkm.and.-h.gt.zkp) then
        !<<<<<Considering Banks, Breakwaters or Seawalls<<<<<
        ! Upwind method
        if (q.eq.ZERO) then  
            find_depth = max(etap,etam) + h
        elseif (q.gt.ZERO) then
            find_depth = etam + h
        elseif (q.lt.ZERO) then
            find_depth = etap + h
        endif
    else !<<<<When no Banks, Breakwaters or Seawalls are present<<<<<
        if ( q.eq.ZERO.or.                                                    &
            ! Inludes effect when fluid is disconnected for overtopping
            (-h.eq.zkp.and.etam.lt.zkp).or.(-h.eq.zkm.and.etap.lt.zkm) ) then
            find_depth = max(etap,etam) - max(zkm,zkp)
        else
            find_depth = h + LI_1(etam,etap,dxm,dxp)
        endif
    endif
    if (find_depth.lt.Dmin) find_depth = ZERO
    !
endfunction find_depth            
