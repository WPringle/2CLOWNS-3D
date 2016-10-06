!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                       SUBROUTINE: LEAPFROG_FLUX                         &&!
!&&           THIS SUBROUTINE CALCULATES THE MOMENTUM FLUXES USING          &&!
!&&              THE LEAPFROG METHOD VALID ON NON-UNIFORM GRIDS             &&!
!%%                                                      - PRINGLE (2016)   &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! 
subroutine leapfrog_flux(L,LN)
    use type_list, only: Layer
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    integer,intent(in) :: LN
    !----------------- Temporary variables used in the routine ----------------
    integer :: i, j, nn, ii, nnw
    !
    !========== Main Loop of Calculation: =====================================
    ! - Fully parallelized with OpenMP
    ! - Functions are used to calculate derivates and 
    !   are found at the bottom of this file
    ! - If a wall is present at the cell boundary and the water level is above 
    !   the wall then an empirical overtopping formula (Homma's) is used for 
    !   the flux calculation instead of the usual momentum equations
    ! - Friction is considered using Manning's n formula and is calculated 
    !   semi-implicitly.
    ! - We can calculate either the linear or nonlinear shallow water
    !   equations in cartesian or spherical coordinates
    ! - If desired, for the linear shallow water equations we can optimise the 
    !   numerical dispersion to mimic the frequency dispersion of the linear 
    !   boussinesq equations (set DISP = 1). 
    !==========================================================================
!$omp parallel do private(nn,i,j,ii) 
!
    MOM_LOOP: do nnw = 1, L%inwet
    !
        i = L%inw(1,nnw) ; j = L%inw(2,nnw) ; nn = L%mn(i,j)
        if (i.lt.L%is.or.i.ge.L%ie.or.j.lt.L%js.or.j.ge.L%je) cycle
    !
        DIM_LOOP: do ii = 1, L%DIM !Loop over the number of dimensions 
            ! Get the necessary momentum for the cell boundary
            call get_mom(L,LN,nn,i,j,ii)
        enddo DIM_LOOP
    !
    enddo MOM_LOOP
!$omp end parallel do  
endsubroutine leapfrog_flux
!    
subroutine get_mom(L,LN,nn2,i2,j2,ii2)
    use type_list, only: Layer 
    use interface_list, only: UPW_COR 
    use arrays, only: ix
    use variables, only: ZERO, HALF, ONE, TWO, TW4, TWOTHIRD, TWFTH,          &
                         COMFLOW, SUBFLOW, R, TW7, g, LNUM, nstart
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    integer,intent(in) :: LN, nn2, i2, j2, ii2
    !----------------- Temporary variables used in the routine ----------------
    integer :: jj, nnpx, nnmx, nnmmx, nnpy, nnppy, nnmy, nnmmy
    integer :: nnpymx, nnmxpy, nnppx, nnmxmy, nnmypx
    real*8  :: hxxx, hxyy, h1, h2, MU, MV, GU, FX, CDCUN_4_2
    !
    ! Dont need to calc fluxes within next layer
    if (L%coupling_dim.eq.2.and.(nstart.ne.0.or.LN.ne.LNUM)) then
        if (i2.ge.L%is2.and.i2.lt.L%ie2.and.j2.ge.L%js2.and.j2.lt.L%je2) then 
            ! Cycle inside of the next layer domain
            if (L%mn2(ii2,i2,j2).le.0) return
        endif
    endif

    if (ii2.eq.1) then
        if (i2.eq.L%is) return !Outermost boundaries (calculated in BC)                         
        nnmx = L%mn(i2-1,j2) ; nnpx = L%mn(i2+1,j2)
    elseif (ii2.eq.2) then
        if (j2.eq.L%js) return !Outermost boundaries (calculated in BC)         
        nnmx = L%mn(i2,j2-1) ; nnpx = L%mn(i2,j2+1)
    endif
    if (L%Dn(ii2,nn2).eq.ZERO) then
        L%Qn(ii2,nn2) = ZERO ; return !Set momentum to zero for zero depth
    endif
!<<<<<<<<<<<<<<<<<<<<<<<<< If not cycled above... <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!================= Calculate Momentum in the ii2-direction ====================
!
    if ( (-L%H(ii2,nn2).gt.L%ZK(nn2).and.-L%H(ii2,nn2).gt.L%ZK(nnmx)).or.     &
        ! New conditions added for overtopping of larger breakwater
        (-L%H(ii2,nn2).eq.L%ZK(nn2).and.L%ETAn(nnmx).lt.L%ZK(nn2)).or.        &
        (-L%H(ii2,nnmx).eq.L%ZK(nnmx).and.L%ETAn(nn2).lt.L%ZK(nnmx)) ) then
!Special Case with wall
!<<<<<<<<<<< If a seawall is present at the cell boundary <<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<< use empirical weir formula to calculate flux <<<<<<<<<<<<<<<<<<<<<
        h1 = L%ETAn(nnmx) + L%H(ii2,nn2) ; h2 = L%ETAn(nn2) + L%H(ii2,nn2)
        if (h1.gt.h2) then                    
            if (h2.gt.TWOTHIRD * h1) then
                !Submerged overflow Eqn.
                L%Qn(ii2,nn2) = SUBFLOW * h2 * sqrt( TWO * g * ( h1 - h2 ) )
            else
                !Complete overflow Eqn.
                L%Qn(ii2,nn2) = COMFLOW * h1 * sqrt( TWO * g * h1 )
            endif
        elseif (h2.gt.h1) then
            if (h1.gt.TWOTHIRD*h2) then
                !Submerged overflow Eqn.
                L%Qn(ii2,nn2) = - SUBFLOW * h1 * sqrt( TWO * g * ( h2 - h1 ) )
            else
                !Complete overflow Eqn.
                L%Qn(ii2,nn2) = - COMFLOW * h2 * sqrt( TWO * g * h2 )
            endif
        else
            L%Qn(ii2,nn2) = ZERO
        endif              
        return !return here <<<<
    endif
!Normal Case 
!<<<<<<<<<<<< Calculate flux from momentum equation <<<<<<<<<<<<<<<<<<<<<<<<<<<
!------- Get the additional adjacent cell numbers -----------------------------
    if (L%SWE.eq.1.or.L%DISP.eq.1) then 
        jj = L%DIM - ii2 + 1
        if (ii2.eq.1) then
            nnppy = L%mn(i2,j2+2) ; nnmmx  = L%mn(i2-2,j2)
            nnpy  = L%mn(i2,j2+1) ; nnmy   = L%mn(i2,j2-1)
            nnmmy = L%mn(i2,j2-2) ; nnpymx = L%mn(i2-1,j2+1) 
            if (i2.ne.L%ie) nnppx = L%mn(i2+2,j2) ; if (i2.eq.L%ie) nnppx = 0
        elseif (ii2.eq.2) then
            nnppy = L%mn(i2+2,j2) ; nnmmx  = L%mn(i2,j2-2)
            nnpy  = L%mn(i2+1,j2) ; nnmy   = L%mn(i2-1,j2)
            nnmmy = L%mn(i2-2,j2) ; nnpymx = L%mn(i2+1,j2-1) 
            if (j2.ne.L%je) nnppx = L%mn(i2,j2+2) ; if (j2.eq.L%je) nnppx = 0
        endif
    endif
!
!------- Calculate the pressure gradient term ---------------------------------
!
    if (L%DISP.eq.2.and.                                                      &
        L%ETAn(nnpx).gt.L%ZK(nnpx).and.L%ETAn(nnmmx).gt.L%ZK(nnmmx)) then
        ! Fourth-order difference
        GU = g * L%Dn(ii2,nn2) * TW4 * ( TW7 * (L%ETAn(nn2) - L%ETAn(nnmx))   &
                              + L%ETAn(nnmmx) - L%ETAn(nnpx) ) * L%RXc(ii2,nn2)    
    else
        ! Second-order difference
        GU = g * L%Dn(ii2,nn2) * (L%ETAn(nn2) - L%ETAn(nnmx)) * L%RXc(ii2,nn2)
    endif
!  
!------ Adjust the pressure gradient term to include explicit disperson -------
!       correction if desired for linear shallow water equations --------------
    if (L%DISP.eq.1) then
        ! Ignore near the boundaries
        if ( i2.gt.L%is.and.i2.lt.L%ie.and.                                   &
            ( L%Dim.eq.1.or.(j2.gt.L%js.and.j2.lt.L%je) ) ) then
            ! Get additional cell numbers
            if (ii2.eq.1) then
                nnmxpy = L%mn(i2-1,j2+1) ; nnmxmy = L%mn(i2-1,j2-1)
                nnmypx = L%mn(i2+1,j2-1)
            elseif (ii2.eq.2) then
                nnmxpy = L%mn(i2+1,j2-1) ; nnmxmy = L%mn(i2-1,j2-1)
                nnmypx = L%mn(i2-1,j2+1)
            endif
            ! Calculate the third-order derivatives of the free surface -------
            call CDSUN_4_2(hxxx,hxyy,L%ETAn(nnmmx),L%ETAn(nnmx),L%ETAn(nn2),  &
                    L%ETAn(nnpx),L%ETAn(nnmxpy),L%ETAn(nnpy),L%ETAn(nnmxmy),  &
                    L%ETAn(nnmy),L%RX(ii2,nn2),L%DIM)
            ! Ignore dispersion effects near the shoreline & lower bound
            if (L%Dn(ii2,nnmx).le.ZERO.or.L%Dn(ii2,nnpx).le.ZERO) hxxx = ZERO
            if (L%Dn(ii2,nnmy).le.ZERO.or.L%Dn(ii2,nnpy).le.ZERO) hxyy = ZERO
            ! Update the pressure gradient term
            GU = GU + g * L%H(ii2,nn2)                                        &
               * TWFTH * ( L%APA(ii2,nn2) * hxxx + L%GMA(ii2,nn2) * hxyy )                    
        endif
    endif
!
!--- Calculate the Coriolis forcing term and add to pressure gradient ---------
    if (L%GEOC.eq.1) then
        if (ii2.eq.1) then
            GU = GU - L%Qs(jj,nn2) * HALF * ( L%sinY(j2) + L%sinY(j2+1) )
        elseif (ii2.eq.2) then
            GU = GU + L%Qs(jj,nn2) * L%sinY(j2)
        endif
    endif
!
!---------- Calculating the nonlinear advection terms   -----------------------                                    
    ADVECTION: if (L%SWE.eq.1) then       
!
!------- Advection term in the ii2 x ii2 direction ----------------------------
        ! Upwind method with correction for truncation error terms
        MU = UPW_COR(L%Q(ii2,nnmx) ,L%Q(ii2,nnmx),L%D(ii2,nnmx),              &
                     L%Q(ii2,nn2)  ,L%Q(ii2,nn2) ,L%D(ii2,nn2),               &
                     L%Q(ii2,nnpx) ,L%Q(ii2,nnpx),L%D(ii2,nnpx),              &
                     L%D(ii2,nnmx) ,L%D(ii2,nnpx),L%D(ii2,nnmx),L%D(ii2,nnpx),&
                     L%RX(ii2,nnmx),L%RX(ii2,nn2),L%DISP,L%H(ii2,nn2),g,L%DT)  
!------------ Advection cross term in the ii2 x jj direction ------------------
        if (L%DIM.eq.2) then !2D Calc only
            ! Upwind method with correction for truncation error terms
            MV = UPW_COR(L%Qs(jj,nnmy) ,L%Q(ii2,nnmy),L%D(ii2,nnmy),          &
                         L%Qs(jj,nn2)  ,L%Q(ii2,nn2) ,L%D(ii2,nn2),           &
                         L%Qs(jj,nnpy) ,L%Q(ii2,nnpy),L%D(ii2,nnpy),          &
                         L%D(jj,nn2)   ,L%D(jj,nnpy) ,L%D(jj,nnmx),           & 
                         L%D(jj,nnpymx),L%RX(jj,nnmy),L%RX(jj,nn2),           &
                         L%DISP,L%H(ii2,nn2),g,L%DT) 
        else
            MV = ZERO
        endif
!
!       !>>>>> Final Calculation for SWE Equations>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!
!=========== Calculating the new f term in the ii2-direction for SWE ==========
        L%Qn(ii2,nn2) = GU + MU + MV  
             
    else !>>>>> Final Calculation for LSW Equations>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!
!=========== Calculating the new f term in the ii2-direction for LSW ==========
        L%Qn(ii2,nn2) = GU
             
    endif ADVECTION
!     
endsubroutine get_mom
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                    SUBROUTINE: DISP_ADJUST                              &&!
!&&         THIS SUBROUTINE ADJUSTS THE MOMENTUM FLUXES CALCULATED IN       &&!
!&&       LEAPFROG_FLUX FOR DISPERSION CAUSED BY THE NON-HYDROSTATIC        &&!
!&&        PRESSURE VIA POISSON PRESSURE CONSIDERATIONS  - PRINGLE (2016)   &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
! Reference:                                                                  !
! Shigihara, Fujima (2007). Adequate numerical scheme for dispersive          !
! wave theory for tsunami simulation and development of  new numerical        ! 
! algorithm. JSCE Coastal and Environmental Engineering, 63 (1), 51-66        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine disp_adjust(L,LN)
    use type_list, only: Layer
    use arrays, only: ix
    use variables, only: ZERO, eps0, HALF, A, B, g
    use interface_list, only: CDCUN_2_1
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    integer,intent(in)        :: LN
    !----------------- Temporary variables used in the routine ----------------
    integer :: yes, i, j, ii, nn, nnpx, nnmx, nnw, ITER = 1000
    real*8,dimension(1:L%indisp) :: PSI
    real*8 :: MXT, hxx, DISP, FX, Utot, SLOPE, SLOPEm, BI
    real*8 :: Break_Condition
    real*8,parameter :: BE = HALF, PARM = 5d0
!
if (L%DISP.eq.2) then ! Only if dispersion required
    YES = 0
!============== FIRST CALCULATE THE RHS OF THE MATRIX EQUATION ================
!$omp parallel do private(i,j,ii,nn,nnpx,nnmx,hxx,MXT,SLOPE,SLOPEm) 
!
    PSI_LOOP: do nnw = 1,L%indisp
        !
        i = L%IND(1,nnw) ; j = L%IND(2,nnw) ; nn = L%MN(i,j)
        !
        !================ Calculate all the values of RHS  ====================
        PSI(nnw) = ZERO !Set PSI to zero as default value
        !Cycle in dry regions
        if (L%ETAn(nn)-L%ZK(nn).lt.L%Dmin) cycle                 
        hxx = ZERO ; MXT = ZERO ; SLOPE = L%ETAn(nn) / -L%ZK(nn)
        do ii = 1,L%DIM !Loop over the number of dimensions
            if (ii.eq.1) then
                ! Cycle for first couple of cells for stability
                if (i.eq.L%is.or.i.eq.L%ie-1) cycle
                nnpx = L%mn(i+1,j); nnmx = L%mn(i-1,j)
            elseif (ii.eq.2) then
                ! Cycle for first couple of cells for stability
                if (j.eq.L%js.or.j.eq.L%je-1) cycle
                nnpx = L%mn(i,j+1); nnmx = L%mn(i,j-1)
            endif
            if (L%BC == 1) then
                ! Ignore disp terms where local momentum gradients are high
                ! (Roeber & Cheung (2012) 
                SLOPEm = max(L%Q(ii,nnpx) - L%Q(1,nn),ZERO )                  &
                       * L%RX(ii,nn) / sqrt( g*( L%ETAn(nn)-L%ZK(nn) ) )
                SLOPE = L%ETAn(nn) / -L%ZK(nn)
                ! Ignore disp terms where eta/h gradients are high 
                ! (Tonelli and Petti, 2009)
                ! Updated to include criteria on slopes 
                ! where flat bed limit is not valid
                BI = Break_Condition(L%H(ii,nnmx),L%H(ii,nn),L%DX(ii,nn),     &
                                     HALF * ( L%Q(ii,nnpx) + L%Q(ii,nn) ) )
                if (i.gt.L%is.and.i.lt.L%ie-1.and.                            &
                    (L%DIM.eq.1.or.(j.gt.L%js.and.j.lt.L%je-1)).and.          &
                     SLOPE.ge.BI.and.L%BREAK(nnw).eq.0) then
                    ! Turn off disp terms when gradients and depths are high
                    write(6,*) 'Breaking initiated at i =', i, 'SLOPE=', SLOPE
                    L%BREAK(nnw) = 1
                elseif(SLOPEm.lt.BE.and.SLOPE.lt.BI.and.L%BREAK(nnw).eq.1) then
                    ! Turn on disp terms when gradients and depths are small
                    write(6,*) 'Breaking ended at i =', i, 'SLOPE=', SLOPE
                    L%BREAK(nnw) = 0
                endif
                if (L%BREAK(nnw).eq.1) then
                    hxx = ZERO ; MXT = ZERO ; exit
                endif
            endif
            !When considering the nonlinear shallow water equations..
            if (L%SWE.eq.1) then
                !Considering already calculated guess value from SWE..
                MXT = MXT + L%RX(ii,nn) * ( L%Qn(ii,nnpx) - L%Qn(ii,nn) )
            endif
            if (nnpx.eq.0) nnpx = nn ; if (nnmx.eq.0) nnmx = nn 
            if (L%ZK(nnmx).gt.ZERO.or.L%ETAn(nnmx).eq.L%ZK(nnmx)) nnmx = nn
            if (L%ZK(nnpx).gt.ZERO.or.L%ETAn(nnpx).eq.L%ZK(nnpx)) nnpx = nn
            !Second derivate of ETA with ii^2
            hxx = hxx + CDCUN_2_1(L%ETAn(nnmx),L%RXc(ii,nn),L%ETAn(nn),       & 
                                  L%RXc(ii,nnpx),L%ETAn(nnpx),2)
        enddo
        !
        !-------------------- Evaluate RHS ------------------------------------
        if (L%SWE.eq.0) then
            !<<<< Calculation for LSW <<<<<<<<<<<<<<<
            PSI(nnw) = (A - B) * g * L%ZK(nn) * L%ZK(nn) * hxx
        elseif (L%SWE.eq.1) then
            !<<<< Calculation for SWE <<<<<<<<<<<<<<<
            PSI(nnw) = A * -L%ZK(nn) * MXT - B * g * L%ZK(nn) * L%ZK(nn) * hxx
        endif
        if (YES.eq.0) then
            if (PSI(nnw).ne.ZERO) YES = 1 
        endif
    !
    enddo PSI_LOOP
!$omp end parallel do
!
!========= Call the Bi-CGSTAB MATRIX solver for PSI ============================                    
    ! We input PSI as the RHS of the matrix equation 
    ! and the solution vector PSI is given on return 
    if (YES == 1) then
        call PCGPME_WILL(L%AL,4,L%indisp,L%LL,L%DIAG,L%indisp,                 &
                         64,PSI,eps0,PARM,ITER) 
    endif
!
endif
!========= ADJUST FLUX USING THE SPATIAL DERIVATIVE OF PSI =====================
!
!$omp parallel do private(i,j,ii,nn) 
!
ADJ_LOOP: do nnw = 1,L%inwet 
    !
    i = L%inw(1,nnw) ; j = L%inw(2,nnw) ; nn = L%mn(i,j)
    if (i.lt.L%is.or.i.ge.L%ie.or.j.lt.L%js.or.j.ge.L%je) cycle
    ! Find the local momentum gradients
    do ii = 1,L%DIM !Loop over the number of dimensions
        ! Get the necessary momentum for the cell boundary
        call mom_final(L,LN,PSI,nn,i,j,ii)   
    enddo
!
    enddo ADJ_LOOP
!$omp end parallel do 
!=========== END OF SUBROUTINE ================================================
!
end subroutine disp_adjust
    
subroutine mom_final(L,LN,PSI,nn2,i2,j2,ii2)
    use variables, only: VSMALL, SVNTHIRD, ONE, EIGHTH, HALF, ZERO,           &
                         g, nstart, LNUM 
    use type_list, only: Layer
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    real*8,dimension(1:L%indisp),intent(in) :: PSI
    integer,intent(in) :: LN, nn2, i2,j2, ii2
    !----------------- Temporary variables used in the routine ----------------
    integer :: jj,nnd,nndm,nnmx
    real*8 :: DISP,FX,Utot
    !
    ! Dont need to calc fluxes within next layer
    if (L%coupling_dim.eq.2.and.(nstart.ne.0.or.LN.ne.LNUM)) then        
        if (i2.ge.L%is2.and.i2.lt.L%ie2.and.j2.ge.L%js2.and.j2.lt.L%je2) then 
            ! Return inside of the next layer domain
            if (L%mn2(ii2,i2,j2).le.0) return   
        endif
    endif
    if (L%Dn(ii2,nn2).eq.ZERO) return ! return for zero depths
    if (ii2.eq.1) then
        if (i2.eq.L%is) return !return at lower boundary
        nnmx = L%mn(i2-1,j2) ; jj = L%DIM 
    elseif (ii2.eq.2) then
        if (j2.eq.L%js) return !return at lower boundary
        nnmx = L%mn(i2,j2-1) ; jj = 1
    endif
    ! Return for flux calculated overtopping weir (already calculated)
    if ( (-L%H(ii2,nn2).gt.L%ZK(nn2).and.-L%H(ii2,nn2).gt.L%ZK(nnmx)).or.     &
        ! New conditions added for overtopping of larger breakwater
        (-L%H(ii2,nn2).eq.L%ZK(nn2).and.L%ETAn(nnmx).lt.L%ZK(nn2)).or.        &
        (-L%H(ii2,nnmx).eq.L%ZK(nnmx).and.L%ETAn(nn2).lt.L%ZK(nnmx)) ) return 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !  Calculating the bed-friction term 
    if (L%MAN(nn2).ne.ZERO.and.L%MAN(nnmx).ne.ZERO                            &
        .and.L%D(ii2,nn2).ne.ZERO) then
        if (L%DIM.eq.1) Utot = abs(L%Q(ii2,nn2))
        if (L%DIM.eq.2) Utot = sqrt(L%Q(ii2,nn2)**2 + L%Qs(jj,nn2)**2)
        FX = EIGHTH * g * ( ( L%MAN(nn2) + L%MAN(nnmx) )**2 ) * Utot          &
           / L%D(ii2,nn2)**(SVNTHIRD)
    else
        FX = ZERO
    endif
    ! Now update the new fluxes
    if (L%DISP.eq.2) then
        nnd = L%mnd(i2,j2)
        if (ii2.eq.1) nndm = L%mnd(i2-1,j2)
        if (ii2.eq.2) nndm = L%mnd(i2,j2-1) 
        if ( nnd.eq.0.or.nndm.eq.0.or.( LN.eq.1.and.                          &
            ( (ii2.eq.1.and.(i2.eq.L%is+1.or.i2.eq.L%ie-1) ).or.              &
            (ii2.eq.2.and.(j2.eq.L%js+1.or.j2.eq.L%je-1))) ) ) then
            ! Cycle on the first cell inside the domain
            DISP = ZERO
        else
            ! Get the contributuon due to dispersion
            DISP = L%H(ii2,nn2) * L%RXc(ii2,nn2) * (PSI(nnd) - PSI(nndm))  
        endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !   Get new flux SWE + DISP  !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L%Qn(ii2,nn2) = ( L%Q(ii2,nn2) * ( ONE - L%dt * FX )                  &
                       + L%dt * ( DISP - L%Qn(ii2,nn2) ) ) / (ONE + L%dt * FX )
    else
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !   Get new flux SWE         !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L%Qn(ii2,nn2) = ( L%Q(ii2,nn2) * ( ONE - L%dt * FX )                  &
                        - L%dt * L%Qn(ii2,nn2) ) / (ONE + L%dt * FX )
    endif
    ! Limit the minimum allowable discharge 
    if (abs(L%Qn(ii2,nn2)).lt.VSMALL) L%Qn(ii2,nn2) = ZERO
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endsubroutine mom_final
!    
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                        SUBROUTINE: CALCMNSUB                            &&!
!&&            THIS SUBROUTINE INTERPOLATES THE FLUXES USING FOUR           &&!
!&&       SURROUNDING NODES (ONLY FOR 2D CALCULATION)     - PRINGLE (2016)  &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
subroutine calcMNsub(L,LN)
    use type_list, only: Layer
    use arrays, only: ix
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    integer,intent(in)        :: LN
    !----------------- Temporary variables used in the routine ----------------
    integer :: nn, i, j, ii, nnw
    !---------- Set fluxes equal at boundaries --------------------------------
    if (LN.eq.1) call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%Q,L%DIM,2)
!$omp parallel do private(i,j,ii,nn)
!
SUB_LOOP: do nnw = 1,L%inwet 
    !
    i = L%inw(1,nnw) ; j = L%inw(2,nnw) ; nn = L%mn(i,j)
    if (i.lt.L%is.or.i.ge.L%ie.or.j.lt.L%js.or.j.ge.L%je) cycle
    ! Find the local momentum gradients
    do ii = 1,L%DIM !Loop over the number of dimensions
        ! Get the 4-point average of the momentum
        call av_4_point(L,nn,i,j,ii)
        ! Get 4-point average of momentum at boundaries
        if (L%NF(i+ix(ii,1),L%ks,j+ix(ii,2)).eq.-1.and.                       &
            L%NFB(i+ix(ii,1),L%ks,j+ix(ii,2)).ne.0.and.                       &
            L%NFB(i+ix(ii,1),L%ks,j+ix(ii,2)).ne.200) then
            call av_4_point(L,L%mn(i+ix(ii,1),j+ix(ii,2)),                    &
                            i+ix(ii,1),j+ix(ii,2),ii)
        endif    
    enddo
!
enddo SUB_LOOP
!$omp end parallel do
call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%Qs,L%DIM,2)
!
endsubroutine calcMNsub
!    
subroutine av_4_point(L,nn2,i2,j2,ii2)
    use type_list, only: Layer
    use variables, only: HALF
    use interface_list, only: LI_1
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    integer,intent(in) :: nn2, i2, j2, ii2
    !----------------- Temporary variables used in the routine ----------------
    integer :: jj, nnpxmy, nnpx, nnmy
    
    if (ii2.eq.1) then
        jj = 2 ; nnpx = L%mn(i2+1,j2)
        nnmy = L%mn(i2,j2-1) ; nnpxmy = L%mn(i2+1,j2-1)
    elseif (ii2.eq.2) then
        jj = 1 ; nnpx = L%mn(i2,j2+1)
        nnmy = L%mn(i2-1,j2) ; nnpxmy = L%mn(i2-1,j2+1) 
    endif
!
!------------ Interpolate the surrounding fluxes to get Qsub ------------------
    if ( (ii2.eq.1.and.j2.ne.L%js-1).or.(ii2.eq.2.and.i2.ne.L%is-1) ) then
        if ( (ii2.eq.1.and.i2.eq.L%ie).or.(ii2.eq.2.and.j2.eq.L%je) ) then
            ! Two point interpolation on the upper boundary
            L%Qs(ii2,nn2) = LI_1(L%Q(ii2,nnmy),L%Q(ii2,nn2),                  &
                                 L%DX(jj,nnmy),L%DX(jj,nn2))   
        elseif ((ii2.eq.1.and.i2.eq.L%is-1).or.                               &
                (ii2.eq.2.and.j2.eq.L%js-1)) then
            ! Two point interpolation on the lower boundary
            L%Qs(ii2,nn2) = LI_1(L%Q(ii2,nnpxmy),L%Q(ii2,nnpx),               &
                                 L%DX(jj,nnmy),L%DX(jj,nn2))        
        else
            !Four point interpolation inside calculation domain
            L%Qs(ii2,nn2) = LI_1( HALF*(L%Q(ii2,nnmy) + L%Q(ii2,nnpxmy)),     &
                                  HALF*(L%Q(ii2,nn2) + L%Q(ii2,nnpx)),        &
                                  L%DX(jj,nnmy),L%DX(jj,nn2) ) 
        endif
    endif
endsubroutine av_4_point
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&         FUNCTION: UPW_COR -> 1st-order Upwind with correction           &&!
!&&            On a Staggered Uniform/Non-Uniform Grid in 1-D.              &&! 
!&&             This function returns values of u.fx from two               &&!
!&&                 surrounding nodes on a structured grid                  &&!   
!&&                                                        - PRINGLE (2016) &&!  
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!    
function UPW_COR(xm,ym,dm,x,y,d,xp,yp,dp,dym,dyp,                             &
                 dymxm,dypxm,rxm,rx,DISP,h,g,dt)
    use variables, only: ZERO, HALF, ONE
    use interface_list, only: wdiv
    implicit none
    real*8             :: UPW_COR
    integer,intent(in) :: DISP 
    real*8,intent(in)  :: xm, ym, x, y, xp, yp, dm, d, dp, dym, dyp,          &
                          dymxm, dypxm, rxm, rx, h, g, dt
    real*8             :: Q, Qm, Qp, Cr
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !   First-order upwind method with adjustment for truncation error terms  !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    if (d.eq.ZERO) then
        UPW_COR = ZERO; return
    endif
    Q  = wdiv(x  * y  , d )
    if (x.ge.ZERO) then !>>>>>>>>>> Positive flow >>>>>>>>>>>>>>>>>>>>>>>>>>>>!
        if (dm.eq.ZERO) then
            UPW_COR = ZERO; return
        endif
        Qm = wdiv(xm * ym , dm)    
        if (x.ne.y.or.ym.ne.xm) then
            ! For cross terms: 
            ! Cannot transfer flow across boundary where there is a wall etc.
            if (dym.eq.ZERO.or.dymxm.eq.ZERO) then
                UPW_COR = ZERO ; return
            endif
        endif    
        if (DISP.eq.2.and.h.gt.ZERO) then
            ! For implicit scheme,
            ! Calculation with back-substitution of truncation errors
            Qp = wdiv(xp * yp , dp) 
            if (dm.eq.ZERO.or.d.eq.ZERO.or.dp.eq.ZERO.or.                     &
               (dp-d)*(d-dm).lt.ZERO) then ! Eliminating unsmooth regions
                Cr = ZERO
            else
                Cr = HALF * (ONE - sqrt(g * h) * rxm * dt)
            endif
            UPW_COR = ( Q - Qm + Cr * ( Qp - Q + Qm - Q ) ) * rxm 
        else
            ! For normal explicit scheme
            UPW_COR = ( Q - Qm ) * rxm 
        endif
    else                !>>>>>>>>>> Negative flow >>>>>>>>>>>>>>>>>>>>>>>>>>>>!
        if (dp.eq.ZERO) then
            UPW_COR = ZERO; return
        endif
        Qp = wdiv(xp * yp , dp) 
        if (x.ne.y.or.yp.ne.xp) then
            ! For cross terms: 
            ! Cannot transfer flow across boundary where there is a wall etc.
            if (dyp.eq.ZERO.or.dypxm.eq.ZERO) then
                UPW_COR = ZERO ; return
            endif
        endif  
        if (DISP.eq.2.and.h.gt.ZERO) then
            ! For implicit scheme,
            ! Calculation with back-substitution of truncation errors
            Qm = wdiv(xm * ym , dm)    
            if (dp.eq.ZERO.or.d.eq.ZERO.or.dm.eq.ZERO.or.                     &
               (dp-d)*(d-dm).lt.ZERO) then ! Eliminating unsmooth regions
                Cr = ZERO
            else
                Cr = HALF * (ONE - sqrt(g * h) * rx * dt)
            endif
            UPW_COR = ( Qp - Q - Cr * ( Qp - Q + Qm - Q ) ) * rx
        else
            ! For normal explicit scheme
            UPW_COR = ( Qp - Q ) * rx
        endif
    endif
    !
    endfunction UPW_COR
!    
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&            FUNCTION:  CDCUN_2_1 -> Central Discretization               &&!
!&&   on a Collated Uniform/Non-Uniform Grid to 2nd-order accuracy in 1-D.  &&! 
!&&    This function returns values of fx (iflag=1) or fxx (iflag=2) from   &&!
!&&     three surrounding nodes on a structured grid for uniform or         &&!
!&&      non-uniform values of dx                         - PRINGLE (2016)  &&!                  
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!  
real*8 function CDCUN_2_1(fm,rxm,f,rxp,fp,iflag)
    use variables, only: TWO, HALF, ZERO, ONE
    use interface_list, only: LI_1
    implicit none
    real*8,intent(in) :: fm, fp, rxm, rxp
    integer,intent(in) :: iflag
    real*8 :: f, dxm, dxp
    !
    if (rxp.eq.rxm) then
        !!!Uniform grid size!!!
        if (iflag.eq.1) then     !Solve for fx
            CDCUN_2_1 = ( fp - fm ) * rxm
        elseif (iflag.eq.2) then !Solve for fxx
            CDCUN_2_1 = ( fp - f + fm - f ) * rxm * rxm
        else
            write(6,*) 'Cannot calculate derivate of order:',iflag
            return
        endif
    else
        !!!Non-uniform grid size!!!
        dxm = ONE / rxm ; dxp = ONE / rxp
        if (iflag.eq.1) then     !Solve for fx
            f = LI_1(fm,fp,dxm,dxp)
            dxm = HALF * dxm ; dxp = HALF * dxp
            CDCUN_2_1 = ( fp * dxm**2 - fm * dxp**2 + f * ( dxp**2 - dxm**2 ))&
                                   / ( dxm * ( dxp * dxp + dxm * dxp ) )
        elseif (iflag.eq.2) then !Solve for fxx
            CDCUN_2_1 = TWO * (fp * dxm - f * ( dxp + dxm ) + fm * dxp )      &
                                   / ( dxm * ( dxp * dxp + dxm * dxp ) )
        else
            write(6,*) 'Cannot calculate derivate of order:',iflag
            return
        endif
    endif
    endfunction CDCUN_2_1
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                         SUBROUTINE: CDSUN_4_2                           &&!
!&&            Central Discretization on a Staggered Uniform Grid to        &&!
!%%                         4th-order accuracy in 2-D.                      &&! 
!&&        This subroutine returns values of fxxx and fxyy for uniform      &&!
!&&            or non-uniform values of dx & dy            - PRINGLE (2016) &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! 
subroutine CDSUN_4_2(fxxx,fxyy,xmm,xm,xp,xpp,xmyp,xpyp,xmym,xpym,rx,dim)
    use variables, only: HALF, ZERO, THREE
    implicit none
    real*8,intent(in)  :: xmm, xm, xp, xpp, xmyp, xpyp, xmym, xpym, rx
    integer,intent(in) :: dim
    real*8,intent(out) :: fxxx, fxyy
    real*8  :: CDCUN_2_1, c, d, dx, dx1, dx2
    integer :: info
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !                              SOLVE FOR FXXX                             !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    fxxx = ( xpp - xmm + THREE * ( xm - xp ) ) * rx 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !                              SOLVE FOR FXYY                             !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !Set to zero and return if one-dimension calculation..
    if (DIM.eq.1) then
        fxyy = ZERO
        return
    endif
    fxyy = (xpyp - xp + xpym - xp + xm - xmyp + xm - xmym) * rx
    !
    endsubroutine CDSUN_4_2 
!    
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                      FUNCTION:  Break_Condition                         &&!
!&&      Finds the condition for breaking based on the local slope          &&!
!&&                                                        - PRINGLE (2016) &&!  
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! 
real*8 function Break_Condition(hm,hp,dx,q)
    use variables, only: ZERO, ONE
    implicit none
    real*8,intent(in) :: hm, hp, q, dx
    real*8 :: S
    !-------------------------------------------------------------------------!
    !Reference: Grilli, et al. (1997). Breaking Criterion and                 !
    !           Characteristics for Solitary Waves on Slopes                  !
    !-------------------------------------------------------------------------!
    ! --> For rundown the criteria will become the normal one on zero slope
    S = max(ZERO,(hm - hp) * sign(ONE,q) / dx)
    !
    Break_Condition = 0.75 + 25d0 * S - 112d0 * S * S + 3870d0 * S * S * S
endfunction Break_Condition