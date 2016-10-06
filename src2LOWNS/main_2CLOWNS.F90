!%%%%%%%%%%%%%%%%%%%%%%% FILE: main_2CLOWNS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The main file where the calculation program is run from
#if defined(DRI) || defined(HENDO)
#define A_ a2
#define AY_ a2
#define AZ_ a2
#define FB_ fbd
#else
#define A_ a
#define AY_ af
#define AZ_ a
#define FB_ fb
#endif
program main
    call read_control_file
    call read_control_file2D
    call data_input
    call initial_setting
    call main_loop
    print *,'seijo shuryo'
    print *,'(End of Calculation)'
end program main
!    
subroutine initial_setting
    use arrays, only: lf, un, u, fp, f, a, alloc_quvwf, dim3d, fc, drift, nff,&
                      L
    use variables, only: t, t0, tw, tw0, dtw, dt, inne2, dt0, dtm, TWO, ONE,  &
                         is, ie, js, je, ko_out, g, ndri, ZERO, DWR, twl,     &
                         LNUM, HALF, nstart, IM
    implicit none
    integer :: LN
    !
    lf(0,0:2) = [1,0,0]; lf(1,0:2) = [0,1,0] ; lf(2,0:2) = [0,0,1]
    t = t0 ;
    twl = t + ONE/real(DWR)
    if (tw0.lt.t) then
        tw = t + dtw !最初のデータ出力時間の設定
    else
        tw = tw0     !最初のデータ出力時間の設定
    endif
    if (IM.eq.1.or.IM.eq.3) then
    ! Initialise 3D variables
    if (is.eq.ie-1.or.js.eq.je-1) then
        write(6,*) '2DV Simulation!' ; dim3d = .false.
        ko_out = 0
    else
        write(6,*) '3D Simulation!'  ; dim3d = .true.
    endif
    dtm = dt
    fc(0:2)   = [ZERO,ZERO,g] ! z方向負に重力
    call alloc_quvwf(inne2)
    call set_dx
    call set_nflag
    call init_matrix
#ifdef DAKU
    call set_turbidity
#endif
#ifdef RAN
    call set_turbulence
#endif
#if defined(DRI) || defined(NORMAL)
    call set_object
    if(drift) call set_drift

    call pre_surface
#endif
    call set_nflag
    endif
    ! Initialise 2D variables
    if (IM.eq.1.or.IM.eq.2) call init2D
    ! Initialise local variables
    call local_initial_process
    ! Initialise "previous" with "new" variables
    if (IM.eq.1.or.IM.eq.3) then 
        un = u; fp = f ; nff%fp = nff%f ; nff%bp = nff%b ; a(:,0) = ZERO
        if (drift) call drift_pressure   
    endif
    ! Set number of threads
!$  call set_ThreadsNumber
    if (t.ne.ZERO) return
    !
    if (IM.eq.1.or.IM.eq.3) call set_dt3d
    ! Adjust dt so that they are integer multiples of each other
    ! (unless just 3D calc.)
    if (IM.eq.1.or.IM.eq.2) call int_adj_dt(0)
    ! Initialise the pressure distribution (hydrostatic)
    if (IM.eq.1.or.IM.eq.3) call set_p
    !Find the water level at time n+1/2 (if u /= 0 it has effect)
    do LN = LNUM,1,-1 !Loop over the layers
        !Need to convert to half time-step dt 
        if (L(LN)%IC(2:2).ne.'T') L(LN)%DT = HALF * L(LN)%DT
        !-- Calculate continuity, depth, F and U
        call continuity(LN,0)
        !Reset dt back to original value
        if (L(LN)%IC(2:2).ne.'T') L(LN)%DT = TWO * L(LN)%DT
    enddo
    nstart = 1
    call data_output
    nstart = 0
end subroutine initial_setting
!    
subroutine data_update
    use variables, only: inns, inne
    use arrays, only: u, un, fp, f, nff, r, r1, c, cn,                        &
                      rho, rhon, sp, spn, prop, pro
    implicit none
    integer :: nn
!$omp parallel 
!$omp do
    do nn = inns,inne
        u(0:2,nn) = un(0:2,nn)
        fp(nn) = f(nn) 
        nff(nn)%fp = nff(nn)%f
        nff(nn)%bp = nff(nn)%b
#ifdef RAN
        r(nn)%k = r1(nn)%k
        r(nn)%e = r1(nn)%e
        prop(nn) = pro(nn)
#endif
#ifdef NORMAL
        sp(nn,0:2) = spn(nn,0:2)
#endif
    enddo
!$omp end do
#ifdef DAKU
!$omp workshare
    c = cn
    rho = rhon
!$omp end workshare
#endif
!$omp end parallel
    call data_update_local
!
end subroutine data_update
    
subroutine main_loop
    use arrays, only: L
    use variables, only: t, nkai, n, icong, ZERO, ONE, LNUM, nstart, MWR,      &
                         cr_ll, cr_sl, cr_ml, cr_bl, cr_ul, dt, IM, inflow, HALF
    implicit none
    integer :: LN3D = 0
    
    if (IM.eq.3) LN3D = 1
    do n = 1, nkai
       
        if (IM.eq.1) then
            ! Decide whether to start 3D calc. or not in hybrid
            if (inflow%wavetype.eq.'GT'.and.t - HALF*L(1)%dt.le.inflow%T_end) then
                ! Don't start 3D before transient rupture ends
            else
                if (nstart.eq.0) then
                    call calc3D_start(L(LNUM))
                    !When we decide to start 3D calc set LN3D to 1
                    if (nstart.ne.0) LN3D = 1
                    if (nstart.ne.1) then
                        call initial_3D
                        call set_nflag
                    endif
                endif
            endif
        endif
        ! Set DT forward
        if (IM.eq.3) then
            t = t + DT
        else
            t = t + L(1)%dt
        endif
        ! Call main recursive time loop
        call inner_loop_recursive(1,LN3D,1)
        
        if (icong >= 10000) goto 231
        ! Check 2D Courant
        if (IM.eq.1.or.IM.eq.2) then
            if (mod(n,mwr).eq.1) call check_courant
        endif
        ! Print out local info
        call local_information

        if (IM.eq.1.and.nstart.ne.0) then
            call check_timestep( cr_ll, cr_sl, cr_ml, cr_bl, cr_ul,           &
                                 ceiling(100d0 * DT / L(1)%DT) )
        
            call int_adj_dt(1) ! If dt has changed in check_timestep, alter dt
                               ! to be an integer multiple of the next layer  
        elseif (IM.eq.3) then
            call check_timestep( cr_ll, cr_sl, cr_ml, cr_bl, cr_ul, 100 )
        endif

        231 continue
    
        call data_output

        if (icong >= 1000) exit
 
    enddo
end subroutine main_loop 
    
recursive subroutine inner_loop_recursive(LN,LN3D,DT_RATIO)
    use arrays, only: drift, L
    use variables, only: icon, icong, dt, LNUM, ONE
    implicit none
    integer,intent(in) :: LN ,LN3D, DT_RATIO
    integer            :: LNN , DT_COUNT, DT_RATIO_UP
    real*8             :: DT_CF
    !
    if (icong >= 10000) return
!-----------------------------------------------------------------------------!
!   THIS SUBROUTINE LOOPS OVER THE INNER LAYERS IN A RECURSIVE MANNER         ! 
!          WHEN THE TIME STEPS BETWEEN EACH DOMAIN ARE DIFFERENT              !
!-----------------------------------------------------------------------------!
    Layer_Iter: do DT_COUNT = 1, max(DT_RATIO - 1,1)
        DT_CF = real(DT_COUNT) / real(DT_RATIO)
        Layer_loop_Up: do LNN = LN , LNUM + LN3D
        if (LNN.le.LNUM) then
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !                      2DH calculation                                !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !    
        !---------- Boundary condition in the outermost layer -----------------
            if (LNN.eq.1) call boundary_condition(2,ONE)
        !
        !---------- Calculate Msub, Nsub at n ---------------------------------
        !           Only for 2D Calc              
            if (L(LNN)%DIM.eq.2) call calcMNsub(L(LNN),LNN)
        !
        !---------- Momentum Flux Calc. at n+1 --------------------------------
            call leapfrog_flux(L(LNN),LNN)
        !
        !---------- Implicit Dispersion adjustment if specified ---------------
            call disp_adjust(L(LNN),LNN)
        !
        !---------- Linearly interpolate the fluxes at outer layer ------------
        !---------- along inner layer boundary; spatially and in time ---------
            if (LNN.gt.1) then
                if (LNN.eq.LN) call interp_outer_flux(L(LNN-1),L(LNN),DT_CF)  
                if (LNN.gt.LN) call interp_outer_flux(L(LNN-1),L(LNN),ONE)  
            endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !                       3D calculation                                !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        else
        !---------- Linearly interpolate the fluxes at outer layer ------------
        !---------- along inner 3D layer boundary; spatially and in time ------
#ifdef BUGWIN
            write(6,*) '3D-> DCOUNT=',DT_COUNT,'DRATIO=',DT_RATIO,DT
#endif
            if (LNN.eq.LN) call boundary_condition(3,DT_CF)
            if (LNN.gt.LN) call boundary_condition(3,ONE)
 
            call cal_momentum
            if (icon.ne.0) then 
                icong = 10000
                return
            endif

            call set_surface_velocity(2)

#ifdef DAKU
            call cal_rhon 
#endif
#ifdef RAN
            call cal_turbulence
#endif
            call cal_vof

#ifdef NORMAL
            call cal_surface_normal
#endif
            call surf_boundary_condition

            if (drift) call drift_pressure
        
            call local_process
    
            call cal_courant
            
            call data_update
        endif 
        !---------- If LNN is not in the innermost layer... -------------------
        if (LNN.lt.LNUM + LN3D) then
            !Find ratio of DT between outer and inner layer
            if (LNN.lt.LNUM) DT_RATIO_UP = int( L(LNN)%DT / L(LNN+1)%DT )
            if (LNN.eq.LNUM) DT_RATIO_UP = int( L(LNN)%DT / DT )
            if (DT_RATIO_UP.gt.1) then
            !- When the time steps are different between inner and ------------
            !- outer regions we need to calc continuity and flux again --------
            !- in a recursive manner ------------------------------------------
                call inner_loop_recursive(LNN+1,LN3D,DT_RATIO_UP)
            endif
        endif
        enddo Layer_loop_Up
        ! If in 2DH layer we need to calculate continuity
        if (LN.le.LNUM) then
            do LNN = LNUM,LN,-1 
                !-- Calculate continuity, depth, F and U
                call continuity(LNN,LN3D)
            enddo
        endif
    enddo Layer_Iter
end subroutine inner_loop_recursive     
!    
subroutine continuity(LN,LN3D)
    use variables, only: LNUM, HALF, ONE, t, inflow
    use arrays, only: L 
    implicit none
    integer,intent(in) :: LN, LN3D
    real*8             :: CF
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Calls the subroutines for layer interpolation, continuity, depth,       !
    !              velocity and F calculation in that order                   !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !
    !-- Interpolate the flux & depth from the inner layer to current layer ----
    if (L(LN)%coupling_dim.eq.2) then
        if (LN.lt.LNUM)                call interp_inner_flux(L(LN),L(LN+1),LN) 
        if (LN3D.eq.1.and.LN.eq.LNUM)  call interp_inner_flux(L(LN),L(LN),LN) 
    endif
!-- Update ZK or ETA for the transient tsunami condition
!-- This is slow and not parallelised but since transient condition is short consider OK
    if (inflow%wavetype(2:2).eq.'T'.and.t - HALF*L(LN)%DT.le.inflow%T_end) then
        if (L(LN)%IC.eq.'GT') then
            call calc_IC_or_GC('GT',1,LN)
            call Update_ETA_with_GC(L(LN),LN)  
        else
            call calc_IC_or_GC('IT',1,LN)
        endif
    endif
!-- Use continuity to calculate new water level at n+1/2 -----------------------
    ! Get the adjustment coefficient to ensure water level 
    ! is at correct time inside the coupling region
    if (LN.eq.LNUM) CF = HALF           ! Always 0.5 for 3D - 2D
    if (LN.lt.LNUM) CF = HALF * (ONE - L(LN+1)%DT / L(LN)%DT) 
                                        ! 0 <= CF < 0.5 for 2D - 2D 
    call leapfrog_cont(L(LN),LN,CF)
!
!-- Calc. depths at n+1/2 and interpolate with n-1/2 for n --------------------
    if (L(LN)%SWE.eq.1) call calc_depth(L(LN),LN)
!
!-- Calculating -> F, NF at n-1/2 and U, FLUX at n ----------------------------
    call calc_FUV(L(LN),LN)
!
endsubroutine
    
subroutine surf_boundary_condition

    call set_surface_velocity(2)
#ifdef RAN
    call turbulence_surface_boundary
#endif

end subroutine

subroutine data_input
      use variables, only: IM
      use interface_list, only: pritime2
      use arrays, only: date_time_start, date_time_old
      implicit none
    
      ! Read 2D data
      if (IM.eq.1.or.IM.eq.2) call read_data2D
      ! Read 3D data
      if (IM.eq.1.or.IM.eq.3) then
        call read_data
        call read_fdata
        call read_rdata
#ifdef NORMAL
        call read_nsdata
#endif

#ifdef DAKU
        call read_cdata
#endif
        call read_o2data
        call read_dridata
      endif
      call read_local_data
      
      call pritime2(date_time_start,date_time_old)

end subroutine data_input
!
subroutine data_output
    use interface_list, only: pritime2
    use arrays, only: date_time_start, date_time_old, drift, L
    use variables, only: t, t0, dt, tw, n, nkai, icong, t_end, ko_t, title,   &
                         nx1, nx2, ko_out, dtw, ko_tmp, testa, dwr, twl, ONE, &
                         ZERO, title2D, nstart, inflow, IM, HALF
    implicit none
    
    ! Output data at desired locations, formats etc
    ! with high frequency (every dwr timesteps)
    if (t-twl.gt.-1.0d-5) then
        call local_data_output(1)
        twl = twl + ONE/real(DWR)
    endif
    ! Output fdata at t = ZERO in some limited cases  
    if (t.eq.ZERO) then
        if (drift) then
            call set_title
            call output_fdata
        endif
        if (inflow%wavetype.eq.'GC') then
            call set_title2D
            call outputF2D(title2D)
        endif
    endif
    !
    if (t-tw.gt.-1.0d-5.or.n.eq.nkai.or.icong.gt.100.or.t.ge.t_end.or.t.eq.ZERO) then
        if (tw.ge.t0.or.ko_t.eq.-1.or.n.eq.nkai) then  !出力if
#ifdef NORMAL
            if (ko_out.eq.0.or.ko_t.eq.0.or.ko_t.eq.-1) then
#endif
                if (IM.eq.1.or.IM.eq.2) then
                    call set_title2D
                    call output2D('NOW')
                    if (inflow%wavetype.eq.'GT'.and.t - dtw < inflow%T_end) then
                        call outputF2D(title2D)
                    endif
                endif   
                if ((IM.eq.1.and.nstart.ne.0).or.IM.eq.3) then
                    call set_title
                    call owari(title,nx2) ; call hajime(title,nx1)
!$omp parallel sections
!$omp section
                    call output_data
#ifdef RAN
!$omp section
                    call output_rdata
#endif
#ifdef DAKU
!$omp section
                    call output_cdata_local
#endif
#ifdef NORMAL
!$omp section
                    call output_nsdata
!$omp section
                    call output_sfdata
#endif
!$omp section
                    if (drift) call output_drdata
!$omp end parallel sections
                endif
#ifdef NORMAL
            else
                if ((IM.eq.1.and.nstart.ne.0).or.IM.eq.3) then
                    call output_sfdata
                endif
            endif
#endif
            call pritime2(date_time_start,date_time_old)
        endif                                    !出力if
        if (tw.ge.t_end.or.icong.gt.100) then      
          if (IM.eq.1.or.IM.eq.2) call output2D('MAX')
          if (IM.eq.1.or.IM.eq.3) call write_max_data_MAT
          write(6,*) ' tw tend,icong=',tw, t_end, icong
          stop
        endif
        if (t-tw.gt.-1.0d-5) then
            tw = tw + dtw ; write(6,*) ' next tw=', tw
        endif
    endif
endsubroutine  

subroutine pritime2(date_time_start,date_time_old)
implicit none
integer::date_time_old(1:8),date_time_start(1:8)
integer::date_time(8)
real::s_int,s_start,tarray(2)
CHARACTER*8 idate
CHARACTER*10 ymd,ctime,zone
CHARACTER*255::fm="(1h ,'*** ',a4,'-',a2,'-',a2,' ',a2,':',a2,':',a2,&
', elapse=',f8.1,'[min] interval=',f8.2,'[s]')"

CALL DATE_AND_TIME(ymd,ctime,zone,date_time)
CALL dtime(tarray,s_int)
CALL etime(tarray,s_start)
s_start = s_start/60d0 ! Convert to minutes 

write(6,fm) ymd(1:4),ymd(5:6),ymd(7:8),ctime(1:2),ctime(3:4),ctime(5:6),s_start,s_int
   
date_time_old(1:8) = date_time(1:8)

end subroutine pritime2 
