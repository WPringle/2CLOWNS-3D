!%%%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-2CLOWNS-local.F90 %%%%%%%%%%%%%%%%%%%%%%%%%
! The local module for generating and processing wave phenomena
subroutine read_local_data
    use variables
    use arrays, only: BC
    implicit none
    real*8       :: Chardepth
    integer      :: na1
    character*60 :: dummyc

    if (IM.eq.1.or.IM.eq.2) read(5,*) BC(2) ! 2D BC
    if (IM.eq.1.or.IM.eq.3) then
        read(5,*) BC(3)        ! 3D BC
        read(5,*) h_in         ! Is the elevation of the initial free-surface
                               ! relative to the arbitrary 3D - datum       
    endif
    call read_Wavedata
#ifdef RAN
    if (IM.eq.1.or.IM.eq.3) read(5,*) Chardepth, TI_b, nut_b
#endif
    read(5,'(A255)') output_file
#include "2CLOWNS_banner.inc"
#ifdef RAN
    if (IM.eq.1.or.IM.eq.3) then
    !--------------------------------------------------------------------------
    ! From: Lin, P & Liu, P L F, (1998).
    !       A numerical study of breaking waves in the surf zone
    akini = HALF * Chardepth * g * TI_b**2 
    eini  = cm * akini**2 / (nut_b * nu)  ! set coefficients prior to calc
    !                                     ! before being corrected here
    if (akini.lt.akmin) akini = akmin
    if (eini.lt.emin)   eini  = emin
    if (ce2i*eini/akini.gt.TEN.or.ce2i*eini/akini.lt.TENTH) then
        write(6,*) 'Initial epsilon and k ill-posed. e_i= ',eini,'k_i= ',akini
        stop
    endif
    ! Get upper limit based on the current cell average velocity 
    ! (i.e. turbulent intensity is limited to one)
    r_lim%k = HALF * Chardepth * g
    r_lim%e = sqrt(Chardepth * g**3)
    !--------------------------------------------------------------------------
    endif
#endif
end subroutine read_local_data
    
subroutine local_initial_process
    use variables, only: ZERO, LNUM, t, IM
    implicit none
    
    if (IM.eq.1.or.IM.eq.3) then
        !Get 3Dfree surface levels and initial water depths
        call calc_FSL
        call calc_IWD
        call calc_FLUX
        ! Change free surface level and velocity if specified
        if (t.eq.ZERO) then
            call calc_IC_or_GC('IC',1,LNUM+1)
            ! Change ground level if specified (just works like IC currently)
            call calc_IC_or_GC('GC',1,LNUM+1)
            call set_nflag
#ifdef PRINT_ZK
            call print_IWD ! Only use to get new IWD to make chikei again
            stop
#endif
        endif
    endif
    !initialise the output folders & DAT files
    call local_data_output(0)
end subroutine local_initial_process  
        
subroutine local_information  
    use variables, only: n, mwr, t, courant, velmax, dt, LNUM, nstart, IM
    use arrays, only: L
    implicit none
    integer            :: LN
    character*256 :: fm  = '(1h ," t = ",f10.4, " [s]")', &
    fm1 = '(1h ,"LN = ",I2," cran = ",f7.4," velmax = ",e12.5," dt = ",e12.5)'
    !
    if (mod(n,mwr).eq.1) then
        write(6,fm) t
        if (IM.eq.1.or.IM.eq.2) then
        do LN = 1, LNUM
            write(6,fm1) LN, L(LN)%CRAN, L(LN)%VMAX, L(LN)%DT
        enddo
        endif
        if ((IM.eq.1.and.nstart.ne.0).or.IM.eq.3) then
            write(6,fm1) LNUM + 1, courant, velmax, dt
        endif
    endif    
end subroutine local_information 


subroutine data_update_local
    call calc_FSL
    call calc_FLUX
end subroutine data_update_local
    
subroutine local_process

end subroutine local_process
    
function ubndS(ax,uc,dxc,dxb)
    use interface_list, only: ustar, wdiv
    use variables, only: karman, HALF, ks_r, ZERO, ONE, TWO
    implicit none
    real*8 :: ubndS, ust
    real*8,intent(in) :: ax, uc, dxc, dxb
    !
    ubndS = ZERO
    if (ax.eq.ZERO) ubndS = uc   ; return
    if (uc.eq.ZERO) ubndS = ZERO ; return
    if (ks_r(0).ge.ZERO) then
        !Wall function,: ks_r = 0: smooth wall, ks_r > 0: rough wall
        !Get the friction velocity from log-law
        ust = ustar( uc, HALF * ax * dxc,ks_r(0))
        !Now find the boundary condition value from the finite 
        !difference formula (see Lemos, C. (1992), A simple numerical
        !technique for turbulent flows with free surfaces)
        ubndS = uc - sign(ONE,uc) * ust * (dxc + dxb) / (karman * dxc)
    else
        ! Guess some small friction, e.g:
        ! Free slip -> ks_r = -1.0
        ! Half-slip -> ks_r = -2.0
        ! No-slip   -> ks_r = -3.0
        ubndS = uc * ( ks_r(0) + TWO ) 
    endif
endfunction
    
function ubndS2(ax,uc,dxc,dxb) !!!!@•¨‘Ì‹«ŠE‚Å‚Ì—¬‘¬Ý’è
    use interface_list, only: ustar, wdiv
    use variables, only: karman, HALF, ks_r, ZERO, ONE, TWO
    implicit none
    real*8 :: ubndS2, ust
    real*8,intent(in) :: ax, uc, dxc, dxb
    ubndS2 = ZERO
    if (ax.eq.ZERO) ubndS2 = uc   ; return
    if (uc.eq.ZERO) ubndS2 = ZERO ; return
    if (ks_r(10).ge.ZERO) then
        !Wall function,: ks_r = 0: smooth wall, ks_r > 0: rough wall
        !Get the friction velocity from log-law
        ust = ustar( uc, HALF * ax * dxc,ks_r(10))
        !Now find the boundary condition value from the finite 
        !difference formula (see Lemos, C. (1992), A simple numerical
        !technique for turbulent flows with free surfaces)
        ubndS2 = uc - sign(ONE,uc) * ust * (dxc + dxb) / (karman * dxc)
    else
        ! Guess some small friction, e.g.:
        ! Free slip -> ks_r = -1.0
        ! Half-slip -> ks_r = -2.0
        ! No-slip   -> ks_r = -3.0
        ubndS2 = uc * ( ks_r(10) + TWO )
    endif
end function
