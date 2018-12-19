!%%%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-wave-BCIC.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   THIS FILE WAS CREATED TO FACILITATE INITIAL CONDITIONS SUCH AS            !
!    A TSUNAMI SOURCE (SUDDEN CHANGE IN GROUND AND/OR WATER LEVEL,            !
!    AUTOMATING BOUNDARY CONDITIONS GENERALLY ASSOCIATED WITH WAVE LIKE FLOW  !
!           ALONG WITH FACILITATING THE HYBRID COUPLING ALGORITHM             !    
!                                            - BY WILLIAM PRINGLE DEC (2015)  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine read_Wavedata
    use variables, only: g, PI_8, ONE, TWO, ZERO, inflow, inflow_file, GDR
    use interface_list, only: calc_sw_per, calc_ns_per, calc_k
    implicit none
    integer :: nn, i, unit, j, na, ist
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      READ_WAVEDATA: Reads the wave type and properties                      !
!                                                      WILLIAM PRINGLE (2015) !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    !Read the wave type and location of the input
    read(5,'(A2,1x)',advance='no') inflow%wavetype
    !Read the information of the input
    if (inflow%wavetype.eq.'SW') then
        ! Solitary wave
        read(5,*) inflow%h_wave, inflow%h_shoki, inflow%angle
        ! Calculate the wave period from the eqn.
        inflow%period = calc_sw_per(inflow%h_shoki,inflow%h_wave,g)
    elseif (inflow%wavetype.eq.'RW') then
        ! Regular wave
        read(5,*) inflow%h_wave, inflow%h_shoki, inflow%period, inflow%angle  
        ! Calculate the wave number 
        inflow%k = calc_k(inflow%period,inflow%h_shoki,g)  
    elseif (inflow%wavetype(1:1).eq.'A'.or.inflow%wavetype(1:1).eq.'M') then
        !Prescribed data (arbitrary wave)           
        read(5,*) inflow%width, inflow%dt, inflow%angle, inflow%t_adj
    elseif (inflow%wavetype(2:2).eq.'C'.or.inflow%wavetype(2:2).eq.'T') then
        !Initial (water level, velocity or ground) condition
        read(5,'(A2,1x)',advance='no') inflow%order
        if (inflow%order.eq.'EQ') then
            read(5,*) inflow%h_wave, inflow%h_shoki, inflow%x0, inflow%angle
            inflow%v_wave = ZERO
        elseif (inflow%wavetype(2:2).eq.'C') then
            read(5,*) inflow%fmt1
            inflow%h_wave = ZERO
        elseif (inflow%wavetype(2:2).eq.'T') then
            read(5,*) inflow%fmt1, inflow%dt, inflow%T_end
            inflow%h_wave = ZERO
        endif
    elseif (inflow%wavetype.eq.'UF') then
        read(5,*) inflow%h_wave, inflow%h_shoki, inflow%angle
    else
        write(6,*) inflow%wavetype 
        stop 'Invalid wave type (as printed above): Please correct'
    endif
    !Read the inflow_file 
    if (inflow%wavetype.eq.'SW'.or.                                           &
       (inflow%wavetype(2:2).eq.'C'.and.inflow%order.eq.'EQ')) return
    read(5,'(A255)') inflow_file
    if (inflow_file(1:3).eq.'NON') return
    na = len_trim(inflow_file)
    if (inflow%wavetype(1:1).eq.'A'.or.inflow%wavetype(1:1).eq.'M') then
        open(newunit = unit, file = inflow_file, status = 'old')
            !Find number of data points
            do while (.true.)
                read(unit,*,END=482) 
                inflow%size = inflow%size + 1
            enddo 
            482 rewind(unit)
            allocate(inflow%hwave(inflow%size))
            if (inflow%width.gt.1) then
                allocate(inflow%vwave(min(2,inflow%width-1),inflow%size))
                if (inflow%width.eq.4) allocate(inflow%mwave(inflow%size))
            endif
            do i = 1,inflow%size
                if (inflow%width.eq.1) then
                    ! Input AW shape at boundary using LSW 
                    ! to estimate momentum at boundary
                    read(unit,*) inflow%hwave(i)
                elseif (inflow%width.lt.4) then
                    ! Input AW shape at boundary prescribing velocity
                    ! (and hence momentum)
                    read(unit,*) inflow%hwave(i),                             &
                                 ( inflow%vwave(j,i),j = 1, inflow%width-1 )
                elseif (inflow%width.eq.4) then
                    ! Input AW shape at boundary prescribing velocity
                    ! (hence momentum) and movement of the boundary (wavemaker)
                    read(unit,*) inflow%hwave(i),                             &
                                 ( inflow%vwave(j,i), j = 1,2), inflow%mwave(i)
                endif
            enddo
            if (inflow%wavetype(2:2).ne.'E') then
                ! Calculate the main wave period (via FFT)
                inflow%period = calc_ns_per(inflow%hwave,inflow%size,inflow%dt) 
                if (inflow%width.gt.1) inflow%v_wave = maxval(inflow%vwave)
                if (inflow%width.eq.4) inflow%m_wave = maxval(inflow%mwave)   &
                                                     - inflow%mwave(1)
            else
                inflow%v_wave = TWO * TWO * sqrt(PI_8) / inflow%period
            endif
            !Get the max heights and velocities, wavemaker movement
            inflow%h_wave = maxval(inflow%hwave(:))
        close(unit)
    endif       
    !Read Ground Deformation Rate
    read(5,*,iostat=ist) GDR
    if (ist.ne.0) then
        backspace(5) ! if GDR is not set in inputfile.
        GDR = ONE
    endif
!
endsubroutine read_Wavedata
!    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%     Subroutine that gives the initial water profile and velocities       %%!
!%            Based on some given user specified theory                     %%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine init_height_vel(eta,vel,dim,nnn,ii,jj,H,d,x0,y0,direc)
    use variables, only: HCIRC, VSMALL, TWOTHIRD, ZERO, ONE, TWO,             &
                         HALF, THREE, g, THREEQ 
    use arrays, only: L
    implicit none
    !----------------- Temporary variables ------------------------------------
    real*8,intent(in)  :: H, d, y0, direc
    integer,intent(in) :: ii, jj, nnn, dim
    real*8,intent(out) :: eta, vel(dim)
    real*8 :: XX, C, kn, xp, yp, xm, ym, x0, etam
    real*8 :: delta, c2, ksi, a, a1, a2, b
!
    ! Use Boussinesq's first approximation for k
    kn = sqrt(THREEQ * H / d**3 )
    !Solitary wave
    if (x0.eq.-1d4) then
        x0 = TWO * d * atanh(0.995d0) / sqrt(THREE * H/d)
        !Calculate from Synolakis
        !x0 = - (1d0/kn) * acosh(sqrt(20d0))
    endif
    xm = L(1)%X(ii) - x0 
    if (kn*abs(xm).lt.acosh(sqrt(H/VSMALL))) then !Limiting height to above VSMALL
        if (dim.eq.2) vel(2) = ZERO
        xp = L(1)%X(ii+1) - x0 ; xx = HALF * (xp + xm)
        yp = L(1)%Y(jj+1) - y0 ; ym = L(1)%Y(jj) - y0
        if (L(1)%DISP.eq.0) then
            ! Use Boussinesq's first approximation
            C      = sqrt(g * ( H + d ) )
            etam   = H * ( cosh(kn * xm) )**-2
            eta    = H * ( cosh(kn * HALF*(xp+xm)) )**-2
            vel(1) = etam * C / (d + etam)
        else
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
            ! Weakly nonlinear dispersive theory    !
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
            ! We found that to get the appropriate wave height we must specify
            ! a different height via the following function...
		    !delta = H/d 
            delta = 0.2821d0 * ( H / d )**2 + 1.0147d0 * ( H / d ) - 0.0012d0
            ! Get the wave speed using the Pade approximate       
            ! Now get the required coefficients..
            c2 = g * d * (ONE + delta)
		    c  = sqrt(c2)
		    a  = (c2 - g * d ) / c
		    a1 = a * ( c - a ) / g
		    a2 = a * a / g
		    b  = sqrt(THREEQ * a / ( c * d * d ) )
            ! Below are the ones for Nwogu's formulation at zalpha
		    !a1 = d*(c2-g*d)/(3d0*((beta+1d0/3d0)*g*d-beta*c2))
		    !a2 = -d*(c2-g*d)**2/(2d0*g*d*c2)*((beta+1d0/3d0)*  &
     	    !	g*d+2d0*beta*c2)/((beta+1d0/3d0)*g*d-beta*c2)
		    !b  = sqrt((c2-g*d)/(4d0*((beta+1d0/3d0)*g*d**3-    &
     	    !	beta*d**2*c2)))
            ksi = HALF * ( xp + xm ) 
            ! Enter the final solution
            eta    = a1 / cosh(b * ksi)**2 + a2 / cosh(b * ksi)**4
            vel(1) = a / cosh(b * xm)**2
        endif
    else
        eta = ZERO
        vel = ZERO
    endif
    if (direc.eq.HCIRC) vel(1) = - vel(1)
!
endsubroutine init_height_vel
!=============================================================================!
!                            BOUNDARY_CONDITION:                              !
!    Calls boundary conditions at the NORTH, SOUTH, EAST and WEST boundaries  !
!                                                             PRINGLE (2015)  !
!=============================================================================!
subroutine boundary_condition(iflag,DT_CF)
    use variables, only: is, ie, js, je, ks, ke
    use arrays, only: L, dim3d, BC
    implicit none
    !----------------- Temporary variables ------------------------------------
    integer,intent(in) :: iflag
    real*8,intent(in) :: DT_CF
    integer :: iss, iee, jss, jee, kss, kee, DIM
    !
    !Call the automatic boundary condition maker
    if (iflag.eq.2) then
        iss = L(1)%is ; iee = L(1)%ie
        jss = L(1)%js ; jee = L(1)%je
        kss = 3       ; kee = 4
        DIM = L(1)%DIM
    elseif (iflag.eq.3) then
        iss = is ; iee = ie
        jss = js ; jee = je 
        kss = ks ; ; kee = ke
        if (dim3d) DIM = 3
        if (.not.dim3d) DIM = 2 
    else
        stop 'Bad iflag in boundary_condition'
    endif
!$omp parallel sections
!$omp section
    ! North
    if (jss.ne.jee-1) call auto_boundary(BC(iflag)%North,iss,iee,jee,         &
                                         kss,kee,DIM,'N',iflag,DT_CF)

!$omp section
    ! South
    if (jss.ne.jee-1) call auto_boundary(BC(iflag)%South,iss,iee,jss-1,       &
                                         kss,kee,DIM,'S',iflag,DT_CF)
!$omp section
    ! East
    if (iss.ne.iee-1) call auto_boundary(BC(iflag)%East,jss,jee,iee,          &
                                         kss,kee,DIM,'E',iflag,DT_CF)
!$omp section
    ! West
    if (iss.ne.iee-1) call auto_boundary(BC(iflag)%West,jss,jee,iss-1,        &
                                         kss,kee,DIM,'W',iflag,DT_CF)
!$omp end parallel sections    
end subroutine boundary_condition
!
subroutine auto_boundary(bound,iss,iee,jss,kss,kee,DIM,cflag,iflag,DT_CF)
    use interface_list, only: hn_sanders_obc 
    use variables, only: LNUM, ZERO, VSMALL, HCIRC, HALF, ONEHALF, inflow,    &
                         g, t, ks, ke
    use arrays, only: L, mn2d, nff, mn, u, dx, f, IWD, FSL
    implicit none
    !----------------- Temporary variables ------------------------------------
    character*4,intent(in) :: bound
    character*1,intent(in) :: cflag
    real*8,intent(in) :: DT_CF
    integer,intent(in) :: iflag, jss, iss, iee, kss, kee, DIM
    integer :: nn, nb, nbm, nbmm, ii, i, j, ix, iy, k
    real*8 :: wl, wlb, vb(DIM,kss:kee-1), mom
!
!=============================================================================!
!   Calculate the new conditions at the boundary:                             ! 
!   This routine may calculate BCs for inflow conditions,                     !
!   outflow/open conditions, reflective conditions, or a combination          !
!=============================================================================!
!
    wl = -1d4 ; vb = -1d4
    if (cflag.eq.'N') then
        ix =  0 ; iy = -1 ; j = jss
        if (inflow%angle.eq.ONEHALF*HCIRC.and.(LNUM.eq.0.or.iflag.eq.2)) then
            call inputwave(wl,vb,iflag,iy,DIM,jss,kss,kee)
        endif
    elseif (cflag.eq.'S') then
        ix =  0 ; iy =  1 ; j = jss 
        if (inflow%angle.eq.HALF*HCIRC.and.(LNUM.eq.0.or.iflag.eq.2)) then
            call inputwave(wl,vb,iflag,iy,DIM,jss,kss,kee)
        endif
    elseif (cflag.eq.'E') then
        ix = -1 ; iy =  0 ; i = jss 
        if (inflow%angle.eq.HCIRC.and.(LNUM.eq.0.or.iflag.eq.2)) then
            call inputwave(wl,vb,iflag,iy,DIM,jss,kss,kee)
        endif
    elseif (cflag.eq.'W') then
        ix =  1 ; iy =  0 ; i = jss 
        if (inflow%angle.eq.ZERO.and.(LNUM.eq.0.or.iflag.eq.2)) then
            call inputwave(wl,vb,iflag,iy,DIM,jss,kss,kee)
        endif
    endif  
!
    if (wl.eq.-1d4.and.bound.eq.'WALL') return
    
    DO_ALONG_BOUNDARY: do ii = iss,iee-1
!
    if (ix.eq.0) i = ii
    if (iy.eq.0) j = ii
    if (iflag.eq.2) then
        nb   = L(1)%mn(i,j) ; nbm = L(1)%mn(i+ix,j+iy)
        nbmm = L(1)%mn(i+max(ix,0),j+max(iy,0))
    endif
    wlb = wl
    if (wl.eq.-1d4) then
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !BOUNDARY CONDITIONS WITHOUT INFLOW
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (bound(2:2).eq.'P') then
            ! Call the open boundary formula
            mom = get_momentum()
            if (iflag.eq.2) then
                if (L(1)%Dn(abs(iy)+1,nbmm).gt.ZERO.and.                      &
                    L(1)%ETAn(nbm)-L(1)%ZK(nbm).gt.L(1)%Dmin) then
                    wlb = hn_sanders_obc(ZERO,L(1)%H(abs(iy)+1,nbmm),         &
                                              L(1)%ETAn(nbm),mom,g)
                else
                    cycle
                endif
            elseif (iflag.eq.3) then
                wlb = hn_sanders_obc(ZERO,IWD(mn2d(i+ix,j+iy)),               &
                                          FSL(mn2d(i+ix,j+iy)),mom,g)
                if (abs(wlb).lt.VSMALL) wlb = ZERO
                !Calling routine to change the water level
                !to hydrostatic pressure and save as F value
                call wl_to_hp(wlb,FSL(mn2d(i+ix,j+iy)),i,j,ix,iy)
                cycle
            else
                cycle
            endif
            !
        elseif (bound.eq.'FREE') then
            !Call the free outflow subroutine
            if (iflag.eq.2) stop 'Bad: FREE boundary condition in 2D domain'
            call free_outflow(i,j,ix,iy)
            cycle
        elseif (bound(1:2).eq.'2D') then
            if (iflag.eq.2) stop 'Bad: "2D" boundary condition for 2D domain'
            call Interp_2D_bound_to_3D( ii, ix, iy, i, j,                     &
                                        L(LNUM)%DIM, DT_CF, bound(3:4) )
            cycle
        endif 
    else
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !BOUNDARY CONDITIONS WITH INFLOW
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !Adjust this value when considering open boundary condition
        if (bound(2:2).eq.'P') then
            !For this condition we calculate based on a
            !combination of the inflow and the sponge layer
            mom = get_momentum()
            if (iflag.eq.2) then
                if (L(1)%Dn(abs(iy)+1,nbmm).gt.ZERO.and.                      &
                    L(1)%ETAn(nbm)-L(1)%ZK(nbm).gt.L(1)%Dmin) then
                    wlb = hn_sanders_obc(wlb,L(1)%H(abs(iy)+1,nbmm),          &
                                             L(1)%ETAn(nbm),mom,g)
                endif
            elseif (iflag.eq.3) then
                wlb = hn_sanders_obc(wlb,IWD(mn2d(i+ix,j+iy)),                &
                                         FSL(mn2d(i+ix,j+iy)),mom,g)
                if (inflow%wavetype.ne.'RW'.and.t.gt.inflow%period) then
                    if (abs(wlb).lt.VSMALL) wlb = ZERO
                    !Calling routine to change the water level
                    !to hydrostatic pressure and save as F value
                    call wl_to_hp(wlb,FSL(mn2d(i+ix,j+iy)),i,j,ix,iy)
                    cycle
                endif
            endif
        endif 
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !RESETTING THE BOUNDARY WATER LEVEL
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (abs(wlb).lt.VSMALL) wlb = ZERO
    if (iflag.eq.2) then
        L(1)%ETAn(nb) = wlb
        ! Calculated inflow flux
        call calculate_flux(wl,wlb,vb(abs(iy)+1,kss),ix,iy,nbmm,nbm,nb)
        ! For the wavemaker condition
        !if (wl.ne.-1d4) then
        !    if (inflow%wavetype(1:1).eq.'M') call adjust_dx_for_wavemaker
        !endif
    elseif (iflag.eq.3.and.wl.ne.-1d4) then
        call wl_to_uf(wlb,vb,DIM,i,j,ix,iy)
    endif
    !
    enddo DO_ALONG_BOUNDARY   
!
    contains
!
    real*8 function get_momentum()
        implicit none
        integer :: k
        !
        if (iflag.eq.2) then
            get_momentum = L(1)%Qn(abs(iy)+1,nbmm)
        elseif (iflag.eq.3) then
            get_momentum = ZERO
            do k = ks, ke-1
                if (nff(mn(i+ix,k,j+iy))%f.eq.-1) cycle
                if (nff(mn(i+ix,k,j+iy))%f.eq.0) exit
                get_momentum = get_momentum                                   &
                             + u(abs(iy),mn(i+max(ix,0),k,j+max(iy,0)))       &
                             * dx(2,k) * f(mn(i+ix,k,j+iy))
            enddo 
        endif
        if (ix.eq.-1.or.iy.eq.-1) get_momentum = - get_momentum
    endfunction get_momentum
!
endsubroutine auto_boundary
!
subroutine inputwave(wl,vel,iflag,iy,DIM,jss,kss,kee)
    use interface_list, only: set_wl
    use variables, only: VSMALL, ZERO, g, t, ks, ke, HALF, inflow
    use arrays, only: L, x, dx
    implicit none
    integer :: jj
    integer,intent(in) :: iflag, iy, DIM, jss, kss, kee
    real*8,intent(out) :: wl, vel(DIM,kss:kee-1) !Water level and velocity 
    real*8 :: calc_erf_vel
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Gets the current wave height and velocity based on the input wave type  ! 
    !                                            - See UserManual Sec 2.3.2   !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !    
    if (inflow%wavetype.eq.'SW') then 
        ! Calculating with solitary wave
        if (iflag.eq.2) then
            call sw_eta_vel(wl,vel,inflow%period,inflow%h_wave,inflow%h_shoki,&
                            t,ZERO,ZERO,g,3,4,ZERO)
        elseif (iflag.eq.3) then
            call sw_eta_vel(wl,vel,inflow%period,inflow%h_wave,inflow%h_shoki,&
                            t,-HALF*dx(abs(iy),jss),ZERO,g,ks,ke,x(2,ks:ke))
        endif
    elseif (inflow%wavetype.eq.'RW') then 
        ! Calculating with regular wave
        call rw_eta_vel(wl,vel,inflow%period,inflow%h_wave,                   &
                        inflow%h_shoki,inflow%k,t,ZERO,g)
    elseif (inflow%wavetype(1:1).eq.'A'.or.inflow%wavetype(1:1).eq.'M') then
        ! Calculating from data file (arbitrary wave)                           
        wl = set_wl(t+inflow%t_adj,inflow%dt,inflow%size,inflow%hwave,2)
                                      !interpolating polynomial order ^   
        if (inflow%width.gt.1.and.inflow%width.lt.4) then
            do jj = 1,inflow%width-1
                vel(jj,kss) = set_wl(t+inflow%t_adj,inflow%dt,          &
                                     inflow%size,inflow%vwave(jj,:),2)
            enddo
        elseif (inflow%wavetype(2:2).eq.'E') then
            ! Calculate from error function
            vel(1,kss) = calc_erf_vel(inflow%period,t)
        endif   
    elseif (inflow%wavetype.eq.'UF') then
        wl = inflow%h_shoki ; vel = inflow%h_wave 
    endif
    if (abs(wl).lt.VSMALL)     wl  = ZERO
    where (abs(vel).lt.VSMALL) vel = ZERO        
endsubroutine inputwave
!
subroutine calculate_flux(wl,wlb,vb,ix,iy,nbmm,nbm,nb)
    use arrays, only: L
    use variables, only: VSMALL, g, inflow
    implicit none
    integer,intent(in) :: ix, iy, nbmm, nbm, nb
    real*8,intent(in)  :: wl, wlb, vb                 !Water level and velocity 
    integer            :: isign
    !
    if (iy.eq.0) isign = ix
    if (ix.eq.0) isign = iy
    if (vb.eq.-1d4.or.wl.eq.-1d4.or.                                          &
       (wl.ne.-1d4.and.abs(wl-wlb).gt.VSMALL)) then
        ! Open flow condition with or without inflow -> Calculated from LSW
        ! -or- inflow with no velocity data
        L(1)%Qn(abs(iy)+1,nbmm) = L(1)%Q(abs(iy)+1,nbmm) + isign * L(1)%DT    &
                                * g * L(1)%H(abs(iy)+1,nbmm)                  &
                              * (wlb - L(1)%ETAn(nbm)) * L(1)%RXc(abs(iy)+1,nb)
    else
        ! Just inflow condition already given by vb in input file
        if (inflow%wavetype(2:2).eq.'I') then
            L(1)%Qn(abs(iy)+1,nbmm) = vb
            !if (inflow%width.gt.2) then
            !    L(1)%Qn(abs(ix)+1,nb) = vb(2)
            !endif
        else
            ! Just inflow condition -> calculated from vb and wl
            L(1)%Qn(abs(iy)+1,nbmm) = isign * vb                              &
                                    * (wl + L(1)%H(abs(iy)+1,nbmm))
            !if (inflow%width.gt.2) then
            !    L(1)%Qn(abs(ix)+1,nb) = vb(2) * (wl + L(1)%H(abs(ix)+1,nb))
            !endif
        endif
    endif
endsubroutine calculate_flux
!
subroutine wl_to_hp(wl,wlm,i,j,ix,iy) 
    use variables, only: g, ks, ke, h_in, ZERO, ONE
    use arrays, only: nff, x, f, p, fp, mn
    implicit none
    integer,intent(in) :: i, j, ix, iy
    real*8,intent(in) :: wl, wlm
    real*8 :: pdiff
    integer :: k
!=============================================================================!
!     Converts the calculated water level Hydrostatic P and F fraction        !
!=============================================================================!
!
pdiff = g * (wl-wlm)
do k = ks,ke-1
    if (nff(mn(i,k,j))%b.ne.4) nff(mn(i,k,j))%b = 4
    if (wl + h_in.ge.x(2,k+1)) then
        f(mn(i,k,j)) = ONE ; fp(mn(i,k,j)) = ONE
        p(mn(i,k,j)) = p(mn(i+ix,k,j+iy)) + pdiff
    elseif (x(2,k).ge.wl+h_in) then
        f(mn(i,k,j)) = ZERO ; fp(mn(i,k,j)) = ZERO
        p(mn(i,k,j)) = p(mn(i+ix,k,j+iy))
    else
        f(mn(i,k,j))  = ( wl + h_in - x(2,k) ) / ( x(2,k+1) - x(2,k) )
        fp(mn(i,k,j)) = ( wl + h_in - x(2,k) ) / ( x(2,k+1) - x(2,k) )
        p(mn(i,k,j))  = p(mn(i+ix,k,j+iy)) + pdiff
    endif        
enddo
!
endsubroutine wl_to_hp
!
subroutine wl_to_uf(wl,vel,DIM,i,j,ix,iy)
    use variables, only: ZERO, h_in, HALF, cm, ks, ke, ONE, inflow,           &
                         nut_b, nu, TI_b, akini, eini
    use arrays, only: mn2D, IWD, FSL, mn, a, nff, x, un, f, fp, dx, r, r1, d
    use interface_list, only: wdiv
    implicit none
    integer,intent(in) :: i, j, ix, iy, DIM
    real*8,intent(in) :: wl,vel(DIM,ks:ke-1)
    integer :: k, nnmx, nnpx, nn_1, nn_2 
    real*8 :: adj, znow, zadj, c1
    !
    !=========================================================================!
    ! Converts the calculated water level to U and F fraction                 !
    ! Updated to consider input velocity: Nov 10 2014                         !
    ! If we have some input velocity, set F fraction at nf = -1 column        !
    ! If we have no input velocity set F fraction at first cell column        !
    !=========================================================================!
    !
    ! Use adjustment U and z based on: 
    ! Kawasaki, K (2013). Chapter 4 in Computational Wave Dynamics, pg90.
    ! Adjust horizontal velocity to match with calculation
    adj = ( wl + IWD(mn2d(i+ix,j+iy)) )                                       &
        / ( FSL(mn2d(i+ix,j+iy)) + IWD(mn2d(i+ix,j+iy)) )
    do k = ks,ke-1
        nn_2 = mn(i+ix,k,j+iy)
        if (nff(nn_2)%f.eq.-1) cycle
        nnmx = mn(i+max(ix,0),k,j+max(iy,0))
        if (a(abs(iy),nnmx).eq.ZERO) cycle
        nn_1 = mn(i,k,j)
        if (nff(nn_1)%b.ne.2) nff(nn_1)%b = 2
        nnpx = mn(i+ix+max(ix,0),k,j+iy+max(iy,0))
        if (x(2,k+1).le.wl + h_in) then
            znow = HALF * (x(2,k) + x(2,k+1))
            zadj = adj * (znow + IWD(mn2d(i+ix,j+iy))) - IWD(mn2d(i+ix,j+iy))
            if (zadj > znow.and.k.lt.ke-1) then
                c1 = ( zadj - znow ) / ( HALF * (x(2,k+1) + x(2,k+2)) - znow)
                un(abs(iy),nnmx) = (vel(1,k+1) * c1 +                         &
                                    vel(1,k) * (ONE - c1) ) * adj
            elseif (zadj < znow.and.k.gt.ks) then  
                c1 = ( znow - zadj ) / ( znow - HALF * (x(2,k) + x(2,k-1)) )
                un(abs(iy),nnmx) = (vel(1,k-1) * c1 +                         &
                                    vel(1,k) * (ONE - c1) ) * adj  
            else
                un(abs(iy),nnmx) = vel(1,k) * adj
            endif
            un(2,nn_1)       = vel(DIM,k) 
            f(nn_1) = ONE ; fp(nn_1) = ONE
        elseif (x(2,k).ge.wl+h_in) then
            un(abs(iy),nnmx) = ZERO
            un(2,nn_1)       = ZERO
            f(nn_1) = ZERO ; fp(nn_1) = ZERO
        else
            !Surface cell      ! !vel(1,k) * adj
            un(abs(iy),nnmx) = un(abs(iy),mn(i+max(ix,0),k-1,j+max(iy,0)))
            un(2,nn_1)       = vel(DIM,k)
            f(nn_1)  = ( wl + h_in - x(2,k) ) / dx(2,k)
            fp(nn_1) =  f(nn_1) 
        endif
        if (inflow%wavetype.eq.'UF') then
            f(nn_2) = f(nn_1) ; fp(nn_2) = f(nn_1) 
        endif
#ifdef RAN !Set k-eps boundary condition as defined by turbulent intensity
        r(nn_1)%k  = akini !HALF * ( un(abs(iy),nnmx) * TI_b )**2  !max(,r1(nn_2)%k)
        d(nn_1)    = nut_b * nu      !wdiv( r1(nn_1)%k**2 , r1(nn_1)%e ) * cm
        r(nn_1)%e  = eini  !cm * r1(nn_1)%k**2 / d(nn_1) !max(,r1(nn_2)%e)
        r1(nn_1)%k = r(nn_1)%k ; r1(nn_1)%e  = r(nn_1)%e
#endif
    enddo
endsubroutine wl_to_uf
!
subroutine free_outflow(i,j,ix,iy)
    use arrays, only: mn, nff, a, un, r, r1, d, u
    use variables, only: ks, ke, cm, ZERO, TI_b, nut_b, nu, HALF, akini, eini
    use interface_list, only: wdiv
    implicit none
    integer,intent(in) :: i, j, ix, iy
    integer :: nnmx, k, nn_1, nn_2
!
    do k = ks,ke-1
        nn_1 = mn(i,k,j) ; nn_2 = mn(i+ix,k,j+iy)
        if (nff(nn_2)%f.eq.-1) cycle
        nnmx = mn(i+max(ix,0),k,j+max(iy,0))
        if (a(abs(iy),nnmx).eq.ZERO) cycle
        if ((nff(nn_1)%b.eq.2.or.nff(nn_1)%b.eq.3)) then
            !Set the water velocity at boundary equal 
            !to adjacent value for outgoing flow only
            if (ix.eq.-1.or.iy.eq.-1) then
                un(abs(iy),nnmx) = max(ZERO,u(abs(iy),nn_2))
            else
                un(abs(iy),nnmx) = min(ZERO,u(abs(iy),nn_2))
            endif
        endif
#ifdef RAN !Set k-eps boundary condition as defined by turbulent intensity
        r(nn_1)%k  = akini !HALF * ( un(abs(iy),nnmx) * TI_b )**2  !max(,r1(nn_2)%k)
        d(nn_1)    = nut_b * nu      !wdiv( r1(nn_1)%k**2 , r1(nn_1)%e ) * cm
        r(nn_1)%e  = eini  !cm * r1(nn_1)%k**2 / d(nn_1) !max(,r1(nn_2)%e)
        r1(nn_1)%k = r(nn_1)%k ; r1(nn_1)%e  = r(nn_1)%e
#endif
    enddo  
endsubroutine free_outflow
!
!subroutine adjust_dx_for_wavemaker
!    use leapfrog_mod, only: PI_8, A, EPS, TWO, ONE, ZERO, HALF
!    real*8  :: m0, Xs, Xp, wm_pos, set_wl, round, ADD_E, ADD_M, DXM, DXP
!    integer :: jj, iee, nn_1, nn_2, nn_3
!!------------------------------------------------------------------------------!
!!   Adjusts the DX or DY values according to the transformation:               !
!!       T(x) = x` = m0/m(t)*(x - xp(t)) : i.e multiply by m0/m(t)              !
!!   Reference : Orszaghova, Jana (2011). Solitary waves and wave groups at the !
!!               shore. Oxford Uni. Doctorate Thesis.                           ! 
!!-------------------------------------------------------------------------------
!! Set edge of paddle domain to 10 times the maximum wavemaker stroke
!m0 = 10d0 !10d0 * inflow(nn)%m_wave
!! Find the current position of the wavemaker
!if (inflow(nn)%width.eq.3) then
!    wm_pos = set_wl(t,inflow(nn)%dt,inflow(nn)%size,inflow(nn)%mwave,2) 
!    Xp = wm_pos - inflow(nn)%mwave(1)
!elseif (inflow(nn)%wavetype(2:2).eq.'E') then 
!    !Calculate from error function                  
!    wm_pos = erf( TWO * PI_8 * (t - inflow(nn)%t_adj) / inflow(nn)%period )
!    Xp = wm_pos + ONE
!endif
!if (vb(1).eq.ZERO.and.inflow(nn)%m_wave-Xp.lt.eps) return !No need to update
!if (iy.eq.0) then
!! Change DX
!    Xs = L(1)%X(i+max(ix,0))
!    if (ix.eq.1) iee = L(1)%ie-1 ;  if (ix.eq.-1) iee = L(1)%is
!    do jj = i,iee,ix
!        ! Don't change DX if beyond the paddle domain
!        if (abs(L(1)%X(jj) - Xs).gt.m0) exit
!        nn_1 = L(1)%mn(jj,ii) ; nn_2 = L(1)%mn(jj-1,ii) 
!        L(1)%DX(1,nn_1) = (m0 - Xp)*round(L(1)%X(jj+1)-L(1)%X(jj),L(1)%DEC) / m0 
!        L(1)%RX(1,nn_1) = ONE / L(1)%DX(1,nn_1) 
!        L(1)%RXc(1,nn_1)= TWO / (L(1)%DX(1,nn_1) + L(1)%DX(1,nn_2)) 
!    enddo
!    do jj = i+max(ix,0),iee,ix
!        if (abs(L(1)%X(jj) - Xs).gt.m0) exit
!        nn_1 = L(1)%mn(jj,ii) ; nn_2 = L(1)%mn(jj-1,ii) ;nn_3 = L(1)%mn(jj+1,ii)
!        if (L(1)%Dn(1,nn_1).eq.ZERO) cycle
!        ADD_E = (m0 + Xs - HALF * (L(1)%X(jj+1) + L(1)%X(jj)))                 &
!              * L(1)%DT * ix * vb(1) * L(1)%RXc(1,nn_1) / (m0 - Xp)
!        ADD_M = (m0 + Xs - L(1)%X(jj)) * L(1)%DT * ix * vb(1)                  &
!              * L(1)%RX(1,nn_2)  / (m0 - Xp)
!        ! Use first order upwind discretization..
!        if (real(ix)*vb(1).gt.ZERO) then
!            L(1)%ETA(nn_1) = L(1)%ETA(nn_1) - ADD_E                            &
!                           * (L(1)%ETAn(nn_1) - L(1)%ETAn(nn_2))
!            L(1)%Q(1,nn_1) = L(1)%Q(1,nn_1) - ADD_M                            &
!                           * (L(1)%Qn(1,nn_1) - L(1)%Qn(1,nn_2))
!        elseif (real(ix)*vb(1).lt.ZERO) then
!            L(1)%ETA(nn_1) = L(1)%ETA(nn_1) - ADD_E                            &
!                           * (L(1)%ETAn(nn_3) - L(1)%ETAn(nn_1))
!            L(1)%Q(1,nn_1) = L(1)%Q(1,nn_1) - ADD_M                            &
!                           * (L(1)%Qn(1,nn_3) - L(1)%Qn(1,nn_1))
!        endif
!        if (L(1)%ETA(nn_1).lt.L(1)%ZK(nn_1)) L(1)%ETA(nn_1) = L(1)%ZK(nn_1)
!        if (L(1)%DISP.ne.2) cycle
!        if (L(1)%ZK(nn_1).ge.ZERO) cycle
!        ! If using implicit dispersion correcter 
!        ! we need to update the matrix divisors
!        DXM = HALF * ( L(1)%DX(1,nn_2) + L(1)%DX(1,nn_1) )
!        DXP = HALF * ( L(1)%DX(1,nn_3) + L(1)%DX(1,nn_1) )
!        L(1)%DIAG(nn_1) = -ONE - A * TWO * L(1)%ZK(nn_1) * L(1)%ZK(nn_1)       &
!                        * (( DXM + DXP ) / ( DXM * ( DXP * DXP + DXM * DXP ) ) )
!        if (L(1)%ZK(nn_2).lt.ZERO) then
!            L(1)%AL(nn_1,1) = A * TWO * L(1)%ZK(nn_1) * L(1)%ZK(nn_1)          &
!                            * DXP / ( DXM * ( DXP * DXP + DXM * DXP ) )
!        endif
!        if (L(1)%ZK(nn_3).lt.ZERO) then
!            L(1)%AL(nn_1,2) = A * TWO * L(1)%ZK(nn_1) * L(1)%ZK(nn_1)          &
!                            * DXM / ( DXM * ( DXP * DXP + DXM * DXP ) )
!        endif
!    enddo
!elseif (ix.eq.0) then
!! Change DY
!endif
!endsubroutine adjust_dx_for_wavemaker
!
!------------------------------------------------------------------------------
subroutine calc_IC_or_GC(cflag,iflag,LN)
    use arrays, only: L, FSL, IWD, U, mn, x, mn2d, a, eta_in, vel_in
    use variables, only: LNUM, ONE, ZERO, HALF, t, is, ie, js, je, ks, ke,    &
                         inflow, inflow_file, GDR
	use interface_list, only: set_wl , LI_1
    implicit none
    !----------------- Input variables ----------------------------------------
    integer,intent(in) :: LN, iflag
    character*2,intent(in) :: cflag
    !----------------- Temporary variables ------------------------------------
    integer :: nn, i, j, k, ii, nn_out, ix, iy, unit, na, LNN
    integer :: DIM, iss, iee, jss, jee, kss, kee, LT, nn_out1
    real*8  :: Xs, Ys, Xe, Ye, Xc, Yc, t_now, t_pre, eta_now, eta_pre
	real*8  :: c1, c2, c3, c4, change
    character*1 :: ch1
    character*2 :: ch2
    character*255 :: inflow_file_temp
    logical :: exist
!
!=============================================================================!
!   Calculate the new conditions at the boundary:                             ! 
!   This routine may calculate BCs for inflow conditions,                     !
!   outflow/open conditions, reflective conditions, or a combination          !
!=============================================================================!
!
    !Do not need to input IC or GC
    if (LN.le.LNUM) then
        if (L(LN)%IC.eq.'NO') return
        if (inflow%wavetype.ne.cflag.and.L(LN)%IC.ne.cflag) return
    else
        if (inflow%wavetype.ne.cflag) return
    endif
    if (inflow%wavetype(2:2).eq.'T') then
        LT = int(inflow%T_end/inflow%dt) + 1
    else
        LT = 1
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !Read the initial water level & velocity values
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !Check to see if the the inflow file exists for this layer..
    na = len_trim(inflow_file)
    if (LN.lt.10) then
        write(ch1,'(I1)') LN
        inflow_file_temp = inflow_file(1:na-min(na-1,4))//'_'//ch1//'.dat'
    elseif (LN.ge.10) then
        write(ch2,'(I2)') LN
        inflow_file_temp = inflow_file(1:na-min(na-1,4))//'_'//ch2//'.dat'
    endif
    inquire(file=inflow_file_temp,exist=exist)
    !
    if (LN.eq.1.or.exist) then
        if (allocated(eta_in)) deallocate(eta_in)
        if (allocated(vel_in)) deallocate(vel_in)
        DIM = L(LN)%DIM
        iss = L(LN)%is; iee = L(LN)%ie ; kss = ks
        jss = L(LN)%js; jee = L(LN)%je
 	    IFLAG_IF: if (iflag.eq.0) then   
        allocate(L(LN)%eta_in(iss:iee-1,jss:jee-1,LT)) 
        if (inflow%order.eq.'EQ') then  
            allocate(vel_in(DIM,iss:iee-1,jss:jee-1,kss:kss))
            do i = iss,iee-1
                do j = jss,jee-1
                    call init_height_vel(eta_in(i,j),vel_in(:,i,j,kss),DIM,nn,&
                                         i,j,inflow%h_wave,inflow%h_shoki,    & 
                                         inflow%x0(0),inflow%x0(1),           &
                                         kss,kee,inflow%angle)  
                enddo
            enddo
        else
            if (exist) then
                open(newunit = unit,file = inflow_file_temp,status = 'old')
            else
                open(newunit = unit,file = inflow_file,status = 'old')
            endif
			do k = 1,LT
            if (inflow%order(2:2).eq.'P') then
                do j = jss,jee-1
                    if (inflow%fmt1(1:1).eq.'*') then
                        if (inflow%order(1:1).eq.'P') then
                            read(unit,*) (L(LN)%eta_in(i,j,k), i = iss, iee-1 )  
                        elseif (inflow%order(1:1).eq.'N') then
                            read(unit,*) (L(LN)%eta_in(i,j,k), i = iee-1, iss, -1 )   
                        endif
                    else
                        if (inflow%order(1:1).eq.'P') then
                            read(unit,inflow%fmt1) (L(LN)%eta_in(i,j,k),i = iss,iee-1 )  
                        elseif (inflow%order(1:1).eq.'N') then
                            read(unit,inflow%fmt1) (L(LN)%eta_in(i,j,k),i=iee-1,iss,-1) 
                        endif
                    endif
                enddo
            elseif (inflow%order(2:2).eq.'N') then
                do j = jee-1,jss,-1
                    if (inflow%fmt1(1:1).eq.'*') then
                        if (inflow%order(1:1).eq.'P') then
                            read(unit,*) (L(LN)%eta_in(i,j,k), i = iss, iee-1 )  
                        elseif (inflow%order(1:1).eq.'N') then
                            read(unit,*) (L(LN)%eta_in(i,j,k), i = iee-1, iss, -1 )  
                        endif
                    else
                        if (inflow%order(1:1).eq.'P') then
                            read(unit,inflow%fmt1) (L(LN)%eta_in(i,j,k),i = iss,iee-1 )  
                        elseif (inflow%order(1:1).eq.'N') then
                            read(unit,inflow%fmt1) (L(LN)%eta_in(i,j,k),i=iee-1,iss,-1)
                        endif
                    endif
                enddo
            endif
            enddo     
            close(unit)
        endif
        if (GDR /= ONE) then
          eta_in = eta_in * GDR
          write(*,'(A,F5.1)') 'Ground Deformation Rate = ',GDR
        endif 
        
        if (maxval(L(LN)%eta_in).gt.inflow%h_wave) then
            inflow%h_wave = maxval(L(LN)%eta_in)
        endif
        if (cflag(1:1).eq.'I') then
            write(6,*) 'LN=',LN,'Max. Water Level of IC= ',maxval(L(LN)%eta_in),'[m]'  
        elseif (cflag(1:1).eq.'G') then
            write(6,*) 'LN=',LN,'Max. Change in GL= ', maxval(L(LN)%eta_in),'[m]'  
        endif
        if (allocated(vel_in)) then
            if (maxval(vel_in).gt.inflow%v_wave) then
                inflow%v_wave = maxval(vel_in)
            endif
            write(6,*) 'LN=',LN,'Max. Velocity of IC= ', maxval(vel_in),'[m/s]' 
        endif
        endif IFLAG_IF
        if (iflag.eq.0.and.inflow%wavetype(2:2).eq.'T') return
        !Now insert the read heights into the data 
        allocate(eta_in(iss:iee-1,jss:jee-1))
        if (inflow%wavetype(2:2).eq.'T') then
        t_now = min(inflow%T_end,t + HALF*L(LN)%DT)
        t_pre = max(ZERO,t - HALF*L(LN)%DT)
            do i = iss,iee-1
                do j = jss,jee-1  
                    eta_now = set_wl(t_now,inflow%dt,LT,L(LN)%eta_in(i,j,:),1)
                    eta_pre = set_wl(t_pre,inflow%dt,LT,L(LN)%eta_in(i,j,:),1)
                    eta_in(i,j) = eta_now - eta_pre
                enddo
            enddo
        else
            eta_in(:,:) = L(LN)%eta_in(:,:,1) 
        endif
!!$omp parallel do default(none)                                                &
!!$omp shared(L,LN,LNN,LNUM,is,ie,js,je,iss,jss,iee,jee,i                       &
!!$omp        X,eta_in,cflag,FSL,mn2D,vel_in,a,u,mn,IWD,DIM,kss,kee)            &
!!$omp private(j,nn_out,nn_out1,ix,iy,ii)
        do i = iss,iee-1
            do j = jss,jee-1
                if (LN.le.LNUM) then
                    nn_out = L(LN)%mn(i,j)
                    if (nn_out.eq.0) cycle
                endif
                !Initial condition for outermost layer
                if (cflag(1:1).eq.'I') then
                    !Initial water level condition
                    if (LN.le.LNUM) then
                        L(LN)%ETA(nn_out) = L(LN)%ETA(nn_out) + eta_in(i,j)
                    elseif (LN.gt.LNUM) then
                        FSL(mn2d(i,j)) = FSL(mn2d(i,j)) + eta_in(i,j)
                        call FSL_to_F(FSL(mn2d(i,j)),i,j)
                    endif
                    !Initial velocity condition
                    if (allocated(vel_in)) then
                        if (LN.le.LNUM) then
                            L(LN)%U(:,nn_out) = vel_in(:,i,j,kss)
                        elseif (LN.gt.LNUM) then
                            !x vel
                            U(0,mn(i,kss:kee-1,j)) = vel_in(1,i,j,kss:kee-1)
                            !z vel
                            if (DIM.eq.2) then
                                U(2,mn(i,kss:kee-1,j)) = vel_in(2,i,j,kss:kee-1)
                            elseif (DIM.eq.3) then
                                !y vel
                                U(1,mn(i,kss:kee-1,j)) = vel_in(2,i,j,kss:kee-1)
                                !z vel
                                U(2,mn(i,kss:kee-1,j)) = vel_in(3,i,j,kss:kee-1)
                            endif
                            ! Velocity can't be non-zero at closed boundaries
                            where (a(:,mn(i,kss:kee-1,j)).eq.ZERO)            & 
                                   U(:,mn(i,kss:kee-1,j)) = ZERO
                        endif
                    endif
                elseif (cflag(1:1).eq.'G') then 
                    !Change in ground level condition
                    if (L(LN)%IC(1:1).eq.'G') then
                        L(LN)%ZK(nn_out) = min(L(LN)%Z(L(LN)%ke),             &
                                               L(LN)%ZK(nn_out) + eta_in(i,j))
                        !Change in teibou condition
                        do ii = 1,L(LN)%DIM
                            ix = 0 ; iy = 0
                            if (ii.eq.1.and.i.ne.iss) ix = -1
                            if (ii.eq.2.and.j.ne.jss) iy = -1
                            nn_out1 = L(LN)%mn(i+ix,j+iy)
                            if (-L(LN)%H(ii,nn_out).gt.L(LN)%ZKm(nn_out)  &
                                .and.-L(LN)%H(ii,nn_out).gt.L(LN)%ZKm(nn_out1)) then
                                L(LN)%H(ii,nn_out) = max(-L(LN)%Z(L(LN)%ke),  &
                                                    L(LN)%H(ii,nn_out)    &
                                                    - HALF * (eta_in(i,j) &
                                                        + eta_in(i+ix,j+iy)))
                            elseif (-L(LN)%H(ii,nn_out).eq.L(LN)%ZKm(nn_out).or.&
                                    -L(LN)%H(ii,nn_out).eq.L(LN)%ZKm(nn_out1)) then 
                                L(LN)%H(ii,nn_out) = -ONE * max(L(LN)%ZK(nn_out1),L(LN)%ZK(nn_out)) 
                            else
                                L(LN)%H(ii,nn_out) = - LI_1( L(LN)%ZK(nn_out1),L(LN)%ZK(nn_out),&
                                                             L(LN)%DX(ii,nn_out1),L(LN)%DX(ii,nn_out) ) 
                            endif
                        enddo 
                    endif
                endif
            enddo
        enddo
!!$omp end parallel do
        return
    endif
	if (iflag.eq.0.and.inflow%wavetype(2:2).eq.'T') return
    do LNN = LN-1,1,-1
        if (LNN.lt.10) then
            write(ch1,'(I1)') LNN
            inflow_file_temp = inflow_file(1:na-min(na-1,4))//'_'//ch1//'.dat'
        elseif (LNN.ge.10) then
            write(ch2,'(I2)') LNN
            inflow_file_temp = inflow_file(1:na-min(na-1,4))//'_'//ch2//'.dat'
        endif
        inquire(file=inflow_file_temp,exist=exist)
        if (exist) exit
    enddo
    if (inflow%wavetype(2:2).eq.'T') then
        if (allocated(eta_in)) deallocate(eta_in)
        allocate(eta_in(L(LNN)%is:L(LNN)%ie-1,L(LNN)%js:L(LNN)%je-1))
        t_now = min(inflow%T_end,t + HALF*L(LN)%DT)
        t_pre = max(ZERO,t - HALF*L(LN)%DT)
        do i = L(LNN)%is,L(LNN)%ie-1
            do j = L(LNN)%js,L(LNN)%je-1
                if (LN.le.LNUM) then
                    ! Interpolate in 2DH section
                    eta_now = set_wl(t_now,inflow%dt,LT,L(LNN)%eta_in(i,j,:),1)
                    eta_pre = set_wl(t_pre,inflow%dt,LT,L(LNN)%eta_in(i,j,:),1)
                    eta_in(i,j) = eta_now - eta_pre
                else
                    ! Just get total deformation for 3D section
                    eta_in(i,j) = L(LNN)%eta_in(i,j,LT)
                endif
            enddo
        enddo
    endif
    call interp_onto_inner_layer(cflag,LNN,LN)
endsubroutine calc_IC_or_GC
!
subroutine interp_onto_inner_layer(cflag,LNN,LN)
    use arrays, only: L, FSL, IWD, U, mn, x, mn2d, a, eta_in, vel_in
    use variables, only: LNUM, ONE, ZERO, HALF, t, is, ie, js, je, ks, ke, SMALL
    use interface_list, only: LI_1
    implicit none
    !----------------- Input variables ----------------------------------------
    integer,intent(in) :: LN, LNN
    character*2,intent(in) :: cflag
    !----------------- Temporary variables ------------------------------------
    integer :: i, j, k, ii, nn_in, i2, j2, ix, iy
    integer :: DIM, iss, iee, jss, jee, kss, kee, nn_in1
    real*8  :: Xs, Ys, Xe, Ye, Xc, Yc, Xx, Yy
	real*8  :: c1, c2, c3, c4, c1x, c2x, c3y, c4y, change
    !--------------------------------------------
    !  We have to interpolate the data for the  ! 
    !  last layer (LNN) onto the inner layer    !
    !--------------------------------------------
!!$omp parallel do default(none)                                                &
!!$omp shared(L,LN,LNN,LNUM,is,ie,js,je,X,eta_in,cflag,FSL,mn2D,vel_in,u,mn,IWD)&
!!$omp private(j,Xs,Ys,Xe,Ye,iss,jss,iee,jee,Xc,Xx,c1,c2,c3,c4,nn_in,nn_in1,    &
!!$omp         Yc,Yy,change,DIM,kss,kee,ix,iy,c1x,c2x,c3y,c4y)
    do i = L(LNN)%is,L(LNN)%ie-1
        do j = L(LNN)%js,L(LNN)%je-1
            !We know we can cycle at the first layer edge
            if (i.lt.L(LNN)%is2-1.or.i.ge.L(LNN)%ie2.or.                      & 
                j.lt.L(LNN)%js2-1.or.j.ge.L(LNN)%je2) cycle 
            Xs = HALF * ( L(LNN)%X(i) + L(LNN)%X(i+1) )
            Ys = HALF * ( L(LNN)%Y(j) + L(LNN)%Y(j+1) )
            Xe = HALF * ( L(LNN)%X(i+1) + L(LNN)%X(i+2) )
            Ye = HALF * ( L(LNN)%Y(j+1) + L(LNN)%Y(j+2) )
            !For LN > 2 we need to cycle here..
            if (LN.le.LNUM) then
                iss = L(LN)%is; iee = L(LN)%ie
                jss = L(LN)%js; jee = L(LN)%je
                kss = 3; kee = 4;
                if (Xe.lt.L(LN)%X(iss).or.Ye.lt.L(LN)%Y(jss)) cycle
                if (Xs.gt.L(LN)%X(iee).or.Ys.gt.L(LN)%Y(jee)) cycle
            else
                iss = is; iee = ie
                jss = js; jee = je
                kss = ks; kee = ke;
                if (Xe.lt.X(0,iss).or.Ye.lt.X(1,jss)) cycle
                if (Xs.gt.X(0,iee).or.Ys.gt.X(1,jee)) cycle
            endif
            !Loop over the inner layer
            do i2 = iss, iee-1
                if (LN.le.LNUM) then
                    Xc = HALF * ( L(LN)%X(i2) + L(LN)%X(i2+1) )
                    Xx = L(LN)%X(i2)
                elseif (LN.gt.LNUM) then
                    Xc = HALF * ( X(0,i2) + X(0,i2+1) )
                    Xx =  X(0,i2)
                endif
                if (Xc.gt.Xs.and.Xc.le.Xe) then
                    c1 = ( Xe - Xc ) / ( Xe - Xs )  
                    c2 = ( Xc - Xs ) / ( Xe - Xs ) 
                    do j2 = jss, jee-1
                        if (LN.le.LNUM) then
                            nn_in = L(LN)%mn(i2,j2)
                            if (nn_in.eq.0) cycle
                            Yc = HALF * ( L(LN)%Y(j2) + L(LN)%Y(j2+1) )
                            Yy = L(LN)%Y(j2)
                        elseif (LN.gt.LNUM) then
                            Yc = HALF * ( X(1,j2) + X(1,j2+1) )
                            Yy = X(1,j2)
                        endif
                        if (Yc.gt.Ys.and.Yc.le.Ye) then
                            if (LN.gt.LNUM) then
                                if (FSL(mn2d(i2,j2))+IWD(mn2d(i2,j2)).lt.SMALL) cycle
                            endif
                            c3 = ( Ye - Yc ) / ( Ye - Ys )  
                            c4 = ( Yc - Ys ) / ( Ye - Ys ) 
                            !Initial condition for inner layer by bilinear
                            !interpolating the data from the outermost layer
                            change = eta_in(i,j)   * c1 * c3                  &
                                   + eta_in(i+1,j) * c2 * c3                  &
                                   + eta_in(i,j+1) * c4 * c1                  &
                                   + eta_in(i+1,j+1) * c4 * c2
                            if (cflag(1:1).eq.'I') then
                                !Initial water level condition
                                if (LN.le.LNUM) then
                                    L(LN)%ETA(nn_in) = L(LN)%ETA(nn_in)+ change
                                elseif (LN.gt.LNUM) then
                                    FSL(mn2d(i2,j2)) = max(FSL(mn2d(i2,j2))   &
                                                    + change,-IWD(mn2d(i2,j2))) 
                                    call FSL_to_F(FSL(mn2d(i2,j2)),i2,j2)
                                endif
                                if (allocated(vel_in)) then
                                    DIM = ubound(vel_in,1)
                                    !Initial velocity condition
                                    if (LN.le.LNUM) then
                                        L(LN)%U(:,nn_in) =                    &
                                                vel_in(:,i,j,kss)    * c1 * c3&
                                              + vel_in(:,i+1,j,kss)  * c2 * c3&
                                              + vel_in(:,i,j+1,kss)  * c4 * c1&
                                              + vel_in(:,i+1,j+1,kss)* c4 * c2
                                    elseif (LN.gt.LNUM) then
                                        U(0:DIM-1,mn(i2,kss:kee-1,j2)) =      &
                                       vel_in(1:DIM,i,j,kss:kee-1)   * c1 * c3&
                                     + vel_in(1:DIM,i+1,j,kss:kee-1) * c2 * c3&
                                     + vel_in(1:DIM,i,j+1,kss:kee-1) * c4 * c1&
                                     + vel_in(1:DIM,i+1,j+1,kss:kee-1)* c4 * c2
                                    endif
                                endif 
                            elseif (cflag(1:1).eq.'G') then 
                                if (LN.le.LNUM) then
                                    !Change in ground level condition
                                    if (L(LN)%IC(1:1).eq.'G') then
                                        L(LN)%ZK(nn_in) = min(L(LN)%Z(L(LN)%ke),&
                                                        L(LN)%ZK(nn_in) + change)
                                    endif
                                    !Change in teibou condition
                                    do ii = 1,L(LN)%DIM
                                        ix = 0 ; iy = 0
                                        if (ii.eq.1) then
                                            if (i2.ne.iss) ix = -1
                                            c1x = ( Xe - Xx ) / ( Xe - Xs )  
                                            c2x = ( Xx - Xs ) / ( Xe - Xs )
                                            c3y = c3
                                            c4y = c4
                                        elseif (ii.eq.2) then
                                            if (j2.ne.jss) iy = -1
                                            c1x = c1
                                            c2x = c2
                                            c3y = ( Ye - Yy ) / ( Ye - Ys )  
                                            c4y = ( Yy - Ys ) / ( Ye - Ys ) 
                                        endif
                                        nn_in1 = L(LN)%mn(i2+ix,j2+iy)
                                        if (-L(LN)%H(ii,nn_in).gt.L(LN)%ZKm(nn_in).and.&
                                            -L(LN)%H(ii,nn_in).gt.L(LN)%ZKm(nn_in1)) then 
                                            change = eta_in(i,j) * c1x * c3y   &
                                                + eta_in(i+1,j) * c2x * c3y    &
                                                + eta_in(i,j+1) * c4y * c1x    &
                                                + eta_in(i+1,j+1) * c4y * c2x
                                            L(LN)%H(ii,nn_in) =                &
                                                max(-L(LN)%Z(L(LN)%ke),        &
                                                L(LN)%H(ii,nn_in) - change)
                                        elseif (-L(LN)%H(ii,nn_in).eq.L(LN)%ZKm(nn_in).or.&
                                                -L(LN)%H(ii,nn_in).eq.L(LN)%ZKm(nn_in1)) then 
                                            L(LN)%H(ii,nn_in) = -ONE * max(L(LN)%ZK(nn_in1),L(LN)%ZK(nn_in)) 
                                        else
                                            L(LN)%H(ii,nn_in) = - LI_1( L(LN)%ZK(nn_in1),L(LN)%ZK(nn_in),&
                                                                        L(LN)%DX(ii,nn_in1),L(LN)%DX(ii,nn_in) )
                                        endif
                                    enddo  
                                else
                                    ! Just do the same as IC ...
                                    FSL(mn2d(i2,j2)) = max(FSL(mn2d(i2,j2))   &
                                                        + change,-IWD(mn2d(i2,j2))) 
                                    call FSL_to_F(FSL(mn2d(i2,j2)),i2,j2)
#ifdef PRINT_ZK
                                    IWD(mn2d(i2,j2)) = IWD(mn2d(i2,j2))- change
#endif
                                endif
                            endif
                        endif
                    enddo
                endif
            enddo
        enddo
    enddo
!!$omp end parallel do
endsubroutine interp_onto_inner_layer
!
subroutine Update_ETA_with_GC(L,LN)
    use variables, only: VSMALL, ZERO, nstart, LNUM
    use type_list, only: Layer
    implicit none
    type(Layer),intent(inout) :: L
    integer,intent(in) :: LN
    integer :: nn, i , j 
!$omp parallel do private(i,j)
    do nn = 1,L%inne
        i = L%IN(1,nn) ; j = L%IN(2,nn)
        !Remove very small ground levels
        if (abs(L%ZK(nn)).lt.VSMALL) L%ZK(nn) = ZERO
        ! Make new suface by conserving depth
        if (L%coupling_dim.eq.1.or.L%F(nn).eq.ZERO.or.(nstart.eq.0.and.LN.eq.LNUM).or. &
            i.lt.L%is2.or.i.ge.L%ie2.or.j.lt.L%js2.or.j.ge.L%je2) then
            L%ETA(nn) = L%ETA(nn) + L%ZK(nn) - L%ZKm(nn)
        endif
        if (abs(L%ETA(nn)).lt.VSMALL) L%ETA(nn) = ZERO
        if (abs(L%ETA(nn)-L%ZK(nn)).lt.VSMALL) L%ETA(nn) = L%ZK(nn)  
        where (abs(L%H(:,nn)).lt.VSMALL) L%H(:,nn) = ZERO
        L%ZKm(nn) = L%ZK(nn)
        L%ETAn(nn) = L%ETA(nn)     
    enddo
!$omp end parallel do
    !Call boundary conditions for H_X,H_Y
    call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%H,L%DIM,2)
    !Boundary conditions for ZK
    call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%ZK,1,1)
    call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%ZKm,1,1)
    !Boundary condition for ETA
    if (LN.gt.1) then
        call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%ETA,1,1)   
        call equal_bound(L%inne,L%is,L%js,L%ie,L%je,L%mn,L%ETAn,1,1) 
    endif
endsubroutine Update_ETA_with_GC   
!
subroutine initial_3D
    use variables, only: LNUM, ks, ke, ZERO
    use arrays, only: L, eta_in, vel_in, a, x, nff, mn
    use interface_list, only: wdiv
    implicit none
    integer :: nn, i, j, k, nnpx, ii
    real*8  :: Qflux
    !
    if (allocated(eta_in)) deallocate(eta_in)
    if (allocated(vel_in)) deallocate(vel_in)
    allocate(eta_in(L(LNUM)%is2-1:L(LNUM)%ie2,L(LNUM)%js2-1:L(LNUM)%je2))
    allocate(vel_in(L(LNUM)%DIM+1,L(LNUM)%is2-1:L(LNUM)%ie2,L(LNUM)%js2-1:L(LNUM)%je2,ks:ke))
    do i = L(LNUM)%is2-1,L(LNUM)%ie2
        do j = L(LNUM)%js2-1,L(LNUM)%je2
            nn = L(LNUM)%MN(i,j)
            ! Transfer ETA 
            eta_in(i,j) = L(LNUM)%ETAn(nn)
            Qflux = ZERO
            do ii = 1,L(LNUM)%DIM !Loop over the number of dimensions
                if (ii.eq.1) nnpx = L(LNUM)%MN(i+1,j)
                if (ii.eq.2) nnpx = L(LNUM)%MN(i,j+1)
                !Calculating one component of the flux derivative
                Qflux = Qflux + L(LNUM)%RX(ii,nn) * ( L(LNUM)%Qn(ii,nnpx) - L(LNUM)%Qn(ii,nn) )
            enddo
            do k = ks,ke-1
                if (X(2,k+1).le.L(LNUM)%ZK(nn)) then
                    vel_in(:,i,j,k) = ZERO 
                    vel_in(L(LNUM)%DIM+1,i,j,k+1) = ZERO
                elseif (X(2,k).lt.eta_in(i,j)) then
                    ! Transfer VEL
                    vel_in(1:L(LNUM)%DIM,i,j,k) = L(LNUM)%U(1:L(LNUM)%DIM,nn)
                    ! Vertical velocity
                    vel_in(L(LNUM)%DIM+1,i,j,k+1) = wdiv( -Qflux * (X(2,k+1) - L(LNUM)%ZK(nn)) , &
                                                          L(LNUM)%ETAn(nn) - L(LNUM)%ZK(nn) )
                else
                    vel_in(1:L(LNUM)%DIM,i,j,k)   = ZERO
                    vel_in(L(LNUM)%DIM+1,i,j,k+1) = ZERO
                endif
            enddo
        enddo
    enddo
    call interp_onto_inner_layer('IC',LNUM,LNUM+1)
endsubroutine initial_3D
!
subroutine FSL_to_F(wl,i,j)
    use variables, only: ZERO, ONE, h_in, ks, ke
    use arrays, only: mn, x, f, fp, nff, fb, dx 
#ifdef NORMAL
    use arrays, only: nor, sp, spn, sfno, sfnoi
#endif
    !----------------- Input variables ----------------------------------------
    real*8,intent(in) :: wl
    integer,intent(in) :: i, j
    !----------------- Temporary variables ------------------------------------
    integer :: k, kk, nn_1, nn_s
    real*8 :: bot
!=============================================================================!
!   Converts the calculated water level to F fraction                         !
!=============================================================================!
!
    do k = ks,ke-1
        nn_1 = mn(i,k,j)
        if (nff(nn_1)%f.eq.-1) cycle
        if (x(2,k+1).le.wl+h_in) then
            !Fluid cell
            nff(nn_1)%f = 1   ; nff(nn_1)%b = 0  
            f(nn_1)     = ONE ; fp(nn_1)    = ONE
        elseif (x(2,k).ge.wl+h_in) then
            !Air cell
            nff(nn_1)%f = 0    ; nff(nn_1)%b = 0  
            f(nn_1)     = ZERO ; fp(nn_1)    = ZERO
        else
            bot = x(2,k+1) - fb(nn_1) * dx(2,k)
            f(nn_1)  = (wl+h_in-bot) / dx(2,k) / fb(nn_1)
            fp(nn_1) = (wl+h_in-bot) / dx(2,k) / fb(nn_1) 
            !Surface cell
            if (f(nn_1).le.ZERO) cycle
            nff(nn_1)%f = 2 ; nff(nn_1)%b = -3 
#ifdef NORMAL
            ! Adjust normal information
            if (sfno(nn_1).eq.0) then
                nn_s = 0
                do kk = ks,ke-1
                    if (nff(mn(i,kk,j))%f.eq.-1) cycle
                    if (sfno(mn(i,kk,j)).ne.0) then
                        nn_s = sfno(mn(i,kk,j))
                        sfno(mn(i,kk,j)) = 0
                        exit
                    endif
                enddo
                if (nn_s.ne.0) then
                    sfno(nn_1) = nn_s ; sfnoi(nn_s) = nn_1
                    spn(nn_1,0:1) = spn(mn(i,kk,j),0:1)
                    sp(nn_1,0:1)  = spn(nn_1,0:1)
                endif
            endif
            sp (nn_1,2) = wl + h_in ; spn(nn_1,2) = sp(nn_1,2)
            nor(nn_1,0:2) = [ ZERO, ZERO, ONE ]
            nor(nn_1,3) = - sp(nn_1,2)
#endif
        endif
    enddo
endsubroutine FSL_to_F
!------------------------------------------------------------------------------
#ifdef PRINT_ZK
subroutine print_IWD
    use arrays, only: L, IWD, mn2d
    use variables, only: LNUM, is, ie, js, je
    integer :: u, i, j
!=============================================================================!
!   PRINTS THE NEW WATER DEPTHS OUT TO FILE                                   !
!=============================================================================!
!
    open(newunit=u,file='IWDout_3D.dat',status='new',action='write')
        do j = je-1, js,-1
            write(u,*) ( IWD(mn2d(i,j)), i = is, ie-1 )
        enddo
    close(u)
    if (LNUM.gt.0) then
        open(newunit=u,file='ZKout_2D.dat',status='new',action='write')
            do j = L(LNUM)%je-1, L(LNUM)%js,-1
                write(u,'(9999(F12.6))') ( L(LNUM)%ZK(L(LNUM)%mn(i,j)),       &  
                                           i = L(LNUM)%is, L(LNUM)%ie-1 )
            enddo
        close(u)
    endif
endsubroutine print_IWD
#endif