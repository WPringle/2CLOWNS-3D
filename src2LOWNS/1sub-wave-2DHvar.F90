!%%%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-wave-2DHvar.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutines for converting 3D quantities to 2DH ones
subroutine calc_IWD
    use variables, only: is, ie, js, je, ks, ke, inn2d, h_in
    use arrays, only: x, dx, nff, mn, fb, in2d, IWD
    implicit none
    integer :: nn2d, i, j, k
    !Calculate the initial mean water depth
    do nn2D = 1, inn2d
        i = in2d(0,nn2d) ; j = in2d(1,nn2d)
        do k = ks,ke-1
            if (nff(mn(i,k,j))%f.eq.-1) cycle
            !Subtracts the ground level from the initial free surface level
            IWD(nn2D) = h_in - x(2,k+1) + fb(mn(i,k,j)) * dx(2,k)
            exit
        enddo
    enddo
endsubroutine calc_IWD
!------------------------------------------------------------------------------
subroutine calc_FSL
    use variables, only: inn2d, ks, ke, ZERO, ONE, VSMALL, h_in
    use arrays, only: in2d, nff, mn, f, x, dx, sfno, sp, fb, a,               &
                      IWD, FSL, FSLmax
    implicit none
    integer :: nn2d, i, j, k, nn
    real*8 :: dzz
    !Lets allocate the IWD and FSL matrices..
    if (.not.allocated(IWD)) then
        allocate(IWD(inn2d),FSL(inn2d),FSLmax(inn2d))
        IWD = ZERO ; FSLmax = -1d4 
    endif
    !Calculates the current free surface water level
!$omp parallel do private (i,j,k,dzz,nn)
    do nn2D = 1, inn2d
        i = in2d(0,nn2d) ; j = in2d(1,nn2d) ; FSL(nn2D) = -1d4 
        do k = ke-1,ks,-1
            nn = mn(i,k,j)
            if ( nff(nn)%f.eq.0.or.( nff(nn)%f.eq.-1.and.                     &
                ( nff(nn)%b.eq.0.or.nff(nn)%b.eq.10) ) ) cycle
            if (nff(nn)%f.eq.-1.and.f(nn).gt.ZERO.and.f(nn).lt.ONE)  then
                !At the boundary
                FSL(nn2D) = x(2,k) + f(nn) * dx(2,k)
                exit
            elseif (nff(nn)%f.ge.1) then             
#ifdef NORMAL       !At free surface cell with normal known
                if (sfno(nn).ne.0) then
                    FSL(nn2D) = sp(nn,2); exit
                endif
#endif
                !At free surface cell with no-normal or at fluid cell
                if (a(2,nn).eq.ZERO) then
                    !Water level above some ground
                    dzz = dx(2,k) * ( fb(nn) * ( f(nn) - ONE ) + ONE )
                else
                    !vertical structure or none present
                    dzz = f(nn) * fb(nn) * dx(2,k)
                endif
                FSL(nn2D) = x(2,k) + dzz
                exit
            endif
        enddo
        if (FSL(nn2D).eq.-1d4) then ! No water in cell column: 
            do k = ks,ke-1          ! Set water level equal to bed level
                if (nff(mn(i,k,j))%f.eq.-1) cycle
                !Gets the ground level 
                FSL(nn2D) = x(2,k+1) - fb(mn(i,k,j)) * dx(2,k)
                exit
            enddo
        endif
        FSL(nn2D) = FSL(nn2D) - h_in
        if (abs(FSL(nn2D)).lt.VSMALL) FSL(nn2D) = ZERO
        FSLmax(nn2D) = max(FSL(nn2D),FSLmax(nn2D))
    enddo
!$omp end parallel do
endsubroutine calc_FSL
!------------------------------------------------------------------------------
subroutine calc_FLUX
    use variables, only: t, ks, ke, inn2d, ZERO, ONE, VSMALL
    use arrays, only: dim3d, nff, dx, u, a, in2d, mn, f, Flux_3D, U_3D, U_3Dmax
    use interface_list, only: wdiv, LI_1
    implicit none
    integer :: nn2D, i, j, k, ii, ix, iy, nn, nnmx, DIM
    real*8  :: dxm, dxp, dp, Depth
    ! Calculates the depth-averaged velocities and fluxes
    if (.not.dim3d) DIM = 0
    if (dim3d)      DIM = 1
    if (.not.allocated(Flux_3D)) then
        allocate(Flux_3D(0:DIM,0:inn2D),U_3D(0:DIM,0:inn2D),                  &
                 U_3Dmax(0:DIM,0:inn2D))
        Flux_3D(:,0) = ZERO ; U_3D(:,0) = ZERO ; U_3Dmax = ZERO
        if (t.gt.ZERO) call read_max_data_MAT
    endif
!$omp parallel do private (ii,i,j,k,nn,nnmx,ix,iy,dxm,dxp,dp,Depth) 
    do nn2D = 1,inn2D
        i = in2d(0,nn2D) ; j = in2d(1,nn2D)
        do ii = 0,DIM
            if (ii.eq.0) then
                ix = -1; iy = 0; dxm = dx(ii,i+ix); dxp = dx(ii,i)
            elseif (ii.eq.1) then
                iy = -1; ix = 0; dxm = dx(ii,j+iy); dxp = dx(ii,j)
            endif
            Flux_3D(ii,nn2D) = ZERO; Depth = ZERO
            do k = ks,ke-1
                if (nff(mn(i,k,j))%f.eq.-1.and.                               &
                    nff(mn(i+ix,k,j+iy))%f.eq.-1) cycle
                if (nff(mn(i,k,j))%f.eq. 0.and.                               &
                    nff(mn(i+ix,k,j+iy))%f.eq. 0) exit
                nn = mn(i,k,j) ; nnmx = mn(i+ix,k,j+iy)
                dp = a(ii,nn) * dx(2,k)                                       &
                    * LI_1(min(ONE,f(nnmx)),min(ONE,f(nn)),dxm,dxp)
                Flux_3D(ii,nn2D) = Flux_3D(ii,nn2D) + U(ii,nn) * dp
                Depth = Depth + dp
            enddo
            U_3D(ii,nn2D) = wdiv(Flux_3D(ii,nn2D),Depth)
            if (abs(Flux_3D(ii,nn2D)).lt.VSMALL) Flux_3D(ii,nn2D) = ZERO
            if (abs(U_3D(ii,nn2D)).lt.VSMALL) U_3D(ii,nn2D) = ZERO
            U_3Dmax(ii,nn2D) = max(U_3D(ii,nn2D),U_3Dmax(ii,nn2D))
        enddo
    enddo
!$omp end parallel do
endsubroutine calc_FLUX