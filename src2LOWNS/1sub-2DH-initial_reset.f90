!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&          INIT2D : INITIALISES THE ARRAYS AND PARAMETERS FOR EACH LAYER  &&!
!&&        WRITTEN BY PRINGLE (2014) & UPDATED BY PRINGLE IN JUNE 2015      &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
subroutine init2D
    use variables, only: LNUM, t, g, HCIRC, R, PSI, nstart,                   &
                         ZERO, TEN, ONE, PI_8, TWO, VSMALL, HALF, TWOTHIRD  
    use arrays, only: L
    use interface_list, only: LI_1, round, wdiv, find_layer_interp
    implicit none
    integer :: LN, nn, i, j, nnmx, nnmy, nnpy, nnpx, ii, jj, cmaxn, dec
    real*8  :: Htemp
!
Layer_loop: DO LN = 1,LNUM !Loop over all the layers
!
!------- Convert the read data into actual dimensional quantities -------------
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                Mannings n setting                           !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    if (L(LN)%Mandef.ge.ZERO) then
        !If constant mannings n defined, set here
        L(LN)%MAN = L(LN)%Mandef
    else
        ! Read Mannings n from F data
        ! Boundary conditions for Mannings n
        call equal_bound(L(LN)%inne,L(LN)%is,L(LN)%js,                        &
                         L(LN)%ie,L(LN)%je,L(LN)%mn,L(LN)%MAN,1,1)
    endif
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                      Spherical Coordinate Parameters                        !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !If using Spherical Coordinates, create the latitude multipliers:
    if (L(LN)%GEOC.eq.1) then
        do j = 1,L(LN)%je+1 
            L(LN)%sinY(j) = PSI * sin( L(LN)%Y(j) * PI_8 / HCIRC)
            L(LN)%secY(j) = ONE / cos( L(LN)%Y(j) * PI_8 / HCIRC)
        enddo
    endif
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                      Cell size arrays                                       !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !Find the domain height
    L(LN)%DZ = L(LN)%Z(L(LN)%KE) - L(LN)%Z(L(LN)%KS)
    ! Determine decimal place for rounding. Determined by suggesting that
    ! dx = 10,000 requires no decimal places (is rounded to zero dec) 
    ! dx = 1,000.0 is rounded to one decimal place and so forth
    do dec = 0,10
        !Guess from first cell size
        if ((L(LN)%X(L(LN)%is+1)-L(LN)%X(L(LN)%is))/TEN**(5-dec).gt.ONE) exit
    enddo
    L(LN)%DEC = dec; L(LN)%DX(:,0) = ZERO
    do nn = 1,L(LN)%inne
        i = L(LN)%IN(1,nn) ; j = L(LN)%IN(2,nn)
        !
        !Determine the cell sizes
        L(LN)%DX(1,nn) = round( L(LN)%X(i+1) - L(LN)%X(i) , L(LN)%DEC ) 
        L(LN)%DX(1,0) = max(L(LN)%DX(1,0),L(LN)%DX(1,nn))
        if (L(LN)%DIM.eq.2) then
            L(LN)%DX(2,nn) = round(L(LN)%Y(j+1)-L(LN)%Y(j),dec)
            L(LN)%DX(2,0) = max(L(LN)%DX(2,0),L(LN)%DX(2,nn))
        endif
        !
        !Conversion to radians in spherical coordinates
        if (L(LN)%GEOC.eq.1) L(LN)%DX(:,nn) = L(LN)%DX(:,nn) * PI_8 / HCIRC
        !
        !Find inversion of DX for faster calculation
        if (L(LN)%DX(1,nn).gt.ZERO) then
            nnmx = L(LN)%mn(i-1,j)
            L(LN)%RX(1,nn) = ONE / L(LN)%DX(1,nn)
            if (L(LN)%DX(1,nnmx).gt.ZERO) then
                L(LN)%RXc(1,nn) = TWO / ( L(LN)%DX(1,nn) + L(LN)%DX(1,nnmx) )
            else
                L(LN)%RXc(1,nn) = L(LN)%RX(1,nn)
            endif
        endif
        if (L(LN)%DIM.eq.2) then
            nnmy = L(LN)%mn(i,j-1)
            if (L(LN)%DX(2,nn).gt.ZERO) then
                L(LN)%RX(2,nn) = ONE / L(LN)%DX(2,nn)
                if (L(LN)%DX(2,nnmy).gt.ZERO) then
                    L(LN)%RXc(2,nn) = TWO / (L(LN)%DX(2,nn) + L(LN)%DX(2,nnmy))
                else
                    L(LN)%RXc(2,nn) = L(LN)%RX(2,nn)
                endif
            endif
        endif
        !In spherical coordinates scale with the earth's radius
        if (L(LN)%GEOC.eq.1) then
            j = L(LN)%IN(2,nn)
            L(LN)%RX(1,nn) = L(LN)%RX(1,nn)  * R                              &
                           * HALF * (L(LN)%secY(j) + L(LN)%secY(j+1))
            L(LN)%RX(2,nn) = L(LN)%RX(2,nn)  * R
            L(LN)%RXc(1,nn)= L(LN)%RXc(1,nn) * R * L(LN)%secY(j)
            L(LN)%RXc(2,nn)= L(LN)%RXc(2,nn) * R
        endif
    enddo
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!               Nesting interpolation data                                    !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !Find the integer bounds where the next inner layer is situated
    call find_layer_interp(L(LN),L(min(LN+1,LNUM)),LN)
    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                    Ground level calculation                                 !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !
    !Convert ground level into dimensional quantity
    do nn = 1,L(LN)%inne
        !Calculate the height of the ground ZK using the FB value
        L(LN)%ZK(nn) = L(LN)%DZ * ( ONE - L(LN)%ZK(nn)) + L(LN)%Z(L(LN)%KS)
    enddo
    where (abs(L(LN)%ZK).lt.VSMALL) L(LN)%ZK = ZERO
    !Boundary conditions for ZKm
    call equal_bound(L(LN)%inne,L(LN)%is,L(LN)%js,&
                     L(LN)%ie,L(LN)%je,L(LN)%mn,L(LN)%ZK,1,1)
    !Set old ground elevation equal to current one
    L(LN)%ZKm = L(LN)%ZK
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                   Initial water depth calculation                           !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !Convert initial water depths into dimensional quantity
    do nn = 1,L(LN)%inne
        i = L(LN)%IN(1,nn) ; j = L(LN)%IN(2,nn)
        !Find the initial water depths on the cell boundaries
        do ii = 1,L(LN)%DIM
            if (i.gt.L(LN)%IS-1.and.j.gt.L(LN)%JS-1) then
                if (ii.eq.1) nnmx = L(LN)%MN(i-1,j)
                if (ii.eq.2) nnmx = L(LN)%MN(i,j-1)
                Htemp = L(LN)%DZ * (ONE - L(LN)%H(ii,nn) )                    &
                      + L(LN)%Z(L(LN)%KS)  
                if (L(LN)%H(ii,nn).eq.ONE.or.                                 &
                    Htemp.le.L(LN)%ZK(nn).or.Htemp.le.L(LN)%ZK(nnmx)) then
                    ! Add new condition to treat a slope > 1:1.5 as a wall...
                    !if (abs(L(LN)%ZK(nnmx)-L(LN)%ZK(nn)) *                    &
                    !        L(LN)%RXc(ii,nn) > TWOTHIRD) then
      				if (LN.eq.5.and.ii.eq.1.and.L(LN)%X(i).gt.408100d0.and.L(LN)%X(i).lt.408200d0.and.&
                        L(LN)%Y(j).ge.845010d0.and.L(LN)%Y(j).lt.845320d0.and.&
                        (L(LN)%ZK(nn).eq.-19.86d0.or.L(LN)%ZK(nnmx).eq.-19.86d0)) then
                       L(LN)%H(ii,nn) = -ONE * max(L(LN)%ZK(nnmx),L(LN)%ZK(nn)) 
                        !if (L(LN)%ZK(nnmx) * L(LN)%ZK(nn).lt.ZERO) then
                        !    ! Treat as a wall just for overtopping flow
                        !    ! by adding a small amount so that it enters 
                        !    ! overtopping formula
                        !    L(LN)%H(ii,nn) = L(LN)%H(ii,nn) - L(LN)%Dmin 
                        !endif
                    else
                        ! No teibou or teibou is lower than ground level
                         L(LN)%H(ii,nn) = - LI_1( L(LN)%ZK(nnmx),L(LN)%ZK(nn),&
                                            L(LN)%DX(ii,nnmx),L(LN)%DX(ii,nn) )
                    endif
                else 
                    ! Consider the cell boundary height using teibou data
                    L(LN)%H(ii,nn) = - Htemp
                endif
            endif
        enddo
    enddo    
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                      Change in ground level Calculation                     !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
    !Get the change in ground level if defined as inflow condition
    if (t.eq.ZERO) then
        call calc_IC_or_GC('GC',0,LN)
        call calc_IC_or_GC('GT',0,LN)
    endif
    !Remove very small ground levels
    where (abs(L(LN)%ZK).lt.VSMALL)  L(LN)%ZK = ZERO
    !Boundary conditions for ZK
    call equal_bound(L(LN)%inne,L(LN)%is,L(LN)%js,&
                     L(LN)%ie,L(LN)%je,L(LN)%mn,L(LN)%ZK,1,1)
    !Remove very small depths
    where (abs(L(LN)%H).lt.VSMALL) L(LN)%H = ZERO
    !Initialise wall for 2D outer layer (change initial water depth)
    if (LN.eq.1) call make_wall(L(LN))
    !Call boundary conditions for H_X,H_Y
    call equal_bound(L(LN)%inne,L(LN)%is,L(LN)%js,L(LN)%ie,&
                     L(LN)%je,L(LN)%mn,L(LN)%H,L(LN)%DIM,2)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                      Initial water level setting                            !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
    if (t.eq.ZERO) then
        !Alter ETA based on tide (add more volume through F value)
        call alter_tide_and_remove_water_on_land(L(LN))
    endif
    L(LN)%INDISP = 0; L(LN)%INWET = 0
    do nn = 1,L(LN)%inne
        !Calculate ETA from F
        L(LN)%ETA(nn) = L(LN)%ZK(nn)                                          &
                      + L(LN)%F(nn) * (L(LN)%Z(L(LN)%KE) - L(LN)%ZKm(nn)) 
        if (t.gt.ZERO) then
            if (L(LN)%F(nn).eq.ZERO) then
                L(LN)%ETAn(nn) = L(LN)%ETA(nn)
            endif
        endif
        i = L(LN)%IN(1,nn) ; j = L(LN)%IN(2,nn)
        if (i.lt.L(LN)%is.or.i.ge.L(LN)%ie.or.                                &
            j.lt.L(LN)%js.or.j.ge.L(LN)%je) cycle
        ! Count number of cells that are initially wet or are a shoreline cell
        nnmx = L(LN)%MN(i-1,j) ; nnpx = L(LN)%MN(i+1,j)
        nnmy = L(LN)%MN(i,j-1) ; nnpy = L(LN)%MN(i,j+1)
        if (t.gt.ZERO) then
            if ( L(LN)%NF(i,L(LN)%ks,j).ne.-1.and.(L(LN)%coupling_dim.eq.1    &
                 .or.(LN.eq.LNUM.and.nstart.eq.0)                             &   
                 .or.(L(LN)%DIM.eq.1.and.(i.le.L(LN)%is2.or.i.ge.L(LN)%ie2-1))& 
                 .or.(L(LN)%DIM.eq.2.and.(i.le.L(LN)%is2.or.i.ge.L(LN)%ie2-1  &
                 .or.j.le.L(LN)%js2.or.j.ge.L(LN)%je2-1))) ) then 
                if (L(LN)%ETAn(nn)-L(LN)%ZK(nn).gt.L(LN)%Dmin.or.             &
                    L(LN)%ETAn(nnmx)-L(LN)%ZK(nnmx).gt.L(LN)%Dmin.or.         &
                    L(LN)%ETAn(nnmy)-L(LN)%ZK(nnmy).gt.L(LN)%Dmin.or.         &
                    L(LN)%ETAn(nnpx)-L(LN)%ZK(nnpx).gt.L(LN)%Dmin.or.         &
                    L(LN)%ETAn(nnpy)-L(LN)%ZK(nnpy).gt.L(LN)%Dmin) then
                    L(LN)%INWET    = L(LN)%INWET + 1
                    L(LN)%INW(:,L(LN)%INWET) = [ i, j ]
                endif
            endif
        endif
        if (L(LN)%coupling_dim.eq.2.and.i.ge.L(LN)%is2.and.                   &
            i.lt.L(LN)%ie2.and.j.ge.L(LN)%js2.and.j.lt.L(LN)%je2) then
            if (L(LN)%mn2(0,i,j).eq.0) cycle
        endif 
        ! Count number of cells that have positive h and not in inner layer
        if (L(LN)%ZK(nn).lt.ZERO) then
            L(LN)%INDISP   = L(LN)%INDISP + 1
            L(LN)%MND(i,j) = L(LN)%INDISP 
            L(LN)%IND(:,L(LN)%INDISP) = [ i, j ]
        endif
    enddo
    !In the case when we calculate from t = 0
    if (t.eq.ZERO) then
        !Get the initial free surface change 
        !and velocity for initial condition if reqd.
        call calc_IC_or_GC('IC',0,LN)
        call calc_IC_or_GC('IT',0,LN)
        ! Alter ETA if we have depths smaller than allowable 
        ! initial one and recount wet cells
        do nn = 1,L(LN)%inne 
            i = L(LN)%IN(1,nn) ; j = L(LN)%IN(2,nn)
            if (i.lt.L(LN)%is.or.i.ge.L(LN)%ie.or.                            &
                j.lt.L(LN)%js.or.j.ge.L(LN)%je) cycle
            ! Count number of cells that are initially 
            ! wet or are a shoreline cell
            nnmx = L(LN)%MN(i-1,j) ; nnpx = L(LN)%MN(i+1,j)
            nnmy = L(LN)%MN(i,j-1) ; nnpy = L(LN)%MN(i,j+1)
            if (L(LN)%ETA(nn) - L(LN)%ZK(nn).le.L(LN)%Dmin) then
                L(LN)%ETA(nn) = L(LN)%ZK(nn)
            elseif (L(LN)%ETA(nn)-L(LN)%ZK(nn).gt.L(LN)%Dmin.or.              &
                L(LN)%ETA(nnmx)-L(LN)%ZK(nnmx).gt.L(LN)%Dmin.or.              &
                L(LN)%ETA(nnmy)-L(LN)%ZK(nnmy).gt.L(LN)%Dmin.or.              &
                L(LN)%ETA(nnpx)-L(LN)%ZK(nnpx).gt.L(LN)%Dmin.or.              &
                L(LN)%ETA(nnpy)-L(LN)%ZK(nnpy).gt.L(LN)%Dmin) then
                if ( L(LN)%NF(i,L(LN)%ks,j).ne.-1.and.(L(LN)%coupling_dim.eq.1&
                     .or.(LN.eq.LNUM.and.nstart.eq.0)                         &    
                     .or.(L(LN)%DIM.eq.1.and.                                 &
                         (i.le.L(LN)%is2+1.or.i.ge.L(LN)%ie2-2))              & 
                     .or.(L(LN)%DIM.eq.2.and.                                 &
                         (i.le.L(LN)%is2+1.or.i.ge.L(LN)%ie2-2                &
                     .or.j.le.L(LN)%js2+1.or.j.ge.L(LN)%je2-2))) ) then 
                    L(LN)%INWET    = L(LN)%INWET + 1
                    L(LN)%INW(:,L(LN)%INWET) = [ i, j ]      
                endif
            endif
        enddo
        !Alter ETA if we have very small values
        where (abs(L(LN)%ETA).lt.VSMALL) L(LN)%ETA = ZERO
        !Boundary condition for ETA
        call equal_bound(L(LN)%inne,L(LN)%is,L(LN)%js,L(LN)%ie,               &
                         L(LN)%je,L(LN)%mn,L(LN)%ETA,1,1)
        L(LN)%ETAn = L(LN)%ETA
    endif
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                   Calculating DT from Courant number                        !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
    call set_dt(L(LN))
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                  Calculate the water depths on the boundary                 !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    L(LN)%D = ZERO ; L(LN)%Dn = ZERO
    call calc_depth(L(LN),LN)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                Dispersion parameters for implicit or explicit scheme        !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
    call disp_param(L(LN),0)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                         Calculating U,V into M,N                            !                         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   
!-------- Initialise the arrays that are not read from file -------------------
    L(LN)%Q = ZERO ; L(LN)%Qn = ZERO ; L(LN)%Qs = ZERO
    do nn = 1,L(LN)%inne
        if (t.gt.ZERO) L(LN)%ETA(nn) = L(LN)%ETAn(nn)
        i = L(LN)%IN(1,nn) ; j=L(LN)%IN(2,nn)
        if (i.lt.L(LN)%is.or.i.gt.L(LN)%ie.or.&
            j.lt.L(LN)%js.or.j.gt.L(LN)%je) cycle
        do ii = 1,L(LN)%DIM
            L(LN)%Qn(ii,nn) = L(LN)%U(ii,nn) * L(LN)%D(ii,nn)
        enddo
    enddo
    !Call boundary conditions for M,N
    call equal_bound(L(LN)%inne,L(LN)%is,L(LN)%js,L(LN)%ie,                   &
                     L(LN)%je,L(LN)%mn,L(LN)%Qn,L(LN)%DIM,2)
    L(LN)%Q = L(LN)%Qn
!
ENDDO Layer_loop
!
endsubroutine init2D
    
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!            CALC_FUV : Calculates F & U from ETA & Q for visualization       !
!                       and resets ETA & Q                    PRINGLE (2015)  !
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
subroutine calc_FUV(L,LN)
    use variables, only: ZERO, ONE, HALF, LNUM, nstart
    use arrays, only: ix
    use type_list, only: Layer
    use interface_list, only: wdiv, LI_1
    use omp_lib
    implicit none
!----------------- Input/Output variables -------------------------------------
    type(Layer),intent(inout) :: L
    integer,intent(in) :: LN
    integer :: nn,i,j,ii,jj,nnpx,nnmx
    integer :: mythr, nthr, nn_s, nn_e
    integer, allocatable :: inwet_thr(:),inw_thr(:,:)
    real*8 :: Depth, Vel
! Manual parallel loop: count number of wet cells locally in each 
!                       thread and then collect into global array
nthr = 1 ; mythr = 0
!$omp parallel &
!$omp private(nn,i,j,ii,Vel,depth,nnpx,nnmx,mythr,nn_s,nn_e,inw_thr,jj)
!$omp single
!$ nthr = omp_get_num_threads()
allocate(inwet_thr(-1:nthr-1)) ; inwet_thr = 0
!$omp end single
!$ mythr = omp_get_thread_num()
allocate(inw_thr(2, L%inne / nthr)) 
nn_s  = int(real(mythr * L%inne / nthr))  + 1
nn_e  = int(real(((mythr + 1) * L%inne / nthr))) 
FUV_LOOP: do nn = nn_s,nn_e
    i = L%IN(1,nn) ; j = L%IN(2,nn)
    !
    do ii = 1,L%DIM
        !Reset Q
        L%Q(ii,nn) = L%Qn(ii,nn)
        !Calc velocities at n
        L%U(ii,nn) = wdiv( L%Q(ii,nn) , L%D(ii,nn) )
    enddo
    !Calc F, NF at n+1/2 and recount number of wet cells etc.
    Depth = L%ETAn(nn) - L%ZK(nn)
    if (Depth.gt.L%Dmin) then
        L%F(nn) = min( Depth / (L%Z(L%KE)-L%ZK(nn)) , ONE )
        if (L%NF(i,L%KS,j).ne.-1) then
            ! Add a wet cell
            L%NF(i,L%KS,j) = 2; L%NFB(i,L%KS,j) = -3
            if (L%coupling_dim.eq.1.or.(LN.eq.LNUM.and.nstart.eq.0).or.       &
                (L%DIM.eq.1.and.(i.lt.L%is2+2.or.i.ge.L%ie2-2)).or.           &
                (L%DIM.eq.2.and.(i.lt.L%is2+2.or.i.ge.L%ie2-2.or.             &
                j.lt.L%js2+2.or.j.ge.L%je2-2))) then 
                inwet_thr(mythr) = inwet_thr(mythr) + 1
                inw_thr(:,inwet_thr(mythr)) = [ i, j ]   
            endif
        endif
        !Calc Max Quantities
        L%ETAmax(nn) = max(L%ETAn(nn),L%ETAmax(nn))
        L%Umax(:,nn) = max(abs(L%U(:,nn)),L%Umax(:,nn))
    else
        L%F(nn) = ZERO
        if (L%NF(i,L%KS,j).ne.-1) then
            L%NF(i,L%KS,j) = 0
            do ii = 1,L%DIM
                if (ii.eq.1) then
                    nnpx = L%mn(i+1,j) ; nnmx = L%mn(i-1,j)
                elseif (ii.eq.2) then
                    nnpx = L%mn(i,j+1) ; nnmx = L%mn(i,j-1)
                endif
                if (L%ETAn(nnpx)-L%ZK(nnpx).gt.L%Dmin.or.                     &
                    L%ETAn(nnmx)-L%ZK(nnmx).gt.L%Dmin) then
                     ! Add a shoreline cell in direction of fluid
                    if (L%ETAn(nnpx)-L%ZK(nnpx).gt.L%Dmin) then
                        L%NFB(i,L%KS,j) = ii
                    elseif (L%ETAn(nnmx)-L%ZK(nnmx).gt.L%Dmin) then
                        L%NFB(i,L%KS,j) = -ii
                    endif
                    if (L%coupling_dim.eq.1.or.(LN.eq.LNUM.and.nstart.eq.0).or.&
                       (L%DIM.eq.1.and.(i.lt.L%is2+2.or.i.ge.L%ie2-2)).or.    &
                       (L%DIM.eq.2.and.(i.lt.L%is2+2.or.i.ge.L%ie2-2.or.      &
                        j.lt.L%js2+2.or.j.ge.L%je2-2))) then 
                        inwet_thr(mythr) = inwet_thr(mythr) + 1
                        inw_thr(:,inwet_thr(mythr)) = [ i, j ]   
                    endif
                    exit
                endif
            enddo
        endif
    endif
    !Reset ETA
    L%ETA(nn) = L%ETAn(nn)
enddo FUV_LOOP
!$omp barrier
!Put into global array
if (mythr.eq.0) L%INWET = sum(inwet_thr)
do nn = 1,inwet_thr(mythr)
    L%INW(:,sum(inwet_thr(-1:mythr-1))+nn) = inw_thr(:,nn)
enddo
deallocate(inw_thr)
!$omp end parallel
!
endsubroutine calc_FUV
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!   DISP_PARAM :Determines the coefficients required for either the explicit  !
!               or implicit dispersion correction schemes      PRINGLE (2015) !
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!    
subroutine disp_param(L,iflag)  
    use variables, only: A, B, R, ZERO, ONE, TWO, g
    use type_list, only: Layer 
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    integer,intent(in) :: iflag
    !----------------- Temporary variables ------------------------------------
    integer :: nn,nnw,i,j,nnmx,nnpx,nnwm,nnwp,ii,jj
!
    if (L%DISP.eq.0) return !Return if no dispersion required
    if (L%DISP.eq.2) allocate(L%AL(L%INDISP,4),L%LL(L%INDISP,4),              &
                              L%DIAG(L%INDISP),L%BREAK(L%INDISP))
    if (L%DISP.eq.1) allocate(L%GMA(L%DIM,L%INDISP),L%APA(L%DIM,L%INDISP))
!$omp parallel do private (i,j,ii,jj,nn,nnmx,nnpx,nnwm,nnwp)
DISP_LOOP: do nnw = 1 , L%INDISP
    i = L%IND(1,nnw) ; j = L%IND(2,nnw)
    nn = L%MN(i,j)
    
    DISP_IF: if (L%DISP.eq.1) then
!       Calculate the frequency dispersion optimizing parameters
!       For explicit correction scheme:
        ! See:                                                                !
        ! Cho, Sohn, & Lee (2007). Practical modified scheme of linear        !
        ! shallow-water equations for distant propagation of tsunamis.        !
        ! Ocean Engineering, 34, 1769-1777                                    !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L%GMA(:,nn)=ZERO ; L%APA(:,nn)=ZERO
        do ii = 1,L%DIM
            if (((ii.eq.1.and.i.ge.L%IS.and.i.le.L%IE-1).or.                  &
                (ii.eq.2.and.j.ge.L%JS.and.j.le.L%JE-1)).and.                 &
                    L%H(ii,nn).gt.ZERO) then 
                if (ii.eq.1) jj = L%DIM
                if (ii.eq.2) jj = 1
                L%APA(ii,nn) = min( ( L%H(ii,nn) *  (TWO * TWO * L%H(ii,nn) + &
                                     g * L%dt**2 ) - ONE / L%RX(ii,nn)**2 )   &
                                   * L%RX(ii,nn)**2 , ONE)
                L%GMA(ii,nn) = L%APA(ii,nn) + ONE
            endif
        enddo        
    elseif (L%DISP.eq.2) then
!       Calculate the linear dispersion matrix (constant coefficients)
!       For implicit correction scheme:
        ! See:                                                                !
        ! Shigihara, Fujima (2007). Adequate numerical scheme for dispersive  !
        ! wave theory for tsunami simulation and development of  new numerical! 
        ! algorithm. JSCE Coastal and Environmental Engineering, 63 (1), 51-66!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L%BREAK(nnw) = 0 
        L%LL(nnw,:)  = 0 ; L%AL(nnw,:)  = ZERO ; L%DIAG(nnw) = - ONE
        jj = 0
        do ii = 1,L%DIM
            if (ii.eq.1) then
                !Cycle outside dispersion calculation domain
                if (i.eq.L%is.or.i.eq.L%ie-1) cycle
                nnmx = L%mn(i-1,j)  ; nnpx = L%mn(i+1,j)
                nnwm = L%mnd(i-1,j) ; nnwp = L%mnd(i+1,j)
            elseif (ii.eq.2) then
                !Cycle outside dispersion calculation domain
                if (j.eq.L%js.or.j.eq.L%je-1) cycle
                nnmx = L%mn(i,j-1)  ; nnpx = L%mn(i,j+1)
                nnwm = L%mnd(i,j-1) ; nnwp = L%mnd(i,j+1)
            endif
            ! Find the matrix mapping and coefficient entries
            if (nnwm.ne.0) then
                jj = jj + 1
                L%LL(nnw,jj) = nnwm
                L%AL(nnw,jj) = A * TWO * L%ZK(nn) * L%ZK(nn) / L%DX(ii,nn)    &
                             / ( L%DX(ii,nn) + L%DX(ii,nnmx) ) 
            endif
            if (nnwp.ne.0) then
                jj = jj + 1
                L%LL(nnw,jj) = nnwp
                L%AL(nnw,jj) = A * TWO * L%ZK(nn) * L%ZK(nn) / L%DX(ii,nn)    &
                             / ( L%DX(ii,nn) + L%DX(ii,nnpx) )                  
            endif
        enddo
        L%DIAG(nnw) = L%DIAG(nnw) - sum(L%AL(nnw,:))
        !Adjust for spherical coordinates
        !if (L%GEOC.eq.1) then
        !     L%DIAG(nn) = L%DIAG(nn) * R * R
        !     L%AL(:,nn) = L%AL(:,nn) * R * R
        !endif
    endif DISP_IF
enddo DISP_LOOP
!$omp end parallel do
endsubroutine disp_param
!    
subroutine set_dt(L)
    use variables, only: ZERO, g, cl_set   
    use type_list, only: Layer 
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    !----------------- Temporary variables ------------------------------------
    integer :: nn
    real*8 :: cmax
    !
    ! Do not need to calc dt if not equal to zero
    if (L%dt.ne.ZERO) return 
    ! Decide dt from courant
    nn   = minloc(L%ZK,1)
    if (L%ZK(nn).ge.ZERO) then
        stop 'Max Initial water depth <= 0; Please set dt manually'
    endif
    cmax = sqrt( g * - L%ZK(nn) ) * maxval( L%RX(:,nn) )
    L%dt = L%cl_set / cmax
!
endsubroutine set_dt
!
subroutine check_courant
    use arrays, only: L
    use variables, only: g, LNUM, ZERO
    implicit none
    !----------------- Temporary variables ------------------------------------
    integer :: nn, i, j, LN, nnw
    real*8 :: Cranx, CRAN2D, VMAX2D
    !
    DO LN = 1,LNUM !Loop over all layers
        VMAX2D = ZERO ; CRAN2D = ZERO
!$omp parallel do private(i,j,nn,Cranx) &
!$omp reduction(max:VMAX2D,CRAN2D)
        do nnw = 1,L(LN)%inwet
            i = L(LN)%inw(1,nnw) ; j = L(LN)%inw(2,nnw); nn = L(LN)%mn(i,j)
            !
            if (L(LN)%coupling_dim.eq.2.and.                                  &
                i.ge.L(LN)%is2.and.i.lt.L(LN)%ie2.and.                        &
                j.ge.L(LN)%js2.and.j.lt.L(LN)%je2) cycle   
            !Find the maximum velocity and courant number
            VMAX2D = max( VMAX2D, maxval( abs(L(LN)%U(:,nn)) ) )
            Cranx  = maxval( ( abs(L(LN)%U(:,nn)) + sqrt(g*L(LN)%D(:,nn)) )   &
                               * L(LN)%DT * L(LN)%RX(:,nn) )
            CRAN2D = max( CRAN2D, Cranx )
        enddo
!$omp end parallel do
        L(LN)%VMAX = VMAX2D ; L(LN)%CRAN = CRAN2D
    ENDDO 
endsubroutine check_courant
!
subroutine make_wall(L)
    use variables, only: ZERO, inflow,  ONEHALF, HALF, HCIRC, TWO
    use arrays, only: BC
    use type_list, only: Layer 
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    !----------------- Temporary variables ------------------------------------    
    real*8 :: Hmax
    integer :: nbw, nbe, i, j, nn
    logical :: whin = .false., shin = .false., ehin = .false., nhin = .false.
!=============== Makes a wall at the boundary where there is no input =========
! Finished Aug 22
! Updated Nov 11 2014, considering inflow angle
Hmax = -L%Z(L%ke)
if (inflow%wavetype(2:2).ne.'C') then     
    if (inflow%angle.gt.ONEHALF * HCIRC.or.                                   &
        inflow%angle.lt.HALF *    HCIRC    )                 whin = .true.
    if (inflow%angle.gt.HALF * HCIRC.and.                                     &
        inflow%angle.lt.ONEHALF * HCIRC  )                   ehin = .true.
    if (inflow%angle.gt.ZERO.and.inflow%angle.lt.HCIRC)      shin = .true.
    if (inflow%angle.gt.HCIRC.and.inflow%angle.lt.TWO*HCIRC) nhin = .true.
endif
if (BC(2)%EAST.eq.'WALL'.and..not.ehin) then
    i = L%ie
    do j = L%js,L%je-1
        nbe = L%mn(i,j)
        L%H(1,nbe) = Hmax
    enddo
endif
if (BC(2)%WEST.eq.'WALL'.and..not.whin) then
    i = L%is 
    do j = L%js,L%je-1
        nbw = L%mn(i,j)
        L%H(1,nbw) = Hmax  
    enddo
endif
if (L%DIM.eq.2) then
   if (BC(2)%NORTH.eq.'WALL'.and..not.nhin) then
        j = L%je
        do i = L%is,L%ie-1
            nbe = L%mn(i,j)
            L%H(2,nbe) = Hmax     
        enddo
    endif
    if (BC(2)%SOUTH.eq.'WALL'.and..not.shin) then
        j = L%js
        do i = L%is,L%ie-1
            nbw = L%mn(i,j)
            L%H(2,nbw) = Hmax                      
        enddo
    endif
endif    
!
endsubroutine make_wall 
!
subroutine alter_tide_and_remove_water_on_land(L)
    use variables, only: tide_level
    use type_list, only: Layer 
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    !----------------- Temporary variables ------------------------------------
    integer                            :: nn,i,j,icc,icc0,ii,iit
    integer,allocatable,dimension(:,:) :: nflag
    integer,allocatable,dimension(:)   :: ic,jc,ic0,jc0,istart,jstart
    real*8                             :: zb,th
    !
    allocate(nflag(L%is-1:L%ie,L%js-1:L%je))
    allocate(ic(L%ie+L%je),jc(L%ie+L%je),ic0(L%ie+L%je),jc0(L%ie+L%je))
    allocate(istart(2*L%ie+2*L%je),jstart(2*L%ie+2*L%je))
    ! Set nflag = 1 by default
    nflag = 1 ; istart = 0 ; jstart = 0 ; iit = 0
    ! Make F equal to zero everywhere
    do nn = 1,L%inne
        i = L%IN(1,nn) ; j = L%IN(2,nn)
        if (i.lt.L%is.or.i.ge.L%ie.or.j.lt.L%js.or.j.ge.L%je) cycle
        if (((i.eq.L%is.or.i.eq.L%ie-1).and.L%DIM.eq.1).or.                   &
           ((i.eq.L%is.or.i.eq.L%ie-1.or.j.eq.L%js.or.j.eq.L%je-1)            &
            .and.L%DIM.eq.2)) then
            if (L%ZKm(nn).le.tide_level) then
                ! Set initial cells to start from (on the boundary)
                iit = iit + 1
                istart(iit) = i ; jstart(iit) = j
            endif
        endif
        ! Everywhere where the land is low enough set nflag = 0
        nflag(i,j) = 0
        ! We already have F equal to zero at surface cells
        if (L%NF(i,L%KS,j).eq.0) cycle
        L%F(nn) = 0.0d0
        L%NF(i,L%KS,j) = 0 ; L%NFB(i,L%KS,j) = 0
    enddo
    ! Now add in fluid where required..  
    do ii = 1,iit
        ic0 = 0 ; jc0 = 0
        icc0 = 1
        ic0(1) = istart(ii) ; jc0(1) = jstart(ii)
        if ( nflag(ic0(1),jc0(1)) == 1 ) cycle
        do !Infinite loop until covergence is reached
            icc = 0 
            do nn = 1,icc0
                i = ic0(nn) ; j = jc0(nn)
                zb = L%ZKm(L%mn(i,j))
                ! Set new F value based on tide level
                L%F(L%mn(i,j)) = (tide_level - zb) / (L%Z(L%KE) - zb)
                L%NF(i,L%KS,j) = 2 ; L%NFB(i,L%KS,j) = -3
                ! Set nflag = 1
                nflag(i,j) = 1
                ! Search neighbouring cells for local land or wall height
                if (nflag(i+1,j).eq.0) then
                    th = max(-L%H(1,L%mn(i+1,j)),L%ZKm(L%mn(i+1,j)))
                    if (th.le.tide_level) then
                        icc     = icc + 1 ; nflag(i+1,j) = 1
                        ic(icc) = i + 1   ; jc(icc) = j
                    endif
                endif
                if (nflag(i-1,j).eq.0) then
                    th = max(-L%H(1,L%mn(i,j)),L%ZKm(L%mn(i-1,j)))
                    if (th.le.tide_level) then
                        icc     = icc + 1 ; nflag(i-1,j) = 1
                        ic(icc) = i - 1   ; jc(icc) = j
                    endif
                endif
                ! Cycle for 1D calculation
                if (L%DIM.eq.1) cycle
                if (nflag(i,j+1).eq.0) then
                    th = max(-L%H(2,L%mn(i,j+1)),L%ZKm(L%mn(i,j+1)))
                    if (th.le.tide_level) then
                        icc     = icc + 1 ; nflag(i,j+1) = 1
                        ic(icc) = i       ; jc(icc)      = j + 1
                    endif
                endif
                if (nflag(i,j-1).eq.0) then
                    th = max(-L%H(2,L%mn(i,j)),L%ZKm(L%mn(i,j-1)))
                    if (th.le.tide_level) then
                        icc     = icc + 1 ; nflag(i,j-1) = 1
                        ic(icc) = i       ; jc(icc)      = j - 1
                    endif
                endif
            enddo
            ! Check it convergence reached or not
            if (icc.ne.0) then
                icc0 = icc
                ic0(1:icc0) = ic(1:icc0)
                jc0(1:icc0) = jc(1:icc0)
            else
               ! Convergence reached
               continue 
               exit
            endif
        enddo
    enddo
    deallocate(nflag)
    deallocate(ic,jc,ic0,jc0)
    deallocate(istart,jstart)
endsubroutine alter_tide_and_remove_water_on_land
!    
subroutine equal_bound(inne,is,js,ie,je,mn,X,DIM,DIM2)
    implicit none
    !----------------- Input/Output variables ---------------------------------
    integer,intent(in) :: inne,is,js,ie,je,mn(ie+1,je+1),DIM,DIM2
    real*8,dimension(DIM,0:inne),intent(inout) :: X
    !----------------- Temporary variables used in the routine ----------------
    integer :: nn,i,j,nss,ne,ns,nee,iss,iee,ii
    !    SUBROUTINE: EQUAL_BOUND
    !    Sets the outer bounds of the input arrays
    !    X and Y equal to the adjacent cell
    !
    do nn = 1,DIM
        !Tangential
        if (nn.eq.1) then
            iss = is; iee = ie
        elseif (nn.eq.2) then
            iss = js; iee = je
        endif
!$omp parallel private(nss,ne,ns,nee)
!$omp do
        do ii = iss,iee
            if (nn.eq.1) then
                nss = mn(ii,js-1) ; ne = mn(ii,je-1)                        
                ns = mn(ii,js)    ; nee = mn(ii,je)
            elseif (nn.eq.2) then
                nss = mn(is-1,ii) ; ne = mn(ie-1,ii)
                ns = mn(is,ii)    ; nee = mn(ie,ii)
            endif
            if (nss.ne.0.and.ns.ne.0) X(nn,nss) = X(nn,ns)
            if (nee.ne.0.and.ne.ne.0) X(nn,nee) = X(nn,ne)
        enddo
!$omp end do
!$omp single
        !Normal
        if (nn.eq.1) then
            iss = js-1; iee = je
        elseif (nn.eq.2) then
            iss = is-1; iee = ie
        endif   
!$omp end single
!$omp do
        do ii = iss,iee
            if (nn.eq.1) then
                nss = mn(is-1,ii);ns = mn(is,ii)
                if (DIM2.eq.1) then
                    ne = mn(ie-1,ii);nee = mn(ie,ii)  
                else
                    nee = 0
                endif
            elseif (nn.eq.2) then
                nss = mn(ii,js-1);ns = mn(ii,js);nee = 0  
            endif
            if (nss.ne.0.and.ns.ne.0) X(nn,nss) = X(nn,ns)      
            if (nee.ne.0.and.ne.ne.0) X(nn,nee) = X(nn,ne)
        enddo
!$omp end do    
!$omp end parallel
    enddo 
!
endsubroutine equal_bound    
!
!!    
!#ifdef DDDDD    
!subroutine make_sponge(L)
!    use leapfrog_mod, only: Layer
!    use local_module, only: BC
!    type(Layer),intent(inout) :: L
!    real*8 :: Hmax
!    integer :: nbq,nbe,i,j,nn
!!=============== Finds the coefficients for the sponge layer ===================
!! Finished Oct 22
!if (BC(2)%NORTH(1:1).ne.'S'.and.BC(2)%EAST(1:1).ne.'S'                         &
!   .and.BC(2)%SOUTH(1:1).ne.'S'.and.BC(2)%WEST(1:1).ne.'S') return
!! Set to 1 as default value
!allocate(L%SPNGCF(1:L%inne))
!L%SPNGCF = 1.0d0
!do i=L%is,L%ie-1
!    do j=L%js,L%je-1
!        continue
!    enddo
!enddo
!!if (BC(1:1).eq.'S') then
!!    !Implement simple sponge layer of width equal to 
!!    !number of cells as defined in bound
!!    read(BC(3:4),'(I2)') w
!!    do jj = 0,w-1
!!        nbe = L(1)%mn(i+ix*(jj+1),j+iy*(jj+1))
!!        nbq = L(1)%mn(i+ix*jj+max(ix,0),j+iy*jj+max(iy,0))   
!!        L(1)%ETAn(nbe) = L(1)%ETAn(nbe) * sin( HPI_8 * real(jj+1/w) )**2 
!!        L(1)%Qn(:,nbq) = L(1)%Qn(:,nbq) * sin( HPI_8 * real(jj+1/w) )**2 
!!    enddo
!!endif
!!
!endsubroutine make_sponge 
!#endif    
    
