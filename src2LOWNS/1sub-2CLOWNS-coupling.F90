!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&          FIND_LAYER_INTERP : Gets the value of IS2, JS2, IE2, JE2       &&!
!&&        corresponding to the bounds of the next inner layer and          &&!
!&&  then calculates the mapping and interpolation matrices between layers  &&!
!&&                     BY PRINGLE (AUG 2014) AND UPDATED JUNE & DEC 2015   &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
! Update : June 25th; addition to algorithm for CFb so to discard 
!                     boundary cell being used in input condition
!          Dec  31st; altered to suit new format
subroutine find_layer_interp(L,L_out,LN)
    use variables, only: LNUM, ZERO, ONE, HALF, is, ie, js, je, IM
    use arrays, only: x, dx, mn2d
    use type_list, only: Layer
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    type(Layer),intent(in)    :: L_out
    integer,intent(in)        :: LN
!----------------- Temporary variables ----------------------------------------
    integer :: LOWER_DIM
    integer :: i, j, nn, ii, i2, j2, jj, jj1, jj2, iss, iee, jss, jee, js2, je2
    real*8  :: Xs, Xe, Ys, Ye, Xc, Yc, DXin, DXout, DYin, DYout
    integer,allocatable :: intemp(:,:)
    !
    if (IM.eq.2.and.LN.eq.LNUM) return ! No more layers inside
    !=========================================================================!
    ! GET THE VALUES OF IS2, JS2, IE2, JE2, MN2, IN2, INTER_MAP, INTER_MAPb,  !
    !                INTER_CFb TO MAP BETWEEN THE NEXT INNER LAYER            !
    !=========================================================================!
    !Set default values for is2, js2, ie2, je2 
    !(so that they will never affect the looping when no nesting is required)
    L%ie2 = L%is ; L%je2 = L%js
    L%is2 = L%ie ; L%js2 = L%je
    !
    if (IM.eq.1.and.LN.eq.LNUM) then
        ! Set bounds for the 3D domain
        Xs = X(0,is); Xe = X(0,ie)  
        Ys = X(1,js); Ye = X(1,je)  
    else
        ! Set bounds for the next layer in 2D domain
        Xs = L_out%X(L_out%is) ; Xe = L_out%X(L_out%ie)   
        Ys = L_out%Y(L_out%js) ; Ye = L_out%Y(L_out%je)     
    endif
    ! Need to calculate one more cell in for dispersion scheme 
    !if (L%DISP.eq.2) then
    !    iss = 1 ; jss = 1 ; jee = 1 ; iee = 1
    !else
        iss = 0 ; jss = 0 ; jee = 0 ; iee = 0
    !endif
    ! Find bounds for the next layer
    do i = L%is, L%ie-1
        if (L%X(i+1).gt.Xs) then
            L%is2 = i; exit
        endif
    enddo
    if (L%is2 == L%is) iss = -1
    do i = L%is2, L%ie
        if (L%X(i).ge.Xe) then
            L%ie2 = i; exit
        elseif (i.eq.L%ie) then
            write(6,*) 'IE: Warning 2D area is shorter than 3D...'
            L%ie2 = i
        endif
    enddo
    if (L%ie2 == L%ie) iee = -1
    do j = L%js, L%je-1
        if (L%Y(j+1).gt.Ys) then
            L%js2 = j; exit
        endif
    enddo
    if (L%js2 == L%js) jss = -1
    do j = L%js2, L%je    
        if (L%Y(j).ge.Ye) then
            L%je2 = j; exit
        elseif (j.eq.L%je) then
            write(6,*) 'JE: Warning 2D area is shorter than 3D...'
            L%je2 = j
        endif
    enddo
    if (L%je2 == L%je) jee = -1
    !=========================================================================!
    !         GET MN2, IN2, AND INTER_MAP TO MAP TO THE NEXT INNER LAYER      !
    !=========================================================================!
    TWO_WAY: if (L%coupling_dim.eq.2) then
    ! Loop over bounds and get number of cells within domain and mapping tool
    L%inne2 = 0 ; allocate( intemp(2,(L%ie2-L%is2)*(L%je2-L%js2)),            &
                            L%MN2(0:L%DIM,L%is2:L%ie2-1,L%js2:L%je2-1) )
    ! Set default value of mn2
    L%MN2 = 0
    do i = L%is2, L%ie2-1
        do j = L%js2, L%je2-1
            ! Cycle at non-cells
            if (L%mn(i,j).eq.0) cycle
            ! Input the values
            L%inne2 = L%inne2 + 1
            intemp(1,L%inne2) = i !(Map to normal i
            intemp(2,L%inne2) = j !(j integers)
            L%MN2(:,i,j) = L%mn(i,j)
            if (i.ge.L%is2+iss.and.i.lt.L%ie2-iee.and.(L%DIM.eq.1.or.         &
                (j.ge.L%js2+jss.and.j.lt.L%je2-jee))) L%MN2(0,i,j) = 0             
            ! April 03 changed ge le to gt lt
            if ((i.gt.L%is2+iss.and.i.lt.L%ie2-iee.and.(L%DIM.eq.1.or.        &  
               (j.ne.L%js2.and.j.ne.L%je2-1))).or.L%mn(i-1,j).eq.0) then 
                ! Only calculate flux on boundary
                L%MN2(1,i,j) = 0
            elseif (((iss.eq.-1.and.i.eq.L%is2+1).or.                         &
                    (iee.eq.-1.and.i.eq.L%ie2-1)).and.                        &
                    (L%DIM.eq.1.or.(j.ne.L%js2.and.j.ne.L%je2-1))) then 
                ! No information is passed through East or West boundary     
                L%MN2(1,i,j) = -1
            endif
            if (L%DIM.eq.1) cycle        
            ! Y-flux
            if ((j.gt.L%js2+jss.and.j.lt.L%je2-jee.and.                       &
                 i.ne.L%is2.and.i.ne.L%ie2-1).or.L%mn(i,j-1).eq.0) then
                L%MN2(2,i,j) = 0 
            elseif (((jss.eq.-1.and.j.eq.L%js2+1).or.                         &
                    (jee.eq.-1.and.j.eq.L%je2-1)).and.                        &
                    (i.ne.L%is2.and.i.ne.L%ie2-1)) then 
                ! No information is passed through North or South boundary     
                L%MN2(2,i,j) = -1
            endif
        enddo
    enddo
    !Allocate the size of the I,J mapping instrument 
    allocate(L%IN2(2,L%inne2))
    L%IN2(:,1:L%inne2) = intemp(:,1:L%inne2)
    deallocate(intemp)
    !=========================================================================!
    ! FIND THE MAPPING AND MATRIX FOR INNER LAYER TO OUTER LAYER CONVERSION   !
    !=========================================================================!
    if (LN.eq.LNUM) then
        ! 3D layer
        iss = is; iee = ie
        jss = js; jee = je   
        DXin = minval(DX(0,:), DX(0,:) > ZERO)
        DYin = minval(DX(1,:), DX(1,:) > ZERO)
        if (L%DIM.eq.1) DYin = ONE
    else !2DH layer
        iss  = L_out%is ; iee = L_out%ie
        jss  = L_out%js ; jee = L_out%je
        DXin = L_out%X(L_out%is+1) - L_out%X(L_out%is)
        DYin = L_out%Y(L_out%js+1) - L_out%Y(L_out%js)
        if (L%DIM.eq.1) DYin = ONE
    endif
    DXout = maxval(L%DX(1,:)) 
    if (L%DIM.eq.1) DYout = ONE
    if (L%DIM.eq.2) DYout = maxval(L%DX(2,:))
    !Allocate the size of the interpolation matrix
    ! For nesting between 2DH layers and 2DH and 3D models we need to
    ! map and interpolate the momentum flux
    allocate(L%INTER_MAP(0:L%DIM,nint(DXout/DXin) * nint(DYout/DYin),L%inne2))
             L%INTER_MAP = 0
    !Now find the interpolation matrix for each cell
    do nn = 1,L%inne2 
        ! Get i and j integers in outer layer
        i2 = L%IN2(1,nn); j2 = L%IN2(2,nn)
        jj1 = 0
        ! Get the bounding box
        Xs = L%X(i2)
        Xe = L%X(i2+1) 
        Ys = L%Y(j2)
        Ye = L%Y(j2+1)   
        ! FREE-SURFACE
        do i = iss,iee
            if (LN.eq.LNUM) then
                Xc = HALF * ( X(0,i) +  X(0,i+1) )
            else
                Xc = HALF * ( L_out%X(i) + L_out%X(i+1) )
            endif
            if (Xc.le.Xs) cycle
            if (Xc.ge.Xe) exit
            do j = jss,jee
                if (LN.eq.LNUM) then
                    Yc = HALF * ( X(1,j) +  X(1,j+1) )
                else
                    Yc = HALF * ( L_out%Y(j) + L_out%Y(j+1) )
                endif
                if (Yc.le.Ys) cycle
                if (Yc.ge.Ye) exit
                !Interpolation mapping and coefficients for ETA
                jj1 = jj1 + 1
                if (LN.eq.LNUM) then
                    !3D layer
                    L%INTER_MAP(0,jj1,nn) = mn2d(i,j)
                else!2DH layer
                    L%INTER_MAP(0,jj1,nn) = L_out%mn(i,j)
                endif
            enddo
        enddo
        ! FLUXES
        do ii = 1,L%DIM
            jj1 = 0
            if (ii.eq.1) then
                ! Get the bounding box
                Xs = HALF * (L%X(i2) + L%X(i2-1))
                Xe = HALF * (L%X(i2) + L%X(i2+1)) 
                Ys = L%Y(j2)
                Ye = L%Y(j2+1)          
            elseif (ii.eq.2) then
                ! Get the bounding box
                Xs = L%X(i2)
                Xe = L%X(i2+1)   
                Ys = HALF * (L%Y(j2) + L%Y(j2-1))
                Ye = HALF * (L%Y(j2) + L%Y(j2+1))
            endif
            do i = iss,iee
                if (LN.eq.LNUM) then
                    !3D layer
                    if (ii.eq.1) Xc = X(0,i)
                    if (ii.eq.2) Xc = HALF * ( X(0,i) +  X(0,i+1) )
                else!2DH layer
                    if (ii.eq.1) Xc = L_out%X(i)
                    if (ii.eq.2) Xc = HALF * ( L_out%X(i) + L_out%X(i+1) )
                endif
                if (Xc.le.Xs) cycle
                if (Xc.ge.Xe) exit
                do j = jss,jee
                    if (LN.eq.LNUM) then
                        !3D layer
                        if (ii.eq.1) Yc = HALF * ( X(1,j) +  X(1,j+1) )
                        if (ii.eq.2) Yc = X(1,j)
                    else!2DH layer
                        if (ii.eq.1) Yc = HALF * ( L_out%Y(j) + L_out%Y(j+1) )
                        if (ii.eq.2) Yc = L_out%Y(j)
                    endif
                    if (Yc.le.Ys) cycle
                    if (Yc.ge.Ye) exit
                    !Interpolation mapping and coefficients for Q(ii) in ii-dir
                    jj1 = jj1 + 1
                    if (LN.eq.LNUM) then
                        !3D layer
                        L%INTER_MAP(ii,jj1,nn) = mn2d(i,j)
                    else!2DH layer
                        L%INTER_MAP(ii,jj1,nn) = L_out%mn(i,j)
                    endif
                enddo
            enddo  
        enddo
    enddo
    endif TWO_WAY
    !=========================================================================!
    !  FIND THE MAPPING AND INTERP MATRICES FOR OUTER LAYER CONVERSION ON THE !
    !                       BOUNDARY OF THE INNER LAYER                       !
    !=========================================================================!
    if (LN.eq.LNUM) then
        !3D domain
        if (L%DIM.eq.1) then
            iss = js ; iee = je
        else
            iss = min(js,is); iee = max(ie,je)
        endif
    else  !2D domain
        if (L%DIM.eq.1) then
            iss = L_out%js; iee = L_out%je
        else
            iss = min(L_out%js,L_out%is); iee = max(L_out%ie,L_out%je)
        endif
    endif
    !To interpolate Q in outer layer to Q in inner layer on the boundary
    !INTER_MAPb/CFb: (L%DIM,TWO INTERPOLATING INTEGERS,BOUNDARY LENGTH)
    allocate(L%INTER_MAPb(L%DIM,L%DIM*2,iss:iee),                             & 
             L%INTER_CFb (L%DIM,L%DIM*2,iss:iee)) 
    L%INTER_MAPb = 0 ; L%INTER_CFb = ZERO
    !Loop from East/West to North/South
    do jj = 1 , L%DIM
        if (jj.eq.1) then
            !East/West
            jj1 = 1     ; jj2 = 2
            js2 = L%js2 ; je2 = L%je2
            if (LN.eq.LNUM) then
                !3D domain
                iss = js       ; iee = je
            else!2D domain
                iss = L_out%js ; iee = L_out%je
            endif
        elseif (jj.eq.2) then
            !North/South
            jj1 = 3     ; jj2 = 4
            js2 = L%is2 ; je2 = L%ie2
            if (LN.eq.LNUM) then
                !3D domain
                iss = is       ; iee = ie
            else!2D domain
                iss = L_out%is ; iee = L_out%ie
            endif
        endif
        !Loop over the dimension
        do ii = 1 , L%DIM
            !Loop over the length of the inner layer
            do j = iss , iee
                !Get the center of the flux coordinate in inner layer
                if (LN.eq.LNUM) then
                    !3D domain
                    if (jj.eq.1.and.ii.eq.1) then
                        Xc = HALF * ( X(1,j) + X(1,j+1) )
                    elseif (jj.eq.1.and.ii.eq.2) then
                        Xc = X(1,j)
                    elseif (jj.eq.2.and.ii.eq.1) then
                        Xc = X(0,j)
                    elseif (jj.eq.2.and.ii.eq.2) then
                        Xc = HALF * ( X(0,j) + X(0,j+1) )
                    endif
                else!2D domain
                    if (jj.eq.1.and.ii.eq.1) then
                        Xc = HALF * ( L_out%Y(j) + L_out%Y(j+1) )
                    elseif (jj.eq.1.and.ii.eq.2) then
                        Xc = L_out%Y(j)
                    elseif (jj.eq.2.and.ii.eq.1) then
                        Xc = L_out%X(j)
                    elseif (jj.eq.2.and.ii.eq.2) then
                        Xc = HALF * ( L_out%X(j) + L_out%X(j+1) )
                    endif
                endif
                !Loop over the outer layer
                do j2 = js2-1 , je2+1
                    !Get the surrounding flux coordinates in outer layer
                    if (jj.eq.1.and.ii.eq.1) then
                        Xs = HALF * ( L%Y(j2)   + L%Y(j2+1) )
                        Xe = HALF * ( L%Y(j2+1) + L%Y(j2+2) )   
                    elseif (jj.eq.1.and.ii.eq.2) then
                        Xs = L%Y(j2) ; Xe = L%Y(j2+1)
                    elseif (jj.eq.2.and.ii.eq.1) then
                        Xs = L%X(j2) ; Xe = L%X(j2+1)
                    elseif (jj.eq.2.and.ii.eq.2) then
                        Xs = HALF * ( L%X(j2)   + L%X(j2+1) )
                        Xe = HALF * ( L%X(j2+1) + L%X(j2+2) )                              
                    endif
                    !Allocate the mapping and coefficient
                    !values if we find the right location
                    if (Xe.ge.Xc.and.Xs.lt.Xc) then
                        if (jj.eq.1.and.ii.eq.1) then
                            if (L%NF(L%is,L%ks,j2).eq.-1.and.                 &
                                L%NFB(L%is,L%ks,j2).ne.0) then
                                ! Cell below is along boundary
                                Xe = Xc ; Xs = ZERO
                            elseif (L%NF(L%is,L%ks,j2+1).eq.-1.and.           &
                                    L%NFB(L%is,L%ks,j2+1).ne.0) then
                                ! Cell above is along boundary
                                Xs = Xc ; Xe = ZERO
                            endif
                        endif
                        if (jj.eq.2.and.ii.eq.2) then
                            if (L%NF(j2,L%ks,L%js).eq.-1.and.                 &
                                L%NFB(j2,L%ks,L%js).ne.0) then
                                ! Cell to the left is along boundary
                                Xe = Xc ; Xs = ZERO
                            elseif (L%NF(j2+1,L%ks,L%js).eq.-1.and.           &
                                    L%NFB(j2+1,L%ks,L%js).ne.0) then
                                ! Cell to the right is along boundary
                                Xs = Xc ; Xe = ZERO
                            endif
                        endif
                        L%INTER_MAPb(ii,jj1,j) = j2
                        L%INTER_CFb(ii,jj1,j)  = ( Xe - Xc ) / ( Xe - Xs )  
                        L%INTER_MAPb(ii,jj2,j) = j2 + 1
                        L%INTER_CFb(ii,jj2,j)  = ( Xc - Xs ) / ( Xe - Xs ) 
                        exit
                    endif
                enddo
            enddo
        enddo
    enddo
endsubroutine find_layer_interp    
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&        INT_ADJ_DT : MAKES SURE THAT DT IN ALL LAYERS ARE                &&! 
!&&                     INTEGER MULTIPLES OF EACH OTHER                     &&!
!&&                                                       - PRINGLE (2016)  &&!     
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!  
subroutine int_adj_dt(iflag)
    use variables, only: DT, DTM, LNUM, ZERO, IM
    use arrays, only: L 
    implicit none
    integer,intent(in) :: iflag
    integer :: LN, DTRAT, ix = 0
    real*8  :: DT_in
!   
    if (iflag.eq.0) then
        !Initialisation of DT
        if (IM.eq.1) ix = 1 
        do LN = 1,LNUM + ix
            if (LN.eq.LNUM+ix) exit
            if (LN.eq.LNUM) then
                DT_in = DT
                DTM   = DT
            else
                DT_in = L(LN+1)%DT
            endif
            !We can cycle if DT is already an integer multiple
            if (mod(L(LN)%DT,DT_in).eq.ZERO) cycle
            !Get integer ratio of outer layer over inner layer
            DTRAT  = ceiling( L(LN)%DT / DT_in )
            !Change DT in inner layer based on integer ratio
            if (LN.eq.LNUM) then
                DT  = L(LN)%DT / real(DTRAT)
                DTM = DT
            else
                L(LN+1)%DT = L(LN)%DT / real(DTRAT)
            endif
        enddo
        ! Make sure Courant numbers are OK
        call check_courant
    elseif (iflag.eq.1) then
        !During hybrid calculation when 3D timestep changes
        !We can return if DT is already an integer multiple
        if (mod(L(LNUM)%DT,DT).eq.ZERO) return
        !Get integer ratio of outer layer over inner layer
        DTRAT  = ceiling( L(LNUM)%DT / DT )
        !Change DT in 3D layer based on integer ratio
        DT = L(LNUM)%DT / real(DTRAT)
        !Make sure DTM is not larger than L(LNUM)%DT
        DTM = min(DTM,L(LNUM)%DT)
    endif
!
endsubroutine int_adj_dt
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                       SUBROUTINE: interp_outer_flux                     &&!
!&&       INTERPOLATES THE FLUXES ON THE BOUNDARY OF THE NEXT LAYER         &&!
!&&                       (ONLY FOR 2D - 2D NESTING)                        &&!
!&&        WRITTEN BY PRINGLE (2014) & UPDATED BY PRINGLE IN JUNE, DEC 2015 &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
subroutine interp_outer_flux(L_out,L_in,dt_cf)
    use type_list, only: Layer 
    use variables, only: ZERO, ONE, HALF, V_lim
    implicit none
    !----------------- Input/Output variables ---------------------------------
    real*8,intent(in)         :: dt_cf
    type(Layer),intent(in)    :: L_out
    type(Layer),intent(inout) :: L_in
    !----------------- Temporary variables used in the routine ----------------
    integer :: bound, iss, iee, ib_in, ib_out_m, ib_out
    integer :: nn_in, nn_out_1, nn_out_2, nn_out_m1, nn_out_m2
    integer :: i, j, j1, j2, jj, ii
    real*8  :: c1 , c2, c_out, add, wdiv
! INTER_MAPb/CFb: (L%DIM,TWO INTERPOLATING INTEGERS,BOUNDARY LENGTH)
! ALL BOUNDARIES
do bound = 1,L_out%DIM ! Loop over the dimension (for # of boundaries)
    if (bound.eq.1) then
        ! West and East
        iss = L_in%js ; iee = L_in%je
    elseif (bound.eq.2) then
        ! North and South
        iss = L_in%is ; iee = L_in%ie
    endif
!$omp parallel do private (ii,j1,j2,c1,c2,c_out,jj,add,ib_in,ib_out_m,ib_out, &
!$omp                      nn_in,nn_out_m1,nn_out_m2,nn_out_1,nn_out_2) 
    do j = iss,iee-1
        do ii = 1,L_out%DIM !Loop over the dimension (for flux dimension)
            j1 = L_out%INTER_MAPb(ii,2*bound-1,j)
            c1 = L_out%INTER_CFb(ii,2*bound-1,j)
            j2 = L_out%INTER_MAPb(ii,2*bound,j) 
            c2 = L_out%INTER_CFb(ii,2*bound,j)
            do jj = 1,2 !Loop from west to east or north to south
                if (bound.eq.1) then     
                    if (jj.eq.1) then     !West 
                        ib_in = L_in%is-ii+1 ; ib_out = L_out%is2-ii+1
                        if (ii.eq.2) ib_out_m = ib_out + 1
                    elseif (jj.eq.2) then !East
                        ib_in = L_in%ie      ; ib_out = L_out%ie2
                        if (ii.eq.2) ib_out_m = ib_out - 1
                    endif
                    nn_in    = L_in%mn(ib_in,j)
                    if (ii.eq.2) then
                        nn_out_m1 = L_out%mn(ib_out_m,j1)
                        nn_out_m2 = L_out%mn(ib_out_m,j2)
                    endif
                    nn_out_1 = L_out%mn(ib_out,j1)
                    nn_out_2 = L_out%mn(ib_out,j2) 
                elseif (bound.eq.2) then   
                    if (jj.eq.1) then     !South
                        ib_in = L_in%js+ii-2 ; ib_out = L_out%js2+ii-2
                        if (ii.eq.1) ib_out_m = ib_out + 1
                    elseif (jj.eq.2) then !North
                        ib_in = L_in%je      ; ib_out = L_out%je2
                        if (ii.eq.1) ib_out_m = ib_out - 1
                    endif
                    nn_in    = L_in%mn(j,ib_in)
                    if (ii.eq.1) then
                        nn_out_m1 = L_out%mn(j1,ib_out_m)
                        nn_out_m2 = L_out%mn(j2,ib_out_m)
                    endif
                    nn_out_1 = L_out%mn(j1,ib_out)
                    nn_out_2 = L_out%mn(j2,ib_out) 
                endif
                ! Ignore if the depth in one of the cells is zero, or if there
                ! is a higher ground level than the water level between layers 
                if (L_in%Dn(ii,nn_in).eq.ZERO.or.                             &
                    L_out%Dn(ii,nn_out_1).eq.ZERO.or.                         &
                    L_out%Dn(ii,nn_out_2).eq.ZERO.or.                         &
                    L_in%ZK(nn_in) + L_in%Dmin.gt.&
                    L_out%ETAn(nn_out_1) * c1 + L_out%ETAn(nn_out_2) * c2) then
                    L_in%Qn(ii,nn_in) = ZERO;  cycle
                endif
                ! Getting momentum
                if (dt_cf.eq.ONE) then
                    ! Linear interpolation in space only
                    L_in%Qn(ii,nn_in) = L_out%Qn(ii,nn_out_1) * c1            &
                                      + L_out%Qn(ii,nn_out_2) * c2
                else
                    ! Bi-linear interpolation in space and time
                    L_in%Qn(ii,nn_in) = ( L_out%Qn(ii,nn_out_1) * c1          &
                                  + L_out%Qn(ii,nn_out_2) * c2 ) * dt_cf      &           
                                  + ( L_out%Q(ii,nn_out_1) * c1               &
                                  + L_out%Q(ii,nn_out_2) * c2 ) * (ONE - dt_cf)           
                endif
                ! For tangential momentum to get the correct location we need 
                ! to interpolate in the jj direction too
                if (bound.ne.ii) then
                    c_out = HALF * ( L_out%DX(ii,nn_out_1)                    &
                          + L_in%DX(ii,nn_in) ) / L_out%DX(ii,nn_out_1)
                    if (dt_cf.eq.ONE) then
                        ! Linear interpolation in space only
                        add = L_out%Qn(ii,nn_out_m1) * c1                     &
                            + L_out%Qn(ii,nn_out_m2) * c2
                    else
                        ! Bi-linear interpolation in space and time
                        add = ( L_out%Qn(ii,nn_out_m1) * c1                   &
                            +   L_out%Qn(ii,nn_out_m2) * c2 ) * dt_cf         &           
                            + ( L_out%Q(ii,nn_out_m1)  * c1                   &
                            +   L_out%Q(ii,nn_out_m2)  * c2 ) * (ONE - dt_cf)
                    endif
                    L_in%Qn(ii,nn_in) = L_in%Qn(ii,nn_in) * c_out             &
                                      + add * (ONE - c_out) 
                endif
                ! Limiting input flux (to avoid ridiculous velocities)
                if (wdiv(L_in%Qn(ii,nn_in),L_in%D(ii,nn_in)).gt.V_lim) then
                    L_in%Qn(ii,nn_in) = V_lim * L_in%D(ii,nn_in)
                endif
            enddo
        enddo
    enddo
!$omp end parallel do
enddo                 
endsubroutine interp_outer_flux
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                     SUBROUTINE: interp_inner_flux                       &&!
!&&  INTERPOLATES THE FLUXES ON THE INSIDE INNER LAYER TO THE OUTER LAYER   &&!
!&&              (FOR BOTH 2D - 2D NESTING & 2D - 3D COUPLING)              &&!
!&&        WRITTEN BY PRINGLE (2014) & UPDATED BY PRINGLE IN JUNE, DEC 2015 &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
subroutine interp_inner_flux(L_out,L_in,LN)
    use type_list, only: Layer
    use variables, only: LNUM, ZERO, HALF, VSMALL
    use arrays, only: IX, mn2d, in2d, Flux_3D, IWD, FSL
    implicit none
    !----------------- Input/Output variables ---------------------------------
    integer,intent(in)        :: LN
    type(Layer),intent(in)    :: L_in !Is equal to L_out if inner layer is 3D
    type(Layer),intent(inout) :: L_out
    !----------------- Temporary variables used in the routine ----------------
    integer :: nn, nn_out, i , j , ii, iii, icount
    real*8  :: Qsum, Dsum, Wsum, wdiv
    !
!$omp parallel do private (ii,i,j,nn_out,iii,icount,Qsum,Dsum,Wsum) 
do nn = 1, L_out%inne2 !Loop over the outer layer cells in the inner layer
    ! Mapping to the normal integers in the outer domain
    i = L_out%IN2(1,nn); j = L_out%IN2(2,nn) ; nn_out = L_out%mn(i,j)
    ! Pass through the free-surface, ETA 
    Wsum = ZERO ; icount = 0
    if (LN.eq.LNUM) then
        ! 3D domain 
        do iii = 1,ubound(L_out%INTER_MAP,2)
            if (L_out%INTER_MAP(0,iii,nn).ne.0) then
                if (FSL(L_out%INTER_MAP(0,iii,nn))                            &
                    + IWD(L_out%INTER_MAP(0,iii,nn)).lt.L_out%Dmin) cycle
                icount = icount + 1
                Wsum = Wsum + FSL(L_out%INTER_MAP(0,iii,nn))
            endif
        enddo
    else! 2DH domain
        do iii = 1,ubound(L_out%INTER_MAP,2)
            if (L_out%INTER_MAP(0,iii,nn).ne.0) then
                if (L_in%ETAn(L_out%INTER_MAP(0,iii,nn)) -                    &
                    L_in%ZK(L_out%INTER_MAP(0,iii,nn)).lt.L_out%Dmin) cycle
                icount = icount + 1
                Wsum = Wsum + L_in%ETAn(L_out%INTER_MAP(0,iii,nn))
            endif
        enddo
    endif
    if (icount.ne.0) then
        L_out%ETA(nn_out) = max( Wsum / real(icount) , L_out%ZK(nn_out) )
    else
        L_out%ETA(nn_out) = L_out%ZK(nn_out)
    endif
    if (abs(L_out%ETA(nn_out)).lt.VSMALL) L_out%ETA(nn_out) = ZERO
    L_out%ETAn(nn_out) = L_out%ETA(nn_out)
    do ii = 1, L_out%DIM !Loop over the number of dimensions
        ! Cycle on the input boundary 
        if (ii.eq.1.and.i.eq.L_out%is2) cycle
        if (ii.eq.2.and.j.eq.L_out%js2) cycle
        ! Cycle at small depths in outer layer
        if (L_out%ETAn(nn_out) - L_out%ZK(nn_out).lt.L_out%Dmin.or.           &
            L_out%ETAn(L_out%mn(i-ix(ii,1),j-ix(ii,2))) -                     &
            L_out%ZK(L_out%mn(i-ix(ii,1),j-ix(ii,2))).lt.L_out%Dmin) then
            L_out%Qn(ii,nn_out) = ZERO 
            L_out%D(ii,nn_out)  = ZERO
            cycle
        endif
        ! Get the mean of the fluxes and depths in the inner layer that 
        ! fit within the bounding box of the outer layer
        Qsum = ZERO ; Dsum = ZERO ; icount = 0 
        do iii = 1,ubound(L_out%INTER_MAP,2)
            if (L_out%INTER_MAP(ii,iii,nn).ne.0) then
                icount = icount + 1
                if (LN.eq.LNUM) then
                    ! 3D domain 
                    Qsum = Qsum + Flux_3D(ii-1,L_out%INTER_MAP(ii,iii,nn))
                    Dsum = Dsum + HALF * (FSL(L_out%INTER_MAP(ii,iii,nn))     &
                                        + IWD(L_out%INTER_MAP(ii,iii,nn))     &
                         + FSL(mn2d(in2d(0,L_out%INTER_MAP(ii,iii,nn))+ii-2,  &
                                    in2d(1,L_out%INTER_MAP(ii,iii,nn))-ii+1)) &
                         + IWD(mn2d(in2d(0,L_out%INTER_MAP(ii,iii,nn))+ii-2,  &
                                    in2d(1,L_out%INTER_MAP(ii,iii,nn))-ii+1))) 
                else! 2DH domain
                    Qsum = Qsum + L_in%Qn(ii,L_out%INTER_MAP(ii,iii,nn))
                    Dsum = Dsum + L_in%D(ii,L_out%INTER_MAP(ii,iii,nn))
                endif
            endif
        enddo
        L_out%Qn(ii,nn_out) = wdiv( Qsum , real(icount) )
        L_out%D (ii,nn_out) = wdiv( Dsum , real(icount) )
        if (L_out%D(ii,nn_out).lt.L_out%Dmin)    L_out%D (ii,nn_out) = ZERO
        if (abs(L_out%Qn(ii,nn_out)).lt.VSMALL)  L_out%Qn(ii,nn_out) = ZERO
    enddo
enddo
!$omp end parallel do
endsubroutine interp_inner_flux
!
!Only subroutines for coupling with 3D model from here on...
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                     SUBROUTINE: calc3D_start                            &&!
!&&     GETS THE MOMENTUM ON THE BOUNDARY OF THE LAST 2DH DOMAIN TO CHECK   &&!
!&&                WHETHER THE 3D CALCULATION SHOULD START OR NOT           &&!
!&&                                             WRITTEN BY PRINGLE (2016)   &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
subroutine calc3D_start(L)
    use interface_list, only: pritime2
    use type_list, only: Layer
    use arrays, only: date_time_start, date_time_old
    use variables, only: n, ZERO, Limit_3D, nstart
    implicit none
    !----------------- Input/Output variables ---------------------------------
    type(Layer),intent(inout) :: L
    !----------------- Temporary variables used in the routine ----------------
    real*8  :: Umax
    integer :: nn, ii, i, j
    !
    Umax = ZERO
!$omp parallel do private (ii,j) reduction(max:Umax) 
    ! Loop over all the cells in the 3D domain
    do i = L%is2, L%ie2
        do j = L%js2, L%je2
            do ii = 1, L%DIM !Loop over the number of dimensions
                ! Cycle inside the domain (only check the boundaries)
                if (ii.eq.1) then
                    if (j.eq.L%je2) cycle
                    if (i.gt.L%is2.and.i.lt.L%ie2) cycle
                elseif (ii.eq.2) then
                    if (i.eq.L%ie2) cycle
                    if (j.gt.L%js2.and.j.lt.L%je2) cycle
                endif
                Umax = max( Umax, abs(L%U(ii,L%mn(i,j))) )
            enddo
        enddo
    enddo
!$omp end parallel do
    !
    if (Umax.gt.Limit_3D) then
        nstart = n
        write(6,*) 'Start 3D Calculation, nstart=',nstart
        call pritime2(date_time_start,date_time_old)
        call flush(6) 
    else
        return
    endif
endsubroutine calc3D_start
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&              SUBROUTINE: Interp_2D_bound_to_3D                           &&!
!&&    CONVERTS THE MOMENTUM AND FREE-SURFACE ON THE BOUNDARY OF THE LAST    &&!
!&&    2DH DOMAIN INTO HORIZONTAL AND VERTICAL VELOCITY DISTRIBUTION IN 3D   &&!
!&&                         WRITTEN BY PRINGLE (2014) & UPDATED JUNE 2015    &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!    
subroutine Interp_2D_bound_to_3D(ii,ix,iy,i,j,DIM,DT_CF,bound)
    use variables, only: LNUM, VSMALL, ZERO, ONE, TWO, HALF, SIXTH, ks, ke,   &
                         h_in, cm, nut_b, nu, TI_b, akini, eini
    use arrays, only: L, mn, Un, U, F, FP, x, nff, dx, a, mn2d, U_3D, IWD,    &
                      r, r1, d 
    use interface_list, only: wdiv, LI_1
    implicit none
    !----------------- Input/Output variables ---------------------------------
    integer,intent(in)     :: ii, ix, iy, i, j, DIM
    real*8,intent(in)      :: DT_CF
    character*2,intent(in) :: bound
    !----------------- Temporary variables used in the routine ----------------
    logical :: exist
    integer :: iii, nn_1, nn_2, nn_3, jj, isign, k, i4, jjj, cc, count
    integer :: nnxm, nnxp, nnxp2, nn_22, ixv(-2:1), iyv(-2:1)
    real*8  :: znow, wll, zkk, U_av, h_av, V_av, V_t, open_rat
    real*8  :: RMSEn, RMSEt, Muln, Mult, DX_R, DX_N, DUDX, D2UDX2
    real*8  :: wl(-2:1), zk(-2:1), c(-1:0)
    real*8,dimension(-1:1) :: m, h, dp, U2D
!=============================================================================!
!   Calculation to assign 2D momentum as velocity distribution in 3D domain   !
!=============================================================================!
!ii -> cell integer, iii -> dimension, ix & iy -> direction
if (iy.eq.0) then
    ! Get the necessary integer numbers ixv & iyv (normal direction)
    if (ix.eq.1) then      ! West
        ixv(-2:1) = [L(LNUM)%is2-2, L(LNUM)%is2-1, L(LNUM)%is2, L(LNUM)%is2+1]
    elseif (ix.eq.-1) then ! East
        ixv(-2:1) = [L(LNUM)%ie2-2, L(LNUM)%ie2-1, L(LNUM)%ie2, L(LNUM)%ie2+1]
    endif
    ! Find the representative cell size of RANS and NSW
    DX_R = dx(0,i+ix) ; isign = ix 
    DX_N = L(LNUM)%DX(abs(iy)+1,L(LNUM)%mn(ixv(0),L(LNUM)%js2))
elseif (ix.eq.0) then
    ! Get the necessary integer numbers ixv & iyv (normal direction)
    if (iy.eq.1) then      ! South
        iyv(-2:1) = [L(LNUM)%js2-2, L(LNUM)%js2-1, L(LNUM)%js2, L(LNUM)%js2+1]
    elseif (iy.eq.-1) then ! North
        iyv(-2:1) = [L(LNUM)%je2-2, L(LNUM)%je2-1, L(LNUM)%je2, L(LNUM)%je2+1]
    endif
    ! Find the representative cell size of RANS and NSW
    DX_R = dx(1,j+iy) ; isign = iy
    DX_N = L(LNUM)%DX(abs(iy)+1,L(LNUM)%mn(L(LNUM)%is2,iyv(0))) 
endif
DUDX = ZERO
DIM_LOOP: do iii = 1, DIM
    if ((iii.eq.1.and.iy.eq.0).or.(iii.eq.2.and.ix.eq.0)) then
        !<<<<<<<<<<<<<<<<<< Normal direction <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        ! Get the necessary integer numbers ixv & iyv (tangential direction) 
        if (iy.eq.0) then
            iyv(-1) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+1,ii)
            iyv( 0) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+2,ii)
        elseif (ix.eq.0) then
            ixv(-1) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+1,ii)
            ixv( 0) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+2,ii)
        endif
        ! Get the interpolation coefficients (tangential direction)
        c(-1) = L(LNUM)%INTER_CFb(iii,abs(iy)*2+1,ii)
        c( 0) = L(LNUM)%INTER_CFb(iii,abs(iy)*2+2,ii) 
        wll = ZERO; zkk = ZERO; h_av = ZERO; U2D = ZERO
        do cc = -1,0
            do jj = -2,1
                if (iy.eq.0) nn_1  = L(LNUM)%mn(ixv(jj),iyv(cc))  
                if (ix.eq.0) nn_1  = L(LNUM)%mn(ixv(cc),iyv(jj))   
                ! Get the water level from extrapolation/interpolation
                ! in time and space 
                ! (We need WL^(n:n+1) but only know WL^(n+1/2) & WL^(n-1/2))
                wl(jj) = L(LNUM)%ETAn(nn_1) * (HALF + DT_CF)                  &
                       + L(LNUM)%ETA(nn_1) * (HALF - DT_CF)
                zk(jj) = L(LNUM)%ZK(nn_1) 
                if (jj.eq.-2) cycle
                ! Find the momentum fluxes from interpolation in time and space
                ! (We need M^(n:n+1) and know M^(n+1) & M^n)
                m(jj) = L(LNUM)%Qn(iii,nn_1) * DT_CF                          &
                      + L(LNUM)%Q(iii,nn_1) * (ONE - DT_CF)
                ! Find the initial water depth from interpolation in space
                h(jj) = L(LNUM)%H(iii,nn_1)
                ! Find the total water depth
                if (-h(jj).gt.zk(jj-1).and.-h(jj).gt.zk(jj)) then
                    ! In case of a wall
                    dp(jj) = max( max(wl(jj-1),wl(jj)) + h(jj), ZERO)
                else
                    ! The normal case 
                    dp(jj) = max( HALF * (wl(jj-1) + wl(jj)) + h(jj), ZERO)
                endif
                if (wl(jj-1)-zk(jj-1).lt.L(LNUM)%Dmin.or.                     &
                    wl(jj)-zk(jj).lt.L(LNUM)%Dmin.or.                         &
                    dp(jj).lt.L(LNUM)%Dmin) dp(jj) = ZERO
                ! Finally find the depth averaged velocity
                U2D(jj) = U2D(jj) + wdiv( m(jj) , dp(jj) ) * c(cc)
            enddo
            if (wl(min(-isign,0))-zk(min(-isign,0)).gt.L(LNUM)%Dmin.and.      &
                wl(min(isign,0))-zk(min(isign,0)).gt.L(LNUM)%Dmin) then
                if (cc.eq.-1.or.h_av.ne.ZERO) then
                    wll = wll + ( wl(min(-isign,0)) * (DX_N + DX_R)           &
                        + wl(min(isign,0)) * (DX_N - DX_R) ) * c(cc)
                    zkk = zkk + ( zk(min(-isign,0)) * (DX_N + DX_R)           &
                        + zk(min(isign,0)) * (DX_N - DX_R) ) * c(cc) 
                    h_av = h_av - HALF * ( zk(min(-isign,0))                  &
                                         + zk(min(isign,0))) * c(cc)
                else
                    wll = wll + wl(min(-isign,0)) * (DX_N + DX_R)             &
                        + wl(min(isign,0)) * (DX_N - DX_R)
                    zkk = zkk + zk(min(-isign,0)) * (DX_N + DX_R)             &
                        + zk(min(isign,0)) * (DX_N - DX_R)
                    h_av = h_av - HALF * (zk(min(-isign,0)) + zk(min(isign,0)))
                endif
            elseif (cc.eq.0.and.c(-1).ne.ZERO) then
                wll = wll / c(-1) ; zkk = zkk / c(-1) ; h_av = h_av / c(-1) 
            endif
        enddo
        !---------------------------------------------------------------------!
        ! Get the water level and bottom required as bc for RANS              !
        !---------------------------------------------------------------------!
        wll = HALF * wll / DX_N 
        zkk = HALF * zkk / DX_N
        ! Return if the difference in 2D depth and 3D depth is too large
        !if (abs(zkk + IWD(mn2d(i+ix,j+iy))) /                                 &
        !   abs(IWD(mn2d(i+ix,j+iy))).gt.0.2d0) return
        !---------------------------------------------------------------------!
        !Get the depth-averaged normal velocity (no normal interpolation reqd)!
        !---------------------------------------------------------------------!
        U_av = U2D(0)
        !---------------------------------------------------------------------!
        ! First-order gradient used in determining vertical velocity profile  !
        !---------------------------------------------------------------------!
        DUDX = DUDX + ( U2D(0) * DX_R + HALF * (U2D(isign) * (DX_N - DX_R)    &
             - U2D(-isign) * (DX_N + DX_R) ) ) / (DX_N * DX_N)
        !---------------------------------------------------------------------!
        ! Second-order difference used in determining horizontal profile      ! 
        !---------------------------------------------------------------------!
        D2UDX2 = (U2D(1) - U2D(0) + U2D(-1) - U2D(0)) / ( DX_N * DX_N )
    else!<<<<<<<<<<<<<<<<<< Tangential direction <<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        do i4 = 0,1 ! Loop for current vel and the next one to get gradient
            if (iy.eq.0) then
                iyv(-1) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+1,ii+i4)
                iyv( 0) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+2,ii+i4)
            elseif (ix.eq.0) then
                ixv(-1) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+1,ii+i4)
                ixv( 0) = L(LNUM)%INTER_MAPb(iii,abs(iy)*2+2,ii+i4)
            endif
            ! Get the interpolation coefficients (tangential direction)
            c(-1) = L(LNUM)%INTER_CFb(iii,abs(iy)*2+1,ii+i4)
            c( 0) = L(LNUM)%INTER_CFb(iii,abs(iy)*2+2,ii+i4) 
            if (iy.eq.0) iyv(-2) = iyv(-1) - 1
            if (ix.eq.0) ixv(-2) = ixv(-1) - 1
            U2D = ZERO
            do jjj = -1,0
                do cc = -1,0
                    do jj = -2,-1
                        if (iy.eq.0) then
                            nn_1  = L(LNUM)%mn(ixv(jjj),iyv(jj))
                            if (nn_1.eq.0) then
                                nn_1  =  L(LNUM)%mn(ixv(jjj),iyv(jj)+cc+1)  
                            endif
                        elseif (ix.eq.0) then
                            nn_1  = L(LNUM)%mn(ixv(jj),iyv(jjj))  
                            if (nn_1.eq.0) then
                                nn_1  =  L(LNUM)%mn(ixv(jj+cc+1),iyv(jjj))     
                            endif
                        endif
                        ! Get the water level from extrapolation/interpolation
                        ! in time and space (We need WL^(n:n+1) but 
                        !                    only know WL^(n+1/2) & WL^(n-1/2))
                        wl(jj) = L(LNUM)%ETAn(nn_1) * (HALF + DT_CF)          &
                               + L(LNUM)%ETA(nn_1) * (HALF - DT_CF)  
                        zk(jj) = L(LNUM)%ZK(nn_1)
                        if (jj.eq.-2) cycle
                        ! Find the momentum fluxes from interpolation in time
                        ! and space (We need M^(n:n+1) and know M^(n+1) & M^n)
                        m(jj) = L(LNUM)%Qn(iii,nn_1) * DT_CF                  &
                              + L(LNUM)%Q(iii,nn_1) * (ONE - DT_CF)
                        ! Find the initial water depth 
                        ! from interpolation in space
                        h(jj)  = L(LNUM)%H(iii,nn_1)
                        ! Find the total water depth
                        if (-h(jj).gt.zk(jj-1).and.-h(jj).gt.zk(jj)) then
                            ! In case of a wall
                            dp(jj) = max( max(wl(jj-1),wl(jj)) + h(jj), ZERO)
                        else
                            ! The normal case 
                            dp(jj) = max( HALF * (wl(jj-1) + wl(jj))          &
                                          + h(jj), ZERO )
                        endif
                        if (wl(jj-1)-zk(jj-1).lt.L(LNUM)%Dmin.or.             &
                            wl(jj)-zk(jj).lt.L(LNUM)%Dmin.or.                 &
                            dp(jj).lt.L(LNUM)%Dmin) dp(jj) = ZERO
                        ! Finally find the depth averaged velocity
                        U2D(jjj) = U2D(jjj) + wdiv( m(jj) , dp(jj) ) * c(cc)
                    enddo
                enddo
            enddo
            !-----------------------------------------------------------------!
            ! Get the depth-averaged tangential velocity                      !
            ! (normal interpolation reqd)                                     !
            !-----------------------------------------------------------------!
            if (i4.eq.0) then
                V_av = HALF * ( U2D(min(-isign,0)) * (DX_N + DX_R)            &
                     + U2D(min(isign,0)) * (DX_N - DX_R) ) / DX_N 
            else
                V_t  = HALF * ( U2D(min(-isign,0)) * (DX_N + DX_R)            &
                     + U2D(min(isign,0)) * (DX_N - DX_R) ) / DX_N 
            endif
        enddo
        !---------------------------------------------------------------------!
        ! First-order gradient used in determining vertical velocity profile  !
        !---------------------------------------------------------------------!
        DUDX = DUDX + ( V_t - V_av ) / DX_R
    endif    
enddo DIM_LOOP 
nnxp2 = mn2d(i+max(ix,0)+ix,j+max(iy,0)+iy) ; nn_22 = mn2d(i+ix,j+iy)
!-----------------------------------------------------------------------------!
! Include correction for large RMSE between depth-averaged and not            !
! for no fluctaution gradient                William Oct 27 2016              ! 
!-----------------------------------------------------------------------------!
if (bound(1:1).eq.'N') then
   RMSEn = ZERO; RMSEt = ZERO; count = 0
   do k = ks,ke-1
        nn_1 = mn(i,k,j) ; nn_2 = mn(i+ix,k,j+iy)
        if (nff(nn_1)%b.eq.0) cycle
        if (nff(nn_1)%b.ne.2) nff(nn_1)%b = 2
        if (nff(nn_2)%f.eq.-1) cycle
        nnxm = mn(i+max(ix,0),k,j+max(iy,0))
        if (a(abs(iy),nnxm).eq.ZERO) cycle
        if (x(2,k).ge.wll + h_in) cycle
        nnxp = mn(i+max(ix,0)+ix,k,j+max(iy,0)+iy)
        count = count + 1
        RMSEn = RMSEn + ( un(abs(iy),nnxp) - U_3D(abs(iy),nnxp2) )**2 
        if (DIM.gt.1) then
            RMSEt = RMSEt + ( un(abs(ix),nn_2) - U_3D(abs(ix),nn_22) )**2
        endif
   enddo
   if (count > 0) then
      RMSEn = sqrt(RMSEn/count)
      Muln  = max( ONE, RMSEn / abs(U_3D(abs(iy),nnxp2)) )    
      if (DIM.gt.1) then
          RMSEt = sqrt(RMSEt/count)
          Mult  = max( ONE, RMSEt / abs(U_3D(abs(ix),nn_22)) )    
      endif 
   else
      Muln  = ONE; Mult = ONE
   endif
endif
!-----------------------------------------------------------------------------!
! Loop over the vertical direction and input the values...                    !
!-----------------------------------------------------------------------------!
do k = ks,ke-1
    nn_1 = mn(i,k,j) ; nn_2 = mn(i+ix,k,j+iy)
    if (nff(nn_1)%b.eq.0) cycle
    if (nff(nn_1)%b.ne.2) nff(nn_1)%b = 2
    if (nff(nn_2)%f.eq.-1) cycle
    nnxm = mn(i+max(ix,0),k,j+max(iy,0))
    if (a(abs(iy),nnxm).eq.ZERO) cycle
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    !   Insert velocities at the boundary   !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Apply no gradient on F
    !f(nn_1) = f(nn_2) ; fp(nn_1) = f(nn_1)
    if (x(2,k).ge.wll + h_in) then                      
        f(nn_1) = ZERO ; fp(nn_1) = ZERO
        un(abs(iy),nnxm) = ZERO ; un(2,nn_1) = ZERO
        if (DIM.gt.1) un(abs(ix),nn_1) = ZERO
    else
        if (x(2,k+1).le.wll + h_in) then
            f(nn_1) = ONE ; fp(nn_1) = ONE    
        else
            f(nn_1)  = ( wll + h_in - x(2,k) ) / dx(2,k)
            fp(nn_1) = f(nn_1)
        endif
        ! Normal velocities
        if (bound(1:1).eq.'U') then     
            ! Uniform velocity distribution (NSW)
            !       u        =   U      
            un(abs(iy),nnxm) = U_av
        elseif (bound(1:1).eq.'Q') then 
            ! Quadratic velocity distribution (Boussinesq)
            znow = HALF * (x(2,k) + x(2,k+1)) + h_in + h_av
            !       u        =  U   + (  1/6*h^2     - 0.5(z + h)^2) d^2Udx^2
            un(abs(iy),nnxm) = U_av + (SIXTH*h_av**2 - HALF*znow**2) * D2UDX2
        elseif (bound(1:1).eq.'N') then 
            ! No gradient in fluctuation from the depth-averaged velocity
            nnxp = mn(i+max(ix,0)+ix,k,j+max(iy,0)+iy)
            nn_3 = mn(i+2*ix,k,j+2*iy)
            ! Find the difference in openness ratio between cell boundaries
            open_rat = min( wdiv( a(abs(iy),nnxp) , a(abs(iy),nnxm) )         &
                          * wdiv( min(ONE,f(nn_3)) + min(ONE,f(nn_2)) ,       &
                                  min(ONE,f(nn_2)) + f(nn_1) ) , TWO )
            !               u    =   U  +  ratio     *       u'
            if (x(2,k+1).le.wll + h_in) then
                un(abs(iy),nnxm) = U_av + open_rat * Muln                     &
                                    * ( un(abs(iy),nnxp) - U_3D(abs(iy),nnxp2))
            else
                nnxp = mn(i+max(ix,0)+ix,k-1,j+max(iy,0)+iy)
                un(abs(iy),nnxm) = U_av + open_rat * Muln                     &
                                    * ( un(abs(iy),nnxp) - U_3D(abs(iy),nnxp2))
            endif
            ! Mathematically since mean(u') = 0 this formulation is
            ! guaranteed to give the correct momentum flux
        endif
        ! Tangential velocity
        if (DIM.gt.1) then    
            if (bound(1:1).eq.'N') then 
                ! No gradient in fluctuation from the depth-averaged velocity
                nnxp = mn(i-abs(iy),k,j-abs(ix))
                nn_3 = mn(i+ix-abs(iy),k,j+iy-abs(ix))
                if (nnxp.eq.0) nnxp = nn_1 ; if (nn_3.eq.0) nn_3 = nn_2
                ! Find the diff in openness ratio between cell boundaries
                open_rat = min( wdiv( min(ONE,f(nn_3)) + min(ONE,f(nn_2)) ,   &
                                      f(nnxp) + f(nn_1) ) , TWO )
                !       v        =   V  +   ratio   *    v'    
                un(abs(ix),nn_1) = V_av + open_rat * Mult                     &
                                    * ( un(abs(ix),nn_2) - U_3D(abs(ix),nn_22))
            else
                ! Uniform velocity distribution (NSW)
                !       v        =  V
                un(abs(ix),nn_1) = V_av
            endif
        endif
        ! Vertical velocity
        if (a(2,nn_2).eq.ZERO) then
            un(2,nn_1) = ZERO
        else
            if (bound(2:2).eq.'N') then !.or.DUDX.eq.0.0d0) then 
                ! No gradient
                un(2,nn_1) = un(2,nn_2) ; u(2,nn_1) = u(2,nn_2)
            else! Linear distribution (NSW or Boussinesq)
                znow = x(2,k) + h_in - zkk
                un(2,nn_1) = - DUDX * znow 
            endif
        endif
    endif
#ifdef RAN !Set k-eps boundary condition as defined by turbulent intensity
    r(nn_1)%k  = akini !HALF * ( un(abs(iy),nnxm) * TI_b )**2  !max(,r1(nn_2)%k)
    d(nn_1)    = nut_b * nu      !wdiv( r1(nn_1)%k**2 , r1(nn_1)%e ) * cm
    r(nn_1)%e  = eini  !cm * r1(nn_1)%k**2 / d(nn_1) !max(,r1(nn_2)%e)
    r1(nn_1)%k = r(nn_1)%k ; r1(nn_1)%e  = r(nn_1)%e
#endif
enddo
endsubroutine Interp_2D_bound_to_3D    
