#ifdef NORMAL
!%%%%%%%%%%%% FILE: 1sub-normal-set_surface_position_and_info.F90 %%%%%%%%%%%%%
! Gets the free surface information such as normals and plane equations
subroutine set_surface_position_and_info(iflag)
! iflag = -1 水面位置nor(:,0:3)を用いて、fを求める。（初期状態）
! iflag =  0 nor(:,0:2)とｆを用いて、nor(:,3)とspnを求める。水塊形状も求める。
! iflag =  1 nor(:,0:3)とｆを用いてて、nor(:,3)とspnを求める。水塊形状は不要。
use interface_list, only: centroid2, PlainAdjust
use variables, only: inns, inne, isfn, isfn2, ncpl00max, ncpl01max, mms, mms2,&
                     t, n, ie, ncpl02max, ncpls0max, ncpls1max, ncpls2max,    &
                     icon, ZERO, je, sdum, VSMALL, HALF, ONE, TINY
use arrays, only: ncpl0, ncpls, nds, pp_s , pp_0 , sfno ,sfnoi, sfno2i, surfp,&
                  sp, spn, dx, f, fp, fb, in, jn, kn, mn, objno, x, nor, sfn, &
                  sfno2, ncplc, nff
use omp_lib
implicit none
integer,intent(in) :: iflag
integer :: sig_limit = 10 !Limit of number of sig. figures in calculation of spn
integer :: nnf
integer :: isfn0
integer :: myhr ! スレッド番号
integer :: i,nn,icon_1,icc,ii,ia,ib,ip1,ip2,ic,jc,kc,ic1,ic2,ic3
integer,allocatable,dimension(:,:) :: ncpl2,ncpl3 ! セル形状情報(作業用)
real*8,allocatable :: pp(:,:) !　セル内形状の頂点
real*8,allocatable :: surfp_tmp(:,:)
real*8,dimension(0:2) :: dxx1 !セルのサイズ（DX,DY,DZ）
real*8,dimension(0:2) :: vt1,vt7  !　セル頂点１と 7の座標
real*8 :: dVV !セル体積
real*8 :: Vseg ! セル内空隙体積
real*8 :: dz, zb !vertical cell size and bottom
!real*8 :: trunc
real*8,dimension(6) :: axyz_sd ! ダミー変数
real*8,dimension(0:2) :: spnn ! セル内水面の図心（作業用）
!
! Set number of surface cells to loop over
if (iflag.eq.0) then
    isfn0 = isfn2
else
    isfn0 = isfn
endif
if (iflag.eq.-1) mms2 = mms

! Reallocate mms based arrays
if (isfn0.ne.mms) then
    mms = isfn0
    if (allocated(sfnoi)) deallocate(sfnoi)
    if (allocated(ncpls)) deallocate(ncpls)
    if (allocated(pp_s))  deallocate(pp_s)
    allocate(ncpls(1:mms,0:ncpls0max,-2:ncpls2max))
    allocate(pp_s(1:mms,0:ncpls1max,0:2))
    allocate(sfnoi(1:mms))
    mms2 = mms
    if (allocated(sfno2i)) deallocate(sfno2i)
    if (allocated(nds))    deallocate(nds)
    if (allocated(surfp))  deallocate(surfp)
    ! Guess amount of cells with extra surfaces
    allocate(sfno2i(mms + max(ie,je)))
    allocate(nds(mms + max(ie,je)))
    allocate(surfp(mms + max(ie,je),ncpls2max,0:2))  
endif
!$omp parallel default(none)                                                  &
!$omp shared(inns,inne,ncplc,t,f,fb,fp,x,dx,iflag,in,jn,kn,isfn0,mms,mms2     &
!$omp ,ncpl00max,ncpl01max,ncpl02max,icon,ncpls0max,ncpls1max,ncpls2max,sdum  &
!$omp ,ncpl0,nds,nff,nor,objno,pp_0,pp_s,spn,sfno,sfno2,sfno2i,sfnoi,sfn,sp   &
!$omp ,surfp,n,ncpls)                                                         &
!$omp private(dz,nn,nnf,Vseg,icc,ii,spnn,ia,dvv,pp,surfp_tmp,zb,ic1,ic2,ic3   &
!$omp ,ncpl3,ncpl2,ip1,ip2,ib,myhr,axyz_sd,dxx1,icon_1,vt1,vt7,ic,jc,kc,i)
!Reset the surface information
!$omp do
do nn = inns, inne
    sfno(nn)    = 0
    spn(nn,0:2) = sdum
enddo
!$omp end do
allocate(ncpl2(0:ncpls0max,-2:ncpls2max),ncpl3(0:ncpls0max,-2:ncpls2max))
allocate(pp(0:ncpls1max,0:2)) ; ncpl2 = 0 ; ncpl3 = 0
!
!===============================================================================
!===============================================================================
!       Loop over all surface cells and get the new surface information        ! 
!===============================================================================
!===============================================================================
!$omp do
SURF_CELL_ALL: do nnf = 1,isfn0   !水面セルのみを対象とする場合
    nn     = sfn(nnf)
    icon_1 = 0 
    ! 水面情報の登録
    sfnoi ( nnf ) = nn         ; sfno2i( nnf ) = 0
    sfno( nn )   = nnf
    ic  = in(nn) ; jc = jn(nn) ; kc = kn(nn)
    dVV = dx(0,ic) * dx(1,jc) * dx(2,kc) * fb(nn)    
!===============================================================================
!All iflag : Get the vertex information for the cell 
!=============================================================================== 
    call get_vertex_values(pp(0:8,0:2),dxx1,nn)
#ifdef OBJECT
#if defined(HENDO) || defined(DRI)
    if (DRINO(nn).ne.0) then
        ncpl2(0:ncpl00max,-2:mcipmx) = ncpl1(DRINO(nn),0:ncpl00max,-2:mcipmx)
        do i = 1, ncpl2(0,1)                  !Adjusting to local coordinates 
            pp(i,0:2) = pp_1(DRINO(nn),i,0:2) - pp(0,0:2)
        enddo
    else &                                      
#endif
    if (OBJNO(nn).ne.0) then
        ic1 = ncpl0(OBJNO(nn),0,0)
        ic2 = ncpl0(OBJNO(nn),0,2)
        ic3 = ncpls1max
        ncpl2(0:ic1,-2:ic2) = ncpl0(OBJNO(nn),0:ic1,-2:ic2)
        do i = 1 , ncpl2(0,1)                 !Adjusting to local coordinates 
            pp(i,0:2) = pp_0(OBJNO(nn),i,0:2) - pp(0,0:2) 
        enddo
    else
        ncpl2(0:6,-2:4) = ncplc(0:6,-2:4) 
        ic1 = 8 ; ic2=4 ; ic3 = ncpls1max
    endif                                       
#endif       
!===============================================================================
!iflag = 1 : Check to see if we can skip calculation of the new spn...
!===============================================================================
    if (iflag.eq.1) then
        if (abs(f(nn)-fp(nn)).lt.VSMALL.and.                                   &
            f(nn).ne.ZERO.and.sp(nn,0).ne.sdum) then
            !If there is little change in the F value...
            spn(nn,:)          = sp(nn,:) - pp(0,:) 
            ! Keep in local coordinates but adjust for bottom height
            if (fb(nn).ne.ONE                                                  &
#ifdef OBJECT
                .and.OBJNO(nn).eq.0                                            &
#endif                
                ) then
                spn(nn,2) = spn(nn,2) + dx(2,kc) * ( ONE - fb(nn) )
            endif
            cycle                         !(Initial spn is in local coordinates)
        endif  
    endif    
!===============================================================================
!水面位置が与えられる場合(iflag=-1) : Calculate F from Vseg (from current sp)
!===============================================================================
    if (iflag.eq.-1) then
        ! Recalculate in local coordinates as equired
        nor(nn,3) = - sum( (sp(nn,0:2) - pp(0,0:2)) * nor(nn,0:2) ) 
        call PlainAdjust(nn,Vseg,nor(nn,0:3),ic1,ic2,ncpl2(0:ic1,-2:ic2),      &
                         axyz_sd,ic3,pp(0:ic3,0:2),ncpls0max,ncpls2max,ncpl3) 
        if (fb(nn).lt.ONE) then !メッシュ内に物体が存在する場合、fを計算しなおす。
            f(nn)  = Vseg / dVV
            fp(nn) = f(nn)
        endif
!==============================================================================
!水面位置を調整する場合(iflag ~= -1) : Adjust the normal position
!==============================================================================
    else 

        if (fb(nn).eq.ONE                                                     &
#ifdef OBJECT
            .or.OBJNO(nn).ne.0                                                &
#endif                
            ) then
            dz = dx(2,kc)
            zb = x(2,kc)
        else
            dz = dx(2,kc) * fb(nn)
            zb = x(2,kc+1) - fb(nn) * dx(2,kc)
        endif
        vt1 = ZERO              ! Calculate in local coordinates in beginning
        vt7 = [dx(0,ic),dx(1,jc),dz] ! dxx of cell
        
        call suimen_teimen_adjust(nn,nff(nn)%f,nff(nn)%b,nff(nn)%fp,vt1,vt7,  &
                                  nor(nn,:),sp(nn,:),spn(nn,:),ncpls1max,pp,  &
                                  f(nn),iflag,Vseg,dVV,ncpl2,ncpls0max,       &
                                  ncpls2max,ncpl3,icon_1,axyz_sd)  
    endif
    !Error: Vseg is zero
    if (icon_1.eq.0.and.Vseg.eq.ZERO) then
        if (iflag.ne.-1) then
            write(6,*) 'Vseg=0.0,nn,t=',nn,t,f(nn),                           &
                       'n,iflag=',n,iflag,icon_1,                             &
                       'nf,nfb=',nff(nn)%f,nff(nn)%b,                         &
                       'ikj=',ic,kc,jc,nor(nn,:),sp(nn,:),spn(nn,:)
        endif
        cycle
    endif

    if (icon_1.eq.2000) then
        ! Keep in local coordinates but adjust for bottom height
        if (fb(nn).ne.ONE                                                     &
#ifdef OBJECT
            .and.OBJNO(nn).eq.0                                               &
#endif                
            ) then
            spn(nn,2) = spn(nn,2) + dx(2,kc) * ( ONE - fb(nn) )
        endif
        cycle !Cycle here if spn has already been found
    endif
!===============================================================================
!All iflag : Evaluate spn from found pp vertices
!=============================================================================== 
    ii   = 0     ! セル内水面数をカウント
    spnn = ZERO ! Set the temporary spn equal to zero
    do i = 1,ncpl3(0,0)
        if (ncpl3(i,-1).ne.2) cycle
        ii = ii + 1
        if (iflag.ne.1) then
            if (ii.eq.1) then
                ! First (and usually) only surface in cell
                do ia = 1,ncpl3(i,0)
                    !Making new free surface pp and adding on minimum values to 
                    !transfer from local to global coordinates
                    surfp(nnf,ia,0:2) = pp(ncpl3(i,ia),0:2) + pp(0,0:2) 
                enddo
                nds(nnf)    = ncpl3(i,0) ! 水面上の頂点数
                sfno2i(nnf) = nn         ! 水面が含まれるセル番号を登録。
            else
!$omp critical
                write(6,*)  'We have an extra surface in the cell..',mms2,ii,nn
                mms2  = mms2 + 1
                do ia = 1,ncpl3(i,0)
                    !Making new free surface pp and adding on minimum values to 
                    !transfer from local to global coordinates
                    surfp(mms2,ia,0:2) = pp(ncpl3(i,ia),0:2) + pp(0,0:2) 
                enddo
                nds(mms2)    = ncpl3(i,0) ! 水面上の頂点数
                sfno2i(mms2) = nn         ! 水面が含まれるセル番号を登録。
!$omp end critical
            endif
        endif
        !For all iflag... calculate new spn
        if (ncpl3(i,0).ge.3) then 
            allocate(surfp_tmp(1:ncpl3(i,0),0:2)) 
            surfp_tmp(1:ncpl3(i,0),0:2) = pp(ncpl3(i,1:ncpl3(i,0)),0:2)
            spnn(0:2) = centroid2( surfp_tmp(1:ncpl3(i,0),0:2) , ncpl3(i,0) )   &
                        + spnn(0:2)
            deallocate(surfp_tmp)
        else if(ncpl3(i,0).eq.2) then
            spnn(0:2) = ( pp(1,0:2) + pp(2,0:2) ) * HALF + spnn(0:2)
        else if(ncpl3(i,0).eq.1) then
            spnn(0:2) = pp(1,0:2) + spnn(0:2)
        else
            continue
        endif
    enddo 
    
    if (ii.ne.0) then
        ! Evaluate spn
        do i = 0,2
            spn(nn,i) = spnn(i) / real(ii)
        enddo
        !
        if (iflag.eq.1) then
            ! If iflag = 1 keep in local cooords but adjust for bottom height
            if (fb(nn).ne.ONE                                                 &
#ifdef OBJECT
                .and.OBJNO(nn).eq.0                                           &
#endif                
                ) then
                spn(nn,2) = spn(nn,2) + dx(2,kc) * ( ONE - fb(nn) )
            endif
        else
            ! If iflag = -1 or 0 revert back to global coordinates 
            spn(nn,:) = spn(nn,:) + pp(0,:)
        endif
        !
        where (abs(spn(nn,:)).lt.TINY) spn(nn,:) = ZERO
        !
        nor(nn,3) = -sum( spn(nn,:) * nor(nn,0:2) )
!
!------- Cycle here for iflag = 1-----------------------------------------------
        if (iflag.eq.1) cycle 
!===============================================================================
! IFLAG = -1 or 0 : Update various arrays
!===============================================================================
        !格子情報の登録
        icc = 0
        do ia = 1,ncpl3(0,0)
            if (ncpl3(ia,0).ne.0) then
                icc = icc + 1
                do ii = -2,ncpl3(ia,0)
                    ncpls(nnf,icc,ii) = ncpl3(ia,ii)
                enddo
            else
                continue
            endif
        enddo
        ncpls(nnf,0,0) = icc          !最終的な面の数
        ip1 = 1 ; ip2 = ncpl3(0,1)
    239 do i = ip1,ip2
            icc = 0
            do ia = 1,ncpls(nnf,0,0)
                do ib = 1,ncpls(nnf,ia,0)
                    if (i.eq.ncpls(nnf,ia,ib)) then
                        icc = icc + 1
                    endif
                enddo
            enddo
            if (icc.eq.0) exit
        enddo
        if (icc.eq.0) then
            do ia = 1,ncpls(nnf,0,0)
                do ib = 1,ncpls(nnf,ia,0)
                    if (ncpls(nnf,ia,ib).gt.i) then
                        ncpls(nnf,ia,ib) = ncpls(nnf,ia,ib) - 1
                    endif
                enddo
            enddo
            pp(i:ip2-1,0:2)   = pp(i+1:ip2,0:2)
            ip2 = ip2-1 ; ip1 = i
            goto 239
        endif
        ncpls(nnf,0,1) = ip2           !最終的な頂点の数              
        do i = 1,ip2               !Add on the origin to map to global coordinates
            pp_s(nnf,i,0:2) = pp(i,0:2) + pp(0,0:2)
        enddo
    else
        !Error: No information
        write(6,*) 'nashi, nn=',nn,ii,ncpl3(0,0),n,Vseg,&
                                ic,kc,jc,f(nn),fb(nn),icon_1
        icon       = 100
    endif
enddo SURF_CELL_ALL
!$omp end do
!===============================================================================
!Deallocate arrays
deallocate(pp)
deallocate(ncpl2)
deallocate(ncpl3)
!$omp end parallel
#ifdef TCHK
if (myhr.eq.0) then
    call LapTime(CT_temp,CT_spn,ET_spn)
endif
#endif
contains

subroutine checkdata
write(611,*) nn,ncpl2(0,1) !,objno(nn)
do ib=1, ncpl2(0,1)
write(611,11) pp(ib,:)
11 format(1h ,e12.5,10(' ',e12.5))
continue
enddo
write(611,*) ncpl2(0,0)
do ib=1, ncpl2(0,0)
write(611,12) ncpl2(ib,-2:ncpl2(ib,0))
12 format(1h ,i10,12(' ',i2))
enddo
write(611,*) 'shape data end'
end subroutine

end subroutine set_surface_position_and_info    
#endif

