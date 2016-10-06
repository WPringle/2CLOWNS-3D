!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-normal.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Updates the normal slope and position based on the new F values
subroutine cal_surface_normal
    use variables, only: icon
    implicit none
    
225 call set_surface_position_and_info(1) ! 新しいfでspnを求める。
    if (icon.eq.100) then
        write(6,*) 'out suimen_adjust3_0'
        icon=0
        call set_nflag
        goto 225
    endif

!   call set_surface_normal ! spnから最小二乗法で新しいnorを求める。
    call set_surface_normal2 ! 最小二乗法を使わない場合

224 call set_surface_position_and_info(0) ! 新しいnorからspnを求め、水塊形状情報を
    if (icon.eq.100) then
!#ifdef BUGWIN
        write(6,*) 'suimen_adjust3 _ icon=',icon; !stop
!#endif
        icon=0
        call set_nflag
        goto 224
    endif
end subroutine
    
function cal_idou2(nfdb,nfux,ux,xx,yy,zz,dxx,dyy,dzz,nn1)
    use interface_list, only: plainadjust
    use variables, only: dt, ncpls0max, ncpls1max, ncpls2max, SMALL, ZERO, ONE
    use arrays, only: cpl, sfno, ncpls, pp_s
    implicit none
    real*8 :: cal_idou2
    integer, intent(in) :: nfdb,nfux,nn1
    real*8, intent(in) :: ux,xx,yy,zz,dxx,dyy,dzz
    integer:: mms1,icc,ip2,ip3,i,icpl2,nmen,iz2,iz4,iz3
    real*8, dimension(0:3) :: surff
    real*8 :: area,Vseg1,auxdt
    real*8, dimension(0:2) :: ss
    integer, dimension(0:30,0:30) :: path
    integer, dimension(0:ncpls0max,-2:ncpls2max) :: ncpl33
    real*8,allocatable,dimension(:,:) :: ppp
!    integer, dimension(0:ncpls0max,-2:ncpls2max) :: ncpl22
    real*8, dimension(6) :: axyz_sd
        
    auxdt = abs(ux)*dt
    if (auxdt.lt.SMALL) then
        cal_idou2 = ZERO
        return
    endif
    
    mms1 = SFNO(nn1)
    area = ONE
    ss(0:2)= ZERO !Change to calculate in local coordinates
    if(nfux.eq.3) then
        area=dxx*dyy;nmen=2;ss(2)=auxdt;!+zz
    else if(nfux.eq.-3) then
        area=dxx*dyy;nmen=1;ss(2)=dzz-auxdt;!+zz
    else if(nfux.eq.2) then
        area=dxx*dzz;nmen=5;ss(1)=auxdt;!yy+
    else if(nfux.eq.-2) then
        area=dxx*dzz;nmen=3;ss(1)=dyy-auxdt;!yy+
    else if(nfux.eq.1) then
        area=dzz*dyy;nmen=4;ss(0)=auxdt;!xx+
    else if(nfux.eq.-1) then
        area=dzz*dyy;nmen=6;ss(0)=dxx-auxdt;!xx+
    endif
    
    surff(0:2)=cpl(nmen,0:2)
    surff(3)=sum(surff(0:2)*ss(0:2))*(-1.0d0)
    path=0
    icc = 0
    ip3 = 0
    ip2 = ncpls(mms1,0,1) !対象断面上の点の総数
    icpl2 = min(ncpls0max,ncpls(mms1,0,0)) !menの総数

    iz2 = min(ncpls2max,max(10,ncpls(mms1,0,2)))
    iz3 = min(ncpls2max,max(10,ip2*2))
    iz4 = min(ncpls0max,max(10,icpl2*2))
    allocate(ppp(0:iz3,0:2))
    !Added Aug 08 2014 for pp in local coordinates
    !Set the origin to transfer local to global coordinates
    ppp(0,0:2) = [xx,yy,zz]
    do i = 1,ip2
        ppp(i,:) = pp_s(mms1,i,:) - ppp(0,0:2)
    enddo
    
    ! Find volume of fluid segment
    call PlainAdjust(nn1,Vseg1,surff(0:3),icpl2,iz2,ncpls(mms1,0:icpl2,-2:iz2),&
                     axyz_sd,iz3,ppp(0:iz3,0:2),iz4,iz3,ncpl33(0:iz4,-2:iz3))
!
    if (Vseg1.eq.ZERO) then
        cal_idou2 = ZERO
    else
        cal_idou2 = dsign(ONE,ux) * Vseg1 / area / dt
    endif

endfunction cal_idou2