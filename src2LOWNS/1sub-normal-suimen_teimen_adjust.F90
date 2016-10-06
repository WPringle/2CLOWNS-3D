#ifdef NORMAL
subroutine suimen_teimen_adjust(nn,nff,nffb,nfpp,vt1,vt7,nrm,spp,sppn,ippmax,pp&
                               ,fnn,iflag,Vseg,dVV,ncpl2,izz1,izz2,ncpl3,icon_1&
                                                                       ,axyz_sd)
    use interface_list, only: PlainAdjust
    use variables, only: t, n, sdum, ncpl00max, ncpl01max, ncpl02max, ic1_0,   &
                         ic2_0, ZERO, Leps, VSMALL, ONE, SMALL, HALF, LITTLE
    use arrays, only: cpl ,OBJNO, pp_0, ncpl0, ncplc
    use omp_lib
    implicit none
!------------------------- In/Out Variables ------------------------------------
integer,intent(in)                            :: nn, nff, nffb, nfpp, iflag,   &
                                                 izz1, izz2, ippmax
integer,intent(out)                           :: icon_1
integer,dimension(0:izz1,-2:izz2),intent(in)  :: ncpl2
integer,dimension(0:izz1,-2:izz2),intent(out) :: ncpl3
real*8,intent(in)                             :: fnn, dVV
real*8,intent(out)                            :: Vseg
real*8,dimension(0:2),intent(in)              :: spp, vt1, vt7
real*8,dimension(0:2),intent(out)             :: sppn
real*8,dimension(0:3),intent(inout)           :: nrm
real*8,dimension(0:ippmax,0:2),intent(inout)  :: pp
!------------------------- Temporary variables ---------------------------------
integer :: ii, jj, ib, ip0, ixn0, iz1, iz2, iz3, icon_error
integer,allocatable    :: ixd(:)
real*8  :: dis_max, dis_min, dif_min, dif_max, Vsea, f_temp, Vse
real*8,dimension(6)    :: axyz_sd 
real*8,allocatable     :: dist(:), Vse_array(:)
#ifdef BUGWIN
integer                :: myhr=1
real*8,dimension(1:60) :: leng
!myhr=omp_get_thread_num()
#endif
!-------------------------------------------------------------------------------
icon_1 = 0 !Once spn is found we set icon_1 to 2000
Vseg = ZERO
axyz_sd = ZERO  !??
f_temp = min(ONE,fnn); !Ensure that f is not beyond 1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! IFLAG=1:    Calculation for new surface (nf = 0 --> nf = 1 or 2 )            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (iflag.eq.1) then 
    !If the cell was an air cell previously we need to guess a new surface
    if (spp(0).eq.sdum.or.nfpp.eq.0.or.sum(nrm(:)**2).eq.ZERO) then                    
        call surface_normal_new(vt1,vt7,nrm(0:3),sppn(0:2),fnn,nff,nffb)
        icon_1 = 2000 ; return !spn found
    endif
endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   When F is very large we need to set the surface to the grid edge           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   added Pringle Sep 08 2014, updated Oct 06 2014
if (fnn.ge.Leps) then 
    nrm(0:2) = ZERO
    if (nffb.ne.0) then
        nrm(abs(nffb)-1) = real ( sign(1,-nffb) )
    else
        !Guess nfb = 2 direction
        nrm(2) = ONE
    endif
endif
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                Calculation where there is no object                          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#ifdef OBJECT
OBJ_IF: if ( OBJNO(nn).eq.0                                                    &
#if defined(HENDO) || defined(DRI)
    .and.DRINO(nn).eq.0                                                        &
#endif
        .or.fnn.ge.Leps) then   
#endif
        !Calculate spn from normal suimen_teimen equation            
        call suimen_teimen(nn,vt1,vt7,nrm(0:2),sppn(0:2),f_temp,1)
        !Get d in normal equation from calculated spn
        nrm(3) = -sum( sppn(:) * nrm(0:2) )
        !If iflag = 1: Surface normal is found and we can return 
        if (iflag.eq.1) then
            icon_1 = 2000; return
        endif
        !If iflag = 0: We need to check the water body shape too
        !Calculate the volume of water in the cell from current normal shape
        iz1 = ncpl2(0,0) ; iz2 = ncpl2(0,2) ; iz3 = ncpl2(0,1) * 2
        call PlainAdjust(nn,Vseg,nrm(0:3),iz1,iz2,ncpl2(0:iz1,-2:iz2),axyz_sd,&
                         iz3,pp(0:iz3,0:2),izz1,izz2,ncpl3(0:izz1,-2:izz2))
        Vsea = abs(Vseg_error())
        if (Vsea.lt.SMALL.and.Vseg.gt.ZERO) then
            return !Calculation is OK, lets return
        else
            !Vseg is zero (bad) or the error of Vseg is too large
            if ( Vseg.eq.ZERO ) return
            !Skip to iteration scheme below used in object calculation
        endif
#ifdef OBJECT
else
    iz1 = ncpl0(OBJNO(nn),0,0) 
    iz2 = ncpl0(OBJNO(nn),0,2) 
    iz3 = min(ncpl0(OBJNO(nn),0,1) *2,ippmax)        
endif OBJ_IF
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                Calculation where there is an object                         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! For iflag = 1 we need to calculate normal in local coordinates
if (iflag.eq.1) nrm(3) = -sum( (spp(:) - pp(0,0:2)) * nrm(0:2) )
ip0 = ncpl2(0,1) !Total number of vertices considering the object
allocate(dist(1:ip0),ixd(1:ip0))
do ii = 1,ip0
    dist(ii) = sum( nrm(0:2) * pp(ii,0:2) ) + nrm(3)
enddo
do ii = 1,ip0
    if(dist(ii).eq.sdum) cycle
    do jj = 1,ip0
        if(ii.eq.jj) cycle
        if (abs(dist(ii)-dist(jj)).lt.VSMALL) then
            dist(jj) = sdum
        endif
    enddo
enddo
call narabekae(dist(1:ip0),ip0,ixd(1:ip0),ixn0,sdum)

if (ixn0.gt.1) then
    allocate(Vse_array(ixn0))
    !===================== 2014/08/13 Following tanaka's version ==============
    do ii = 1,ixn0           
        nrm(3) = -sum( (pp(ixd(ii),:) ) * nrm(0:2) )
        call PlainAdjust(nn,Vseg,nrm(0:3),iz1,iz2,ncpl2(0:iz1,-2:iz2),axyz_sd,&
                         iz3,pp(0:iz3,0:2),izz1,izz2,ncpl3(0:izz1,-2:izz2)) 
        Vse = Vseg_error(); Vsea = abs(Vse) 
        if ( Vsea.lt.SMALL.and.Vseg.gt.ZERO ) return !Calculation OK    
        if ( fnn.ge.Leps.and.Vsea.lt.LITTLE.and.Vseg.gt.ZERO ) return
        Vse_array(ii) = Vse
    enddo
    !Get the max and min d values 
    call find_max_min_t(Vse_array,ixn0,dis_max,dis_min,dif_max,dif_min)
    !Iteration using double-false position method to find Vseg
    call dble_false_positn(icon_error)
    if (icon_error.eq.0) then
        return !Iteration method successful
    else
        !Iteration unsuccessful
        write(6,*) 'nn=', nn,'fnn=', fnn
        write(6,*) 'Vseg=', Vseg, 'Vse=', Vse, 'dVV=',dVV
        write(6,*) 'dis_max=', dis_max, 'dis_min=', dis_min
        write(6,*) 'dif_max=',dif_max,'dif_min=',dif_min
        stop 'Stop in Suimen_teimen_adjust: Too many iterations!'         
    endif
else
    call vol_check        
endif
        
contains

subroutine narabekae(b,na,ia,nb,sdum)
    implicit none
    integer :: na,is1,is2
    real*8,dimension(1:na)    :: a,b
    integer,dimension(1:na)   :: ia
    real*8,intent(in)         :: sdum
    real*8                    :: s1,s2 
    integer :: ic,nb,icc

nb = 0
do ic=1,na
    if(b(ic).ne.sdum) then
        nb=nb+1
        ia(nb)=ic
    a(nb)=b(ic)
    endif
enddo

do
    icc=0
    do ic=1,nb
        do ib=1,nb
            if (ic.ge.ib) cycle
            if (a(ic).gt.a(ib)) then
                icc=icc+1
                s1=a(ic);s2=a(ib)
                a(ic)=s2;a(ib)=s1
                is1=ia(ic);is2=ia(ib)
                ia(ic)=is2;ia(ib)=is1
            endif
        enddo
    enddo
    if(icc.eq.0) exit
enddo
end subroutine narabekae

real*8 function Vseg_error
    Vseg_error = Vseg / dVV - f_temp
end function

subroutine find_max_min_t(diff22,ixn0,d_max,d_min,df_max,df_min)
    implicit none
    integer,intent(in) :: ixn0
    real*8,intent(out) :: d_max,d_min,df_max,df_min
    real*8 :: diff_max,diff_min
    real*8,dimension(ixn0),intent(in) :: diff22
    integer :: dist_p,dist_m

dist_p   = 0     ; dist_m   = 0
diff_max = 1.0d8 ; diff_min = 1.0d8
do ii = 1,ixn0
    if (diff22(ii).gt.ZERO) then
        if (diff22(ii).lt.diff_max) then
            diff_max = diff22(ii)
            dist_p   = ii
        endif
    elseif (diff22(ii).lt.ZERO) then
        if (abs(diff22(ii)).lt.abs(diff_min)) then
            diff_min = diff22(ii)
            dist_m   = ii
        endif
    endif
enddo

if(dist_m.eq.0.or.dist_p.eq.0) then
    write(6,*) 'nn',nn,fnn
    stop 'mondaiari in find_max_min_t calc : Line 208'
endif

d_max  = -sum( (pp(ixd(dist_p),:) ) * nrm(0:2) )
d_min  = -sum( (pp(ixd(dist_m),:) ) * nrm(0:2) )
df_max = diff22(dist_p)
df_min = diff22(dist_m)
!
end subroutine find_max_min_t

subroutine dble_false_positn(icon)
    integer,intent(out) :: icon
    integer :: counter, iter
    real*8 :: eps
    counter = 0; iter = 0; icon = 0
    eps = SMALL + 2.0d0**-53 * max(abs(dis_max),abs(dis_min),ONE)
    do while (iter < 100000) !Loop iterating until suitable d value is found
        iter = iter + 1
        !Iteration using double-false position method
        nrm(3) = dis_max - dif_max * ( dis_max - dis_min )                     &
                                   / ( dif_max - dif_min ) 
        call PlainAdjust(nn,Vseg,nrm(0:3),iz1,iz2,ncpl2(0:iz1,-2:iz2),axyz_sd, &
                         iz3,pp(0:iz3,0:2),izz1,izz2,ncpl3(0:izz1,-2:izz2))
        Vse = Vseg_error(); Vsea = abs(Vse)
        if ((Vsea.lt.eps.or.abs(dis_max-dis_min).lt.VSMALL)                    &
            .and.ncpl3(ncpl3(0,0),0).gt.0.and.Vseg.gt.ZERO) then
            if (iter > 1000) write(6,*) 'iter=',iter,'nn=',nn,'Vse=',Vse
            return !Calculation OK
        elseif (Vse * dif_max.gt.ZERO) then
            !We need to update the upper bounds and reiterate
            dis_max = nrm(3)
            dif_max = Vse       
            counter = counter + 1
            !Adjustment using Illinois algorithm
            if (counter.eq.2) then
                dif_min = HALF * dif_min; counter = 0
            endif
        elseif (Vse * dif_min.gt.ZERO) then
            !We need to update the lower bounds and reiterate
            dis_min = nrm(3)
            dif_min = Vse   
            counter = counter - 1
            !Adjustment using Illinois algorithm
            if (counter.eq.-2) then
                dif_max = HALF * dif_max; counter = 0
            endif
        else
            stop 'Stop in dble_false_positn: Bad IF statement' 
        endif        
    enddo
    icon = 1
endsubroutine dble_false_positn

!The following subroutines are used for checking errors the may arise...
subroutine checkdata
write(611,*) nn,ncpl2(0,1)
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

subroutine checkdata2
write(613,*) nn,ncpl3(0,1)
do ib=1, ncpl3(0,1)
write(613,11) pp(ib,:)
11 format(1h ,e12.5,10(' ',e12.5))
continue
enddo
write(613,*) ncpl3(0,0)
do ib=1, ncpl3(0,0)
write(613,12) ncpl3(ib,-2:ncpl3(ib,0))
12 format(1h ,i10,20(' ',i2))
enddo
write(613,*) 'shape and face data end'
end subroutine 

subroutine vol_check
    use interface_list, only: crossV
    use arrays, only: in, kn, jn
    implicit none
    integer :: i,ib1,ib2,nnp
    real*8 :: VV
    real*8,dimension(0:2)  :: tmpV,tmpV2
    real*8,dimension(1:60) :: s
    VV=ZERO;s=ZERO
do i=1,ncpl2(0,0)
    nnp=ncpl2(i,0)
    if(nnp.eq.0) cycle
    tmpV=ZERO
    do ib=2,nnp-1
        ib1=ncpl2(i,ib)
        ib2=ncpl2(i,ib+1)
        tmpV = crossV(pp(ib1,0:2)-pp(ncpl2(i,1),0:2),    &
                      pp(ib2,0:2)-pp(ncpl2(i,1),0:2))+tmpV
    enddo
    s(i)=sqrt(sum(tmpV**2))*HALF
#ifdef BUGWIN
    leng(i)=sum(tmpV*pp(ncpl2(i,1),0:2))/6.0d0
#endif
    if(ncpl2(i,-1).eq.0) then
        ii=ncpl2(i,-2)
        continue
    else if(ncpl2(i,-1).eq.1) then
        ii=ncpl2(i,-2)
        continue
    else
        continue
    endif
    VV=VV+sum(tmpV*pp(ncpl2(i,1),0:2))/6.0d0
enddo

if (abs(VV-dVV).gt.1.0d-8) then
    tmpV(0:2)  = [ minval(pp(1:ncpl2(0,1),0)),&
                   minval(pp(1:ncpl2(0,1),1)),&
                   minval(pp(1:ncpl2(0,1),2)) ]
    tmpV2(0:2) = [ maxval(pp(1:ncpl2(0,1),0)),&
                   maxval(pp(1:ncpl2(0,1),1)),&
                   maxval(pp(1:ncpl2(0,1),2)) ]
    write(63,*) 'nn=,',nn,',i,k,j=,',in(nn),',',kn(nn),',',jn(nn)
    write(63,*) vv,',',dvv
    do i=1,ncpl2(0,1)
        write(63,661) i,(pp(i,0:2)-tmpV(0:2))/(tmpV2(0:2)-tmpV(0:2))
    enddo
    do i=1,ncpl2(0,0)
        write(63,662) i,ncpl2(i,0:ncpl2(i,0))
    enddo
    do i=1,ncpl2(0,0)
        write(63,663) i,cpl(i,:)
    enddo
    do i=1,ncpl2(0,0)
        write(63,664) i,s(i)
    enddo
    write(6,*) 'adjust3 taiseki'           
        call checkdata 
        call checkdata2 
    stop 'Error in suimen_teimen_adjust: Vol_check'
endif
661 format(1h ,i2,3(',',e14.7))
662 format(1h ,i2,8(',',i2))
663 format(1h ,i2,4(',',e14.7))
664 format(1h ,i2,(',',e14.7))
end subroutine vol_check

end subroutine suimen_teimen_adjust  
    
subroutine surface_normal_new(x,dx,nor1,sp1,ff,nff,nfb1)
    use variables, only: ONE, HALF, ZERO, TINY
    implicit none
    integer,intent(in) :: nff,nfb1
    real*8,intent(in)  :: ff
    real*8,dimension(0:2),intent(in)  :: x, dx
    real*8,dimension(0:2),intent(out) :: sp1
    real*8,dimension(0:3),intent(out) :: nor1
    integer :: nfbb, anfbb
    real*8  :: fff

    if (nff.eq.2) then
        nfbb = nfb1  
    elseif (nff.eq.1) then
        nfbb = -3
    elseif (nff.eq.0) then
        nfbb = -3
    else
        stop 'Error in surface_normal_new: nf = -1'
    endif

    anfbb = abs(nfbb)
    if (nfbb.lt.0) fff = ff
    if (nfbb.gt.0) fff = ONE - ff
    nor1            = ZERO 
    nor1(anfbb - 1) = ONE * sign(1,-nfbb)
    sp1(0:2)        = x(0:2) + dx(0:2) * HALF
    sp1(anfbb - 1)  = x(anfbb - 1) + dx(anfbb - 1) * fff
    where (abs(sp1(:)).lt.TINY) sp1(:) = ZERO
    !Evaluate d in the normal equation
    nor1(3) = -sum( sp1(0:2) * nor1(0:2) )
    
end subroutine surface_normal_new
#endif