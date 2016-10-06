subroutine set_object
use variables,only:npln,ndri,ndri_face,ipl0,mm0
use arrays,only:side,path0,cpl,cplc,ipl,idri_kdc,idri_fc_n,idri_fc_all_n,upl,upl0,side
implicit none
!integer::i,ia,kdc

call set_segment
call set_cpl_upl

#ifdef DDDDDD
if(mm0.gt.0) then
    continue
    if(ndri.gt.0) then
        continue  ! ínå`Ç†ÇËÅ{í«â¡ï®ëÃÇ†ÇË
    else
        continue  ! ínå`Ç†ÇËÅ{í«â¡ï®ëÃÇ»Çµ
        ndri_face=6
        npln=ndri_face+ipl0
        allocate(cpl(0:npln,0:3));cpl(1:6,0:3)=cplc(1:6,0:3)   
        allocate(idri_kdc(1:npln),idri_fc_n(1:npln),idri_fc_all_n(0:ndri,1:ipl0))        
        do ia=1,ipl0
            ndri_face=ndri_face+1
            idri_kdc(ndri_face)=0      !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®î‘çÜÇÃä÷åW
            idri_fc_n(ndri_face)=ia     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®Ç≤Ç∆ÇÃñ î‘çÜÇÃä÷åW
            idri_fc_all_n(0,ia)=ndri_face !
            cpl(ndri_face,:)=(-1.0d0)*upl0(ia,:) !ïYó¨ï®ì‡ë§å¸Ç´ÉxÉNÉgÉã
        enddo
    endif
else
    if(ndri.gt.0) then
        continue  ! ínå`Ç»ÇµÅ{í«â¡ï®ëÃÇ†ÇË
        npln=6+sum(ipl(1:ndri)) 
        allocate(upl(1:ndri,1:maxval(ipl(1:ndri)),1:4)); upl=0.0d0
        do kdc=1,ndri
            call UNV(kdc)
        enddo
        ndri_face=6
        allocate(cpl(0:npln,0:3));cpl(1:6,0:3)=cplc(1:6,0:3)
        ia=maxval(ipl(1:ndri))
        allocate(idri_kdc(1:npln),idri_fc_n(1:npln),idri_fc_all_n(1:ndri,ia))
        
        continue
        do kdc=1,ndri
            do ia=1,ipl(kdc)
                ndri_face=ndri_face+1
                idri_kdc(ndri_face)=kdc     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®î‘çÜÇÃä÷åW
                idri_fc_n(ndri_face)=ia     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®Ç≤Ç∆ÇÃñ î‘çÜÇÃä÷åW
                idri_fc_all_n(kdc,ia)=ndri_face
                cpl(ndri_face,:)=(-1.0d0)*upl(kdc,ia,:) !ïYó¨ï®ì‡ë§å¸Ç´ÉxÉNÉgÉã
            enddo
        enddo
    else
        continue  ! ínå`Ç»ÇµÅ{í«â¡ï®ëÃÇ»Çµ
        npln=6+6
        allocate(cpl(0:npln,0:3));cpl(1:6,0:3)=cplc(1:6,0:3)          
    endif
endif
!


 if(ndri.eq.0) then
    if(mm0.eq.0) then

    else
!        allocate(ipl(0));ipl(0)=ipl0
        ndri_face=6
        npln=ndri_face+ipl0
        allocate(cpl(0:npln,0:3));cpl(1:6,0:3)=cplc(1:6,0:3)   
        allocate(idri_kdc(1:npln),idri_fc_n(1:npln),idri_fc_all_n(0:ndri,1:ipl0))        
        do ia=1,ipl0
            ndri_face=ndri_face+1
            idri_kdc(ndri_face)=0      !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®î‘çÜÇÃä÷åW
            idri_fc_n(ndri_face)=ia     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®Ç≤Ç∆ÇÃñ î‘çÜÇÃä÷åW
            idri_fc_all_n(0,ia)=ndri_face !
            cpl(ndri_face,:)=(-1.0d0)*upl0(ia,:) !ïYó¨ï®ì‡ë§å¸Ç´ÉxÉNÉgÉã
        enddo
    endif
else if(ndri.eq.0) then
    npln=6+sum(ipl(0:ndri)) 
    allocate(cpl(0:npln,0:3));cpl(1:6,0:3)=cplc(1:6,0:3)
#ifndef HENDO
    if(ubound(upl,1).eq.0.and.ndri.ge.0) then
        allocate(upl(1:ndri,1:maxval(ipl(1:ndri)),1:4)); upl=0.0d0
    endif
    do kdc=1,ndri
        call UNV(kdc)
    enddo

    if(ndri.ne.0) then
        ia=max(ipl0,maxval(ipl(1:ndri)))
    else
        ia=ipl0    
    endif
    allocate(idri_kdc(1:npln),idri_fc_n(1:npln),idri_fc_all_n(0:ndri,ia))
    ndri_face=6

    do ia=1,ipl0
        ndri_face=ndri_face+1
        idri_kdc(ndri_face)=0      !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®î‘çÜÇÃä÷åW
        idri_fc_n(ndri_face)=ia     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®Ç≤Ç∆ÇÃñ î‘çÜÇÃä÷åW
        idri_fc_all_n(0,ia)=ndri_face !
        cpl(ndri_face,:)=(-1.0d0)*upl0(ia,:) !ïYó¨ï®ì‡ë§å¸Ç´ÉxÉNÉgÉã
    enddo

    do kdc=1,ndri
        do ia=1,ipl(kdc)
            ndri_face=ndri_face+1
            idri_kdc(ndri_face)=kdc     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®î‘çÜÇÃä÷åW
            idri_fc_n(ndri_face)=ia     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®Ç≤Ç∆ÇÃñ î‘çÜÇÃä÷åW
            idri_fc_all_n(kdc,ia)=ndri_face
            cpl(ndri_face,:)=(-1.0d0)*upl(kdc,ia,:) !ïYó¨ï®ì‡ë§å¸Ç´ÉxÉNÉgÉã
        enddo
    enddo
#endif
#endif
!endif

!allocate(cpl2(npln,0:3))
end subroutine  
    
subroutine set_segment
use arrays,only:cplc,ncplc,men,path0
implicit none
integer::i
allocate(path0(1:12,1:3))  
path0(1,1:2)=(/1,2/)
path0(2,1:2)=(/2,3/)
path0(3,1:2)=(/3,4/)
path0(4,1:2)=(/1,4/)
path0(5,1:2)=(/5,6/)
path0(6,1:2)=(/6,7/)
path0(7,1:2)=(/7,8/)
path0(8,1:2)=(/5,8/)
path0(9,1:2)=(/1,5/)
path0(10,1:2)=(/2,6/)
path0(11,1:2)=(/3,7/)
path0(12,1:2)=(/4,8/)

allocate(cplc(1:6,0:3))
allocate(ncplc(0:6,-2:4))
cplc(1,0:2)=(/0.0d0,0.0d0,-1.0d0/)   !z-
ncplc(1,1:4)=(/1,4,3,2/)

cplc(2,0:2)=(/0.0d0,0.0d0,1.0d0/)  !z+
ncplc(2,1:4)=(/5,6,7,8/)

cplc(3,0:2)=(/0.0d0,-1.0d0,0.0d0/)   !x-
ncplc(3,1:4)=(/1,2,6,5/)

cplc(4,0:2)=(/1.0d0,0.0d0,0.0d0/)  !y+
ncplc(4,1:4)=(/2,3,7,6/)

cplc(5,0:2)=(/0.0d0,1.0d0,0.0d0/)  !x+
ncplc(5,1:4)=(/3,4,8,7/)

cplc(6,0:2)=(/-1.0d0,0.0d0,0.0d0/)   !y-
ncplc(6,1:4)=(/1,5,8,4/)

ncplc(1:6,-1)=0         !ñ ÇÃéÌóﬁÇÕäiéqñ 
ncplc(1:6,0)=4 !äeñ ÇÃì_êîÇÕ4
ncplc(0,2)=4  !ñ ÇÃç≈ëÂì_êîÇÕ4
ncplc(0,0)=6  !ïΩñ ÇÃêîÇÕ6
ncplc(0,1)=8  !ì_ÇÃêîÇÕ8

do i=1,6
ncplc(i,-2)=i
enddo
allocate(men(1:8,1:4))
men(1,1:4)=(/1,2,3,4/)  !Å@men(äiéqñ î‘çÜ,ï¿Ç—)=(/ç\ê¨Ç∑ÇÈï”î‘çÜÇÃ
men(2,1:4)=(/5,6,7,8/)
men(3,1:4)=(/1,5,9,10/)
men(4,1:4)=(/2,6,10,11/)
men(5,1:4)=(/3,7,11,12/)
men(6,1:4)=(/4,8,9,12/)
end subroutine set_segment 
   

subroutine pre_surface
use arrays,only:nff,sfno,surfp,sfnoi,sfno2i,nds,in,jn,kn,sp,spn,pp_s,ncpls,fb,nsdata,cplusp
use variables,only:is,ie,js,je,inne,ncpls0max,ncpls1max,icon,t,ZERO&
    ,ncpl00max,ncpl01max,ncpl02max&
    ,ncpls0max,ncpls1max,ncpls2max,mmsmx,mms
      implicit none
!      integer::ic2
#ifdef BUGWIN
            integer::n2
#endif
      CALL set_nflag
      nff(:)%fp=nff(:)%f  !nfp=nf
      allocate(sfno(0:inne))
      mmsmx=(ie-is)*(je-js)*2
      ncpls0max=ncpl00max+10
      ncpls1max=ncpl01max+10
      ncpls2max=ncpl02max+10
      allocate(surfp(mmsmx,ncpls2max,0:2))
      allocate(SFNOI(mmsmx),SFNO2I(mmsmx))
      allocate(nds(mmsmx))
#ifdef BUGWIN
     n2=2488
      write(6,*) in(n2),kn(n2),jn(n2),nff(n2)%f,fb(n2),nff(n2)%b
      continue
#endif
!      ic1=ubound(ncpl0,1)

      allocate(ncpls(mmsmx,0:ncpls0max,-2:ncpls2max)); ncpls=0
      allocate(pp_s(mmsmx,0:ncpls1max,0:2)); pp_s=0.0d0
!      if(ubound(cpl,1).eq.0) then
!        allocate(cpl(0:12,0:3));cpl=ZERO
!        cpl(1:6,0:3)=cpl0(1:6,0:3)
!      endif
      
      if(t.eq.ZERO.and.nsdata.eq..false.) then
       if(cplusp.eq..false.) then  ! c+p ÇÃéû
        write(6,*) 'shoki_suimen (no nsdata)'
        call set_initial_surface_position
       endif
        write(6,*) 'suimen_adjust3'
        call set_surface_position_and_info(-1);icon=0
        write(6,*) 'set_nflag'
        CALL set_nflag
      else
  226   call set_surface_position_and_info(-1)
        if(icon.eq.100) then
          icon=0
          CALL set_nflag
          goto 226
        endif
        CALL set_nflag
      endif
      sp=spn
      write(6,*) 'end pre_surface'
end subroutine pre_surface
    

subroutine path_find_icc1_2(pp,ipp,surf,icc1,icv1,icv2,paths,spath,ncpl3,icpl2,ireturn)
use variables,only:ncpl00max,ncpl02max
use interface_list,only:crossV
implicit none

interface
    function cal_area_with_sign(nnp,icv,cpls,pp,ipp)
        integer::nnp,ipp    
        doubleprecision::cal_area_with_sign        
        doubleprecision,dimension(0:3)::cpls
        doubleprecision::pp(1:ipp,0:3)
        integer,dimension(1:nnp)::icv
    end function
end interface
integer,dimension(0:ncpl00max,-2:ncpl02max),intent(out)::ncpl3
integer,intent(out)::ireturn
doubleprecision,intent(in)::surf(0:3),pp(1:ipp,0:2)
integer,allocatable::route(:,:),routem1(:,:),routem2(:,:),route2(:,:),route2m1(:,:),icr(:),route2m2(:,:),icv4(:),ictemp(:),icv5(:),ictempc(:,:),ictempc1(:,:)
integer,intent(in)::icv1(1:icc1),icv2(1:icc1),spath,icpl2,ipp
integer::ic4,paths(1:spath,-2:3),ic1,ic2,ic5,ib,ib1,ib2,ib3,i_1,icc,ii,nnp,loop,icc1 !,ic6,ic7,ic8,ic9,ib4,ir2,irm2,ib3p
doubleprecision,allocatable::stemp(:)
doubleprecision::stempp !,tmpV(0:2),dd1

    if(icv2(1).ge.3.and.icv2(2).ge.3) then
        if(icv2(1).eq.5) then
            continue
        else if(icv2(1).eq.4.and.icv2(2).eq.4) then       
            continue
        else
            continue    
        endif
        if(icc1.eq.4) then
            continue
        endif
    
    ic4=maxval(icv2(1:icc1))

    allocate(route(0:ic4,0:spath),routem1(1:ic4,0:spath),routem2(1:ic4,0:spath))
    allocate(route2(1:ic4,0:spath),route2m1(1:ic4,0:spath),route2m2(1:ic4,0:spath))
    allocate(icv4(1:ic4),icv5(1:ic4),ictemp(1:spath),stemp(1:ic4),ictempc(1:10,1:spath),ictempc1(1:10,1:spath))
!     write(6,*) 'koko2'   
!  icv1(1) Ç©ÇÁÅ@icv1(2) Ç‹Ç≈ÇÃÉãÅ[ÉgÇêÙÇ¢èoÇ∑ÅB

if(icc1.eq.2) then
ic1=0  ;ib3=0  ;ic2=1
call find_ptp(icv1(1),icv1(2),ic2,route(1:ic4,0:spath),routem1(1:ic4,0:spath),routem2(1:ic4,0:spath))
if(ic1+ib3.lt.ic4) then
    ic2=ic1+ib3+1
    call find_ptp(icv1(2),icv1(1),ic2,route(1:ic4,0:spath),routem1(1:ic4,0:spath),routem2(1:ic4,0:spath))
endif
else if(icc1.eq.4) then
    call find_ptp4(icc1,icv1(1:icc1),route(0:ic4,0:spath),routem1(1:ic4,0:spath),routem2(1:ic4,0:spath))    
     ic1=0;ib3=route(0,0);ii=0
    do icc=1,ib3
        nnp=route(icc,0)
call cal_inner_area(nnp,route(icc,1:nnp),routem1(icc,1:nnp),routem2(icc,1:nnp),surf,pp,ipp)
         if(nnp.ne.0) then
        ii=ii+1              
        call surface_dir(nnp,route(icc,1:nnp),surf,ncpl3(ii,1:nnp),pp,ipp)
        ncpl3(ii,0)=nnp;ncpl3(0,0)=ii  
         endif
    enddo
    return
endif




ic5=0
if(ic1.le.1) then
    continue
    goto 7861
endif
    

do ib2=1,ic1
    ictemp(1:route(icv4(ib2),0))=route(icv4(ib2),1:route(icv4(ib2),0))
    icc=0
    do i_1=1,ic1
        if(i_1.eq.ib2) cycle
        ib1=route(icv4(ib2),0)
        do ib=route(icv4(i_1),0)-1,2,-1
            ib1=ib1+1
            ictemp(ib1)=route(icv4(i_1),ib)
            continue
        enddo
        stempp=cal_area_with_sign(ib1,ictemp,surf,pp,ipp)
        if(stempp.lt.1.0d-9) then
!            goto 265
        else
            icc=icc+1
            stemp(icc)=stempp
            ictempc(icc,1:2)=(/ib2,i_1/)
            continue
        endif
    enddo
    if(icc.gt.0) then
        allocate(icr(1:icc))
        icr=MINLOC(stemp(1:icc))
        ic5=ic5+1
        ictempc1(ic5,1:2)=ictempc(icr(1),1:2)
        deallocate(icr)
    endif
enddo


7861    ii=0
        icc=0
    do i_1=1,ib3
        icc=icc+1
        route2(icc,0)=route(icv5(i_1),0)-1
        route2(icc,1:route(icv5(i_1),0))=route(icv5(i_1),1:route(icv5(i_1),0))
        route2m1(icc,1:route(icv5(i_1),0))=routem1(icv5(i_1),1:route(icv5(i_1),0))
        route2m2(icc,1:route(icv5(i_1),0))=routem2(icv5(i_1),1:route(icv5(i_1),0))    
        if(sum(route2m1(icc,1:route2(icc,0))).eq.route2(icc,0)) then
            write(6,*) 'shuui ga bu ttai 0'
        endif
       
        nnp=route2(icc,0)
call cal_inner_area(nnp,route2(icc,1:nnp),route2m1(icc,1:nnp),route2m2(icc,1:nnp),surf,pp,ipp)
    if(nnp.ne.0) then    
         ii=ii+1       
        call surface_dir(nnp,route2(icc,1:nnp),surf,ncpl3(ii,1:nnp),pp,ipp)
        ncpl3(ii,0)=nnp;ncpl3(0,0)=ii  
    endif
           
    enddo
    
if(ic5.ge.1) then
do i_1=1,ic5
        ib=ictempc1(i_1,1)
        ib1=route(ib,0)-1
        route2(i_1,1:ib1)=route(ib,1:ib1)
        route2m1(i_1,1:ib1)=routem1(ib,1:ib1)
        route2m2(i_1,1:ib1)=routem2(ib,1:ib1)
        continue
        ib=ictempc1(i_1,2)
    do ib2=route(ib,0),2,-1
        ib1=ib1+1
        route2(i_1,ib1)=route(ib,ib2)
        route2m1(i_1,ib1)=routem1(ib,ib2-1)
        route2m2(i_1,ib1)=routem2(ib,ib2-1)
        route2(i_1,0)=ib1
    enddo
    
        nnp=route2(i_1,0)
    
call cal_inner_area(nnp,route2(i_1,1:nnp),route2m1(i_1,1:nnp),route2m2(i_1,1:nnp),surf,pp,ipp)
      if(nnp.ne.0) then   
        ii=ii+1       
        call surface_dir(nnp,route2(i_1,1:nnp),surf,ncpl3(ii,1:nnp),pp,ipp)
        ncpl3(ii,0)=nnp;ncpl3(0,0)=ii 
      endif  
enddo
endif
    
    ncpl3(0,0)=ii;


    !     write(6,*) 'koko3'
   deallocate(route,routem1,routem2)
   deallocate(route2,route2m1,route2m2)
   deallocate(icv4,icv5,ictemp,stemp)
   
        ireturn=2000
        return       
else
        write(6,*) 'icc1=2'
        ireturn=1000
        return
!        goto 912
!        stop
endif
    
!     write(6,*) 'koko4'
contains
subroutine find_ptp(icv_s,icv_e,ic2,rte,rtem1,rtem2)
integer::icv_s,icv_e,rte(1:ic4,0:spath),rtem1(1:ic4,0:spath),rtem2(1:ic4,0:spath),ic2
do i_1=ic2,ic4
    rte(i_1,1)=icv_s
    icc=1
    if(sum(paths(1:spath,3)).eq.0) then
        exit
    endif
    loop=0
765 do ib=1,spath
        if(paths(ib,3).ne.1) cycle
        if(paths(ib,1).eq.rte(i_1,icc)) then
            rtem1(i_1,icc)=paths(ib,-1);rtem2(i_1,icc)=paths(ib,-2)
            icc=icc+1;rte(i_1,icc)=paths(ib,2);paths(ib,3)=0;exit
        else if(paths(ib,2).eq.rte(i_1,icc)) then
            rtem1(i_1,icc)=paths(ib,-1);rtem2(i_1,icc)=paths(ib,-2)
            icc=icc+1;rte(i_1,icc)=paths(ib,1);paths(ib,3)=0;exit
        endif
    enddo
    if(rte(i_1,icc).eq.icv_e) then
        ic1=ic1+1;icv4(ic1)=i_1
        rte(i_1,0)=icc
    else if(rte(i_1,icc).eq.icv_s) then
        if(icc.ne.1) then
            rte(i_1,0)=icc;ib3=ib3+1;icv5(ib3)=i_1
        else
            rte(i_1,0)=0
        endif
    else
        loop=loop+1
        if(loop.lt.100) then
        goto 765
        else
            write(6,*) 'loop=1000 object-mod'
        continue
        endif
    endif
enddo  
end subroutine find_ptp

subroutine find_ptp4(icvc,icvv1,rte,rtem1,rtem2)
implicit none   
integer::rte(0:ic4,0:spath),rtem1(1:ic4,0:spath),rtem2(1:ic4,0:spath),icvc,icvv1(1:icvc),ia,ib,ic,iev(1:10),ievv,iev1(1:10),ievv1,iev2(1:10),irte !,ic2,icv_s,icv_e,icvv2
doubleprecision::dd1,dd2,dista
irte=0
2315 ievv=1;irte=irte+1
ic=0
do ia=1,spath
    if(paths(ia,3).eq.0) cycle
    ievv1=0
    do ib=1,icvc
    if(paths(ia,1).eq.icvv1(ib)) then
        ievv1=ievv1+1
    endif
    if(paths(ia,2).eq.icvv1(ib)) then
        ievv1=ievv1+1
    endif    
    enddo
    if(ievv1.eq.0) then
        ic=ia;exit
    endif
enddo
if(ic.eq.0) then
    continue
    rte(0,0)=irte-1
    return
endif
ievv=1;rte(irte,ievv)=paths(ic,1);paths(ic,3)=0
rtem1(irte,ievv)=paths(ic,-1);rtem2(irte,ievv)=paths(ic,-2)
ievv=2;rte(irte,ievv)=paths(ic,2)

3214 ievv1=0
do ia=1,spath
    if(paths(ia,1).eq.rte(irte,ievv).and.paths(ia,2).ne.rte(irte,ievv-1)) then
        ievv1=ievv1+1
        iev1(ievv1)=paths(ia,2);iev2(ievv1)=ia
    endif
    if(paths(ia,2).eq.rte(irte,ievv).and.paths(ia,1).ne.rte(irte,ievv-1)) then
        ievv1=ievv1+1
        iev1(ievv1)=paths(ia,1);iev2(ievv1)=ia
    endif
enddo
if(ievv1.eq.1) then
    if(paths(iev2(ievv1),3).eq.0) then
        continue
        irte=irte-1
        goto 2315
    endif
else if(ievv1.eq.2) then
    dd1=dista(pp(iev1(1),:)-pp(rte(irte,1),:))
    dd2=dista(pp(iev1(2),:)-pp(rte(irte,1),:))
    if(dd1.lt.dd2) then
        ievv1=1
    else
        ievv1=2
    endif
else if(ievv1.eq.0) then
    continue
else
    continue
endif
rtem1(irte,ievv)=paths(iev2(ievv1),-1);rtem2(irte,ievv)=paths(iev2(ievv1),-2)
ievv=ievv+1;rte(irte,ievv)=iev1(ievv1) 
if(rte(irte,ievv).eq.rte(irte,1)) then
    rte(irte,0)=ievv-1
continue
goto 2315
else
goto 3214
endif

do ia=1,spath
    if(paths(ia,1).eq.iev(ievv).and.paths(ia,2).ne.iev(ievv-1)) then
        ievv1=ievv1+1
        iev1(ievv1)=paths(ia,2);iev2(ievv1)=ia
    endif
    if(paths(ia,2).eq.iev(ievv).and.paths(ia,1).ne.iev(ievv-1)) then
        ievv1=ievv1+1
        iev1(ievv1)=paths(ia,1);iev2(ievv1)=ia
    endif
enddo
end subroutine find_ptp4
end subroutine path_find_icc1_2   
    

    
doubleprecision function cal_area_with_sign(nnp,icv,cpls,pp,ipp)
use variables, only: ZERO, HALF
use interface_list,only:crossV
implicit none
doubleprecision,dimension(0:2)::tmpV  !,tmpV1
doubleprecision,dimension(0:2),intent(in)::cpls
doubleprecision,intent(in)::pp(1:ipp,0:2)
integer,dimension(1:nnp)::icv
integer,intent(in)::nnp,ipp
integer::ib
tmpV=ZERO
do ib=2,nnp-1
!tmpV1=crossV(pp(icv(ib),0:2)-pp(icv(1),0:2),pp(icv(ib+1),0:2)-pp(icv(1),0:2))
tmpV=crossV(pp(icv(ib),0:2)-pp(icv(1),0:2),pp(icv(ib+1),0:2)-pp(icv(1),0:2))+tmpV
enddo
cal_area_with_sign=dot_product(tmpV,cpls(0:2))*HALF
end function cal_area_with_sign   


doubleprecision function surface_dir2(nnp,icv,cpls,pp,ipp)
use variables, only: ZERO, HALF
use interface_list,only:crossV
implicit none
doubleprecision,dimension(0:2)::tmpV
doubleprecision,dimension(0:2),intent(in)::cpls
doubleprecision::pp(1:ipp,0:2)
!integer,dimension(1:nnp)::ncp
integer,dimension(1:nnp)::icv
!doubleprecision::dd2,dis1 !,dis2
integer,intent(in)::nnp,ipp
integer::ib !,ic
tmpV=ZERO
do ib=2,nnp-1
tmpV=crossV(pp(icv(ib),0:2)-pp(icv(1),0:2),pp(icv(ib+1),0:2)-pp(icv(1),0:2))+tmpV
enddo
surface_dir2=dot_product(tmpV,cpls(0:2))*HALF

end function

subroutine surface_dir(nnp,icv,cpls,ncp,pp,ipp)
use variables, only: ZERO
!doubleprecision,dimension(0:2)::tmpV
implicit none
integer,intent(in)::nnp,ipp
doubleprecision::pp(1:ipp,0:2)
doubleprecision,dimension(0:2),intent(in)::cpls
integer,dimension(1:nnp)::ncp
integer,dimension(1:nnp)::icv
doubleprecision::dd2 !,dis1,dista !,dis2
real*8::surface_dir2
integer::ic !,ib
 
!tmpV=ZERO
!do ib=2,nnp-1
!tmpV=crossV(pp(icv(ib),0:2)-pp(icv(1),0:2),pp(icv(ib+1),0:2)-pp(icv(1),0:2))+tmpV
!enddo
!dd2=dot_product(tmpV,cpls(0:2))
write(6,*) 'ipp=',ipp,'nnp=',nnp,'icv(nnp-1 ; 2 ; 1)',icv(nnp-1),icv(2),icv(1)
dd2=surface_dir2(nnp,icv,cpls,pp,ipp)

 if(dd2.gt.ZERO) then
        ncp(1:nnp)=icv(1:nnp)
 else
        do ic=1,nnp
        ncp(ic)=icv(nnp-ic+1)
        enddo
 endif   
end subroutine

subroutine cal_inner_area(nnp,icv,icvm,icvn,surf,pp,ipp)
use variables, only: ZERO
use arrays,only:upl,upl0,idri_kdc,idri_fc_n
use interface_list,only:crossV
implicit none
!#include "interface.inc"
doubleprecision::pp(1:ipp,0:2),tmpV(0:2)
doubleprecision,dimension(0:3),intent(in)::surf
!integer,dimension(1:nnp)::ncp
integer,dimension(1:nnp)::icv,icvm,icvn
doubleprecision::dd1 !,dis1,dis2
integer,intent(in)::ipp
integer::nnp,ib,ib1,ib2,ib3,ib3p !,ic

     ib1=0;ib2=0;ib3=0
     do ib=1,nnp
     if(icvn(ib).gt.6) ib1=ib1+1
   !  if(icvm(ib).eq.1.and.icvn(ib).eq.0) ib2=ib2+1     !ó†ñ Çä‹ÇﬁèÍçáÅB
   !  if(icvm(ib).eq.0.and.icvn(ib).eq.2) ib3=ib3+1     !ÉZÉãè„ñ Çä‹ÇﬁèÍçáÅB    
     enddo
     
     if(ib1.eq.nnp.or.ib1+ib3.eq.nnp) then
       !  if(icv(1).ne.ncpl3(ii,1)) then
       !  icv(1:nnp)=ncpl3(ii,1:nnp)
       !  icvnt(1:nnp)=icvn(1:nnp)
       !  do ib3=1,nnp-1
       !  icvn(ib3)=icvnt(nnp-ib3)
       !  enddo
       !  endif
         dd1=10.0d0
         do ib3=1,nnp 
             if(icvn(ib3).eq.2) cycle
             ib3p=ib3+1;if(ib3.eq.nnp) ib3p=1
#ifdef OBJECT             
             if(idri_kdc(icvn(ib3)).eq.0) then
             tmpV=crossV( pp(icv(ib3p),0:2)-pp(icv(ib3),0:2) , upl0(idri_fc_n(icvn(ib3)),1:3) )                 
             else
             tmpV=crossV( pp(icv(ib3p),0:2)-pp(icv(ib3),0:2) , upl(idri_kdc(icvn(ib3)),idri_fc_n(icvn(ib3)),1:3) )
             endif
#endif             
             dd1=dot_product(surf(0:2),tmpV)
             if(dd1.lt.ZERO) exit
         enddo
        if(dd1.lt.ZERO) then 
             continue
             nnp=0
       !     cpl22(ii,:)=surf(:)*(-1.0d0)
       !     call surface_dir3(nnp,icv(1:nnp),cpl22(ii,:),ncpl3(ii,1:nnp))
        endif
     else if(ib1+ib2.eq.nnp) then
             continue
       !     cpl22(ii,:)=surf(:)*(-1.0d0)
       !     call surface_dir3(nnp,icv(1:nnp),cpl22(ii,:),ncpl3(ii,1:nnp))     
     endif
    end subroutine
    
    
    
#ifdef DDDDDDD     
subroutine naigai2(cc,ipp,pp,icc1,surf)
use interface_list,only:dista
implicit none
integer,intent(in)::ipp
doubleprecision,intent(in)::cc(0:2),pp(1:ipp,0:2),surf(0:3)
doubleprecision::pc(0:2),cc_pc(0:2),b_c(0:2),a_b(0:2),aa(1:2,1:2),det,bb(1:2,1:2),dd(1:2),d3,de(0:2)
doubleprecision::r1(0:2,0:2),r2(0:2,0:2),rr(0:2,0:2),across(0:2)
integer::icc,ib,ib1
integer,intent(out)::icc1

pc(0)=sum(pp(1:ipp,0))/dfloat(ipp)
pc(1)=sum(pp(1:ipp,1))/dfloat(ipp)
pc(2)=sum(pp(1:ipp,2))/dfloat(ipp)  
d3=sum(surf(0:2)*pc(0:2))+surf(3)
continue

if(surf(0)**2+surf(1)**2.eq.ZERO) then
rr=ZERO
rr(0,0)=1.0d0;rr(1,1)=1.0d0;rr(2,2)=1.0d0
else
r1(0,0)=surf(1)/sqrt(surf(0)**2+surf(1)**2)
r1(0,1)=(-1.0d0)*surf(0)/sqrt(surf(0)**2+surf(1)**2)
r1(0,2)=ZERO
r1(1,0)=surf(0)/sqrt(surf(0)**2+surf(1)**2)
r1(1,1)=surf(1)/sqrt(surf(0)**2+surf(1)**2)
r1(1,2)=ZERO
r1(2,0)=ZERO
r1(2,1)=ZERO
r1(2,2)=1.0d0
de=MATMUL(r1,surf(0:2))
 
r2(0,0)=1.0d0
r2(0,1)=ZERO
r2(0,2)=ZERO
r2(1,0)=ZERO
r2(1,1)=de(2)/sqrt(de(1)**2+de(2)**2)
r2(1,2)=(-1.0d0)*de(1)/sqrt(de(1)**2+de(2)**2)
r2(2,0)=ZERO
r2(2,1)=de(1)/sqrt(de(1)**2+de(2)**2)
r2(2,2)=de(2)/sqrt(de(1)**2+de(2)**2)
!de=MATMUL(r2,de(0:2))
rr= MATMUL(r2,r1)
endif
!de=MATMUL(rr,surf(0:2))
 
    icc=0;icc1=0
    cc_pc=pc-cc
    if(dista(cc_pc).lt.1.0d-10) then
        continue
        icc1=1
        return
    endif
        cc_pc(0:2)=MATMUL(rr,cc_pc(0:2))
    do ib=1,ipp
        ib1=ib+1
    if(ib.eq.ipp) ib1=1
        b_c(0:2)=pp(ib1,0:2)-pp(ib,0:2)
        a_b(0:2)=pp(ib,0:2)-cc
        
        a_b(0:2)=MATMUL(rr,a_b(0:2))
        b_c(0:2)=MATMUL(rr,b_c(0:2))
        aa(1:2,1)=cc_pc(0:1)
        aa(1:2,2)=(-1.0d0)*b_c(0:1)
        det=aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)
        if(abs(det).lt.1.0d-10) cycle
        bb(1,1)=aa(2,2);bb(2,2)=aa(1,1);bb(2,1)=aa(2,1)*(-1.0d0);bb(1,2)=aa(1,2)*(-1.0d0)
        dd=MATMUL(bb,a_b(0:1))/det  

    if(dd(2).lt.1.0d0+1.0d-9.and.dd(2).gt.-1.0d-9.and.abs(dd(1)).lt.1.0d-9) then !ÉIÉìÉâÉCÉì
            icc1=1
            return
    else if(abs(dd(2)).lt.1.0d-9.and.dd(1).ge.ZERO) then
            icc=icc+1
    else if(abs(dd(2)-1.0d0).lt.1.0d-9.and.dd(1).ge.ZERO) then
            continue
    else if(dd(2).gt.ZERO.and.dd(2).lt.1.0d0.and.dd(1).ge.ZERO) then
            icc=icc+1
    endif
        continue
    enddo     
    if(icc.eq.2.or.icc.eq.4) then
        icc1=0
    else if(icc.eq.1.or.icc.eq.3.or.icc.eq.5) then
        icc1=1
        return
    else if(icc.eq.0) then
        write(6,*) 'icc=',icc
        continue
    else
        continue
    endif
#ifdef DDDDDD    
!omp parallel do private(pxmin,pxmax,pymax,pymin,cc_pc,pc,ib,ib1,b_c,a_b,aa,det,dd,icc) schedule(dynamic) 
    do ii=1,ipl(kdc)
    if(ipo(kdc,ii).eq.0) cycle
    
pxmin=minval(vt(1,npo(1,ii,1:ipo(kdc,ii)),1));if(cc(0).lt.pxmin) cycle
pxmax=maxval(vt(1,npo(1,ii,1:ipo(kdc,ii)),1));if(cc(0).gt.pxmax) cycle
pymin=minval(vt(1,npo(1,ii,1:ipo(kdc,ii)),2));if(cc(1).lt.pymin) cycle
pymax=maxval(vt(1,npo(1,ii,1:ipo(kdc,ii)),2));if(cc(1).gt.pymax) cycle

    b_c(1:2)=vt(1,npo(kdc,ii,ib1),1:2)-vt(1,npo(kdc,ii,ib),1:2)
    a_b(1:2)=vt(1,npo(kdc,ii,ib),1:2)-cc

    det=aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)
    bb(1,1)=aa(2,2);bb(2,2)=aa(1,1);bb(2,1)=aa(2,1)*(-1.0d0);bb(1,2)=aa(1,2)*(-1.0d0)
    !ee=MATMUL(bb,aa)/det
    dd=MATMUL(bb,a_b)/det
    !ff1(1:2)=cc_pc*dd(1)
    !ff2(1:2)=a_b(1:2)+b_c(1:2)*dd(2)
    !if(sum((ff1-ff2)**2).gt.1.0d-9) then
    !    continue
    !else
        if(dd(2).le.1.0d0.and.dd(2).ge.ZERO.and.abs(dd(1)).lt.1.0d-9) then !ÉIÉìÉâÉCÉì
            zh=(-upl(1,ii,1)*cc(0)-upl(1,ii,2)*cc(1)-upl(1,ii,4))/upl(1,ii,3)
            continue
            return
        else if(dd(2).lt.1.0d0.and.dd(2).ge.ZERO.and.dd(1).ge.ZERO) then
            icc=icc+1
        endif
    !endif
    enddo

    if(icc.eq.2.or.icc.eq.4) then
        cycle
    else if(icc.eq.1.or.icc.eq.3.or.icc.eq.5) then
        zh=(-upl(1,ii,1)*cc(0)-upl(1,ii,2)*cc(1)-upl(1,ii,4))/upl(1,ii,3)
        return
    else if(icc.eq.0) then
        write(6,*) 'icc=',icc
        continue
    else
        continue
    endif
    enddo    
!omp end parallel do    
#endif
    end subroutine
#endif       
function side_num(ic1,ic2)
integer::side_num
integer,intent(in)::ic1,ic2;
side_num=0; if(ic1.gt.8.or.ic2.gt.8) return
    if(ic1.eq.1.and.ic2.eq.2.or.ic1.eq.2.and.ic2.eq.1) then
        side_num=1
    else if(ic1.eq.2.and.ic2.eq.3.or.ic1.eq.3.and.ic2.eq.2) then
        side_num=2
    else if(ic1.eq.3.and.ic2.eq.4.or.ic1.eq.4.and.ic2.eq.3) then
        side_num=3
    else if(ic1.eq.4.and.ic2.eq.1.or.ic1.eq.1.and.ic2.eq.4) then
        side_num=4
    else if(ic1.eq.5.and.ic2.eq.6.or.ic1.eq.6.and.ic2.eq.5) then    
        side_num=5
    else if(ic1.eq.6.and.ic2.eq.7.or.ic1.eq.7.and.ic2.eq.6) then
        side_num=6
    else if(ic1.eq.7.and.ic2.eq.8.or.ic1.eq.8.and.ic2.eq.7) then
        side_num=7
    else if(ic1.eq.8.and.ic2.eq.5.or.ic1.eq.5.and.ic2.eq.8) then 
        side_num=8
    else if(ic1.eq.1.and.ic2.eq.5.or.ic1.eq.5.and.ic2.eq.1) then
        side_num=9
    else if(ic1.eq.2.and.ic2.eq.6.or.ic1.eq.6.and.ic2.eq.2) then
        side_num=10
    else if(ic1.eq.3.and.ic2.eq.7.or.ic1.eq.7.and.ic2.eq.3) then
        side_num=11
    else if(ic1.eq.4.and.ic2.eq.8.or.ic1.eq.8.and.ic2.eq.4) then    
        side_num=12     
    endif
    end function   


      subroutine set_cpl_upl
      use variables,only:ndri,mm0,npln,ndri_face
      use arrays,only:cpl,cplc,upl,idri_kdc,idri_fc_n,idri_fc_all_n,ipl,rhod
      implicit none
      integer::ia,kdc

      ia=0
      if(ndri.ge.1) ia=sum(ipl(1:ndri))
      npln=6+ia
      allocate(cpl(1:npln,0:3)); cpl(1:6,0:3)=cplc(1:6,0:3)
      if(ndri.eq.0) return

      allocate(upl(1:ndri,1:ia,1:4)); upl=0.0d0

      do kdc=1,ndri
        if(rhod(kdc).gt.0.0d0) call volume(kdc)
        call UNV(kdc)
      enddo

      allocate(idri_kdc(1:npln),idri_fc_n(1:npln),idri_fc_all_n(1:ndri,ia))
      idri_kdc=0; idri_fc_n=0; idri_fc_all_n=0

      ndri_face=6
      do kdc=1,ndri
        do ia=1,ipl(kdc)
          ndri_face=ndri_face+1
          idri_kdc(ndri_face)=kdc     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®î‘çÜÇÃä÷åW
          idri_fc_n(ndri_face)=ia     !ïYó¨ï®ï\ñ í Çµî‘çÜÇ∆ïYó¨ï®Ç≤Ç∆ÇÃñ î‘çÜÇÃä÷åW
          idri_fc_all_n(kdc,ia)=ndri_face
          cpl(ndri_face,:)=(-1.0d0)*upl(kdc,ia,:) !ïYó¨ï®ì‡ë§å¸Ç´ÉxÉNÉgÉã
        enddo
      enddo

      end subroutine


subroutine volume(kdc)
use variables, only: ZERO
use arrays,only:ipl,ipo,vt,npo,upl,md,rhod,GD,ivt
use interface_list,only:crossL
implicit none
doubleprecision,dimension(3)::DL1,DL2,DL,GD0 !,dGD1=(/ZERO,ZERO,ZERO/)
integer,intent(in)::kdc
doubleprecision::S,dv
INTEGER::i,j
doubleprecision::dri_v     
call UNV                
            dri_v=ZERO        
! élñ ëÃÇÃëÃêœÇÃòaÇÊÇËÅCïYó¨ï®ÇÃëÃêœdri_vÇ∆èdêSGD0(:)ÇãÅÇﬂÇÈÅD
            GD0(:)=ZERO
            
            
do i=1,ipl(kdc)
    do j=2,ipo(kdc,i)-1
    DL1(:)=vt(kdc,npo(kdc,i,j),:)-vt(kdc,npo(kdc,i,1),:)
    DL2(:)=vt(kdc,npo(kdc,i,j+1),:)-vt(kdc,npo(kdc,i,1),:)
    S=crossL(DL1,DL2)/2.0d0
    DV=-S*upl(kdc,i,4)/3.0d0
    dri_v=dri_v+DV          
        DL(:)=(vt(kdc,npo(kdc,i,1),:)+vt(kdc,npo(kdc,i,j),:)+vt(kdc,npo(kdc,i,j+1),:))/4.0d0
        GD0(:)=GD0(:)+DL(:)*DV
    enddo
enddo
    if(rhod(kdc).gt.0.0d0) then
        md(kdc)=dri_v*rhod(kdc)
        write(6,*) 'md=',md(kdc)
    endif
          GD0(:)=GD0(:)/dri_v
          if(GD(kdc,1).eq.-1.0d0) then
             GD(kdc,0:2)=GD0(:)    
          else
            do i=1,ivt(kdc)
                vt(kdc,i,:)=vt(kdc,i,:)+(GD(kdc,0:2)-GD0(:))
                where(vt(kdc,i,:).lt.1.0d-8)   vt(kdc,i,:)=ZERO
            enddo
          endif 
          !!!!!!!!!!!!!!!!!!!
          !gd(kdc,3)=gd(kdc,3)-0.03
          !!!!!!!!!!!!!!!!!!!
  !        call rotation_and_move(dGD1,theta(kdc,:))
!#ifdef BUGWIN
        write(6,*) 'V,GD=',dri_v,GD(kdc,0:2)
      !  do i=1,ivt(kdc)
      !    write(6,*) vt(kdc,i,:)
      !  enddo
        continue
!#endif
end subroutine