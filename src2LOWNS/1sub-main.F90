!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-main.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Contains various main subroutines
subroutine set_ThreadsNumber
use variables,only:nthr
use omp_lib
integer::nt,ist

    !   only for Open MP
!$    write(*,'(/,1X,A)') 'Open MP Parallel Calculation'
!$omp parallel
!$    nthr=omp_get_num_threads()
!$omp end parallel
       read(5,*,iostat=ist ) nt       
       if(ist.eq.0.and.nt.ne.0) then
!$    write(*,'(1X,A,I3)') 'Required Threads Number=',nt
       else
!$    write(*,'(1X,A)') 'Required Threads Number=Not Set'   
           nt=nthr
       endif
!$    write(*,'(1X,A,I3)') 'Available Threads Number=',nthr
!$    nthr=min(nthr,nt)
!$    call omp_set_num_threads(NTHR)
!$    write(*,'(1X,A,I3,/)') 'Using Threads Number=',nthr
end subroutine set_ThreadsNumber
    
subroutine set_title
    use variables,only: nkai, t, tw, ko_t, ko_tmp, ko_out, testa, title,      &
                        icong, n, t_end, ZERO
    implicit none
    integer :: n1, n2, nslash, ns, ne
    character(len=255) :: num
    
    ko_t = ko_t + 1
    write(6,*) 'n=',n,'nkai,t,tw,ko_t < ko_tmp=',nkai,t,tw,ko_t,ko_tmp
    call owari(testa,n2) ; call hajime(testa,n1) ; call sura(testa,nslash)
    if (t.lt.ZERO.or.icong.gt.100) then
        title = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'.error'
        ko_t  = -1
        write(6,*) 'start error writing'
    elseif (ko_t.eq.ko_tmp.or.t.ge.t_end.or.n.eq.nkai) then
        call ttonum(t,num,ns,ne)
        title = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'.t'//num(ns:ne)
        ko_t = 0
    elseif (ko_out.eq.1) then ! sfdata を細かく出力する場合
        call ttonum(t,num,ns,ne)
        title = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'.t'//num(ns:ne)
    else
        title = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'.temp'
        write(6,*) 'start temp writing'
    endif
end subroutine set_title     
!    
subroutine set_dx
    use interface_list !, only: round
    use variables, only: ZERO, nthr, ONE
    use arrays, only: le, x, ls, dx
    implicit none
    integer :: j, i, dec
    
    do j = 0,2
        if (x(j,le(j)+1).eq.ZERO) then
            x(j,le(j)+1) = x(j,le(j)) + x(j,le(j)) - x(j,le(j)-1)
        endif       
        do dec = 0,10
            !Guess from first cell size
            if ((x(j,ls(j)+1)-x(j,ls(j)))/10d0**(4-dec).gt.ONE) exit
        enddo
        do i = ls(j)-2,le(j)
            !Use round function in Lib_x64_will to round dx to n dec. places
            dx(j,i) = round( x(j,i+1)-x(j,i) ,dec )
        enddo
        dx(j,le(j)+1) = dx(j,le(j))
    enddo
endsubroutine set_dx
!
subroutine set_p
    use variables, only: is ,js ,ks, ie, je ,ke, g, HALF
    use arrays, only: nff, sp, mn, x, f, p, dx, rhon
    implicit none
    real*8 :: hmax
    integer :: k1,i,j,k,nn
    
    k1 = -10000
    do i=is-2,ie+1
        do j=js-2,je+1
            do k=ke-1,ks,-1
                nn=mn(i,k,j)
                if(nff(nn)%f.eq.2) then
                    k1=k
#ifdef NORMAL
                    hmax=sp(nn,2)
#else             
                    hmax=x(2,k)+(x(2,k+1)-x(2,k))*f(nn)
#endif
                    goto 111
                endif
            enddo       
            if(i.ge.is.and.i.lt.ie.and.j.ge.js.and.j.lt.je.and.k1.ne.-10000) goto 111
            cycle
         
    111     do k=k1,ke
                p(mn(i,k,j))=(hmax-(x(2,k)+x(2,k+1))*HALF)&
#ifdef DAKU
                *rhon(mn(i,k,j))&
#endif             
                *g
            enddo
            do k=k1-1,ks,-1
                p(mn(i,k,j))=p(mn(i,k+1,j))+(dx(2,k+1)+dx(2,k))&
#ifdef DAKU
                *rhon(mn(i,k+1,j))&
#endif                
                *HALF*g
            enddo
        enddo
    enddo
end subroutine set_p
!
function cal_UNV(vtc,ipoc)
use variables, only: ZERO
use interface_list,only:crossL,crossV
implicit none
doubleprecision::cal_UNV(1:4)
integer,intent(in)::ipoc
doubleprecision,intent(in)::vtc(ipoc,1:3)
doubleprecision,dimension(3)::DL1,DL2 !,DL
doubleprecision::DLL,dd1 !,tmpv(0:2)
integer::ii,jj,ierror=0,loop !,ib,kdc

if(ipoc.eq.0) then
      cal_unv=ZERO
      return
endif

  DL1(:)=vtc(2,:)-vtc(1,:)
  DL2(:)=vtc(ipoc,:)-vtc(1,:)
  loop=0
238  call cross(DL1,DL2,cal_unv(1:3))
  DLL=crossL(DL1,DL2)
  cal_unv(1:3)=cal_unv(1:3)/DLL
  cal_unv(4)=-sum( cal_unv(1:3) * vtc(1,1:3) )
 
if(ipoc.eq.3) return

do jj=3,ipoc-1
    dd1=sum(cal_unv(1:3)*vtc(jj,:))+ cal_unv(4)
    if(abs(dd1).gt.1.0d-9) then
        write(6,*) 'ia=',ii,'npo=',jj,dd1
        ierror=1
        cal_unv=-10000.0d0
        !if(loop.eq.100) then
        !write(6,*) 'unv loop=',loop
        return
        !endif
    endif
enddo
#ifdef DDDDD
if(ierror.eq.1) then
 loop=loop+1
  DL1(:)=vtc(3,:)-vtc(2,:)
  DL2(:)=vtc(1,:)-vtc(2,:)
  ierror=0
  goto 238
endif
tmpV=ZERO
do ib=2,ipoc-1
    tmpV=crossV(vtc(ib,:)-vtc(1,:),vtc(ib+1,:)-vtc(1,:))+tmpV
enddo
    dll=dot_product(tmpV,cal_unv(1:3))
if(dll.lt.ZERO) then
    cal_unv(1:4)=cal_unv(1:4)*-ONE
endif
#endif

end function       
    
subroutine UNV(kdc) ! 面の単位法線ベクトルを求める
use interface_list,only:crossL,crossV
use arrays,only:upl,vt,ipl,ipo,npo
use variables,only:iplmax,ndri,ZERO,ONE
implicit none
doubleprecision,dimension(3)::DL1,DL2 !,DL
doubleprecision::DLL,tmpv(0:2),dd1
integer::ii,ib,jj,ierror=0,loop,kdc
if(.not.allocated(upl)) then
    allocate(upl(1:ndri,1:maxval(ipl(1:ndri)),1:4))
    upl=0.0d0
    continue
endif
do ii=1,ipl(kdc)
  if(ipo(kdc,ii).eq.0) cycle
  DL1(:)=vt(kdc,npo(kdc,ii,2),:)-vt(kdc,npo(kdc,ii,1),:)
  DL2(:)=vt(kdc,npo(kdc,ii,ipo(kdc,ii)),:)-vt(kdc,npo(kdc,ii,1),:)
  loop=0
238  call cross(DL1,DL2,upl(kdc,ii,1:3))
  DLL=crossL(DL1,DL2)
  upl(kdc,ii,1:3)=upl(kdc,ii,1:3)/DLL
  upl(kdc,ii,4)=-sum( upl(kdc,ii,1:3) * vt(kdc,npo(kdc,ii,1),1:3) )
 
if(ipo(kdc,ii).ge.4) then
do jj=3,ipo(kdc,ii)-1
    dd1=sum(upl(kdc,ii,1:3)*vt(kdc,npo(kdc,ii,jj),:))+ upl(kdc,ii,4)
    if(abs(dd1).gt.1.0d-9) then
        write(6,*) 'ia=',ii,'npo=',npo(kdc,ii,jj),dd1
        ierror=1
        if(loop.eq.100) then
        write(6,*) 'unv loop=',loop
        stop
        endif
    endif
enddo
if(ierror.eq.1) then
 loop=loop+1
  DL1(:)=vt(kdc,npo(kdc,ii,3),:)-vt(kdc,npo(kdc,ii,2),:)
  DL2(:)=vt(kdc,npo(kdc,ii,1),:)-vt(kdc,npo(kdc,ii,2),:)
  ierror=0
  goto 238
endif
endif

tmpV=ZERO
do ib=2,ipo(kdc,ii)-1
    tmpV=crossV(vt(kdc,npo(kdc,ii,ib),:)-vt(kdc,npo(kdc,ii,1),:),vt(kdc,npo(kdc,ii,ib+1),:)-vt(kdc,npo(kdc,ii,1),:))+tmpV
enddo
    dll=dot_product(tmpV,upl(kdc,ii,1:3))
if(dll.lt.ZERO) then
    upl(kdc,ii,1:4)=upl(kdc,ii,1:4)*(-ONE)
endif
enddo

endsubroutine UNV

subroutine set_dt3d
    use variables, only: cl_set, dt, g, dtmax, ZERO,                          &
                         is, ie, js, je, ks, ke, inn2d
    use arrays, only: x, dx, nff, mn, fb, in2d, dx, in2d, sp
    implicit none
    integer :: nn2d, i, j, k
    real*8 :: cmax, h_in, depth
    !
    ! Do not need to calc dt if not equal to zero
    if (dt.ne.ZERO) return
    ! Decide dt from courant
    ! Use max. initial mean water depth
    cmax = ZERO
    do nn2D = 1, inn2d
        i = in2d(0,nn2d) ; j = in2d(1,nn2d)
        ! Get the free surface level
        h_in = -1d4
        do k = ke-1,ks,-1
            if (nff(mn(i,k,j))%f.eq.2) then
                h_in = sp(mn(i,k,j),2)
                exit
            endif
        enddo
        if (h_in.eq.-1d4) cycle
        ! Get the initial depth
        do k = ks,ke-1
            if (nff(mn(i,k,j))%f.eq.-1) cycle
            !Subtracts the ground level from the initial free surface level
            depth = h_in - x(2,k+1) + fb(mn(i,k,j)) * dx(2,k)
            exit
        enddo
        if (depth.le.ZERO) cycle
        if (is.eq.ie-1) then
            cmax = max( cmax, sqrt( g * depth ) / dx(1,j) )
        elseif (js.eq.je-1) then
            cmax = max( cmax, sqrt( g * depth ) / dx(0,i) )
        else
            cmax = max( cmax, sqrt( g * depth ) / min( dx(0,i),dx(1,j) ) )
        endif
    enddo
    dt = min( cl_set / cmax, dtmax )   
!
endsubroutine set_dt3d
    
subroutine check_timestep(cl1,cl2,cl3,cl4,cl5,ivstep)
    use variables, only: svmax, ivmax, icong, dtmax, dtm, dt,                 &
                         courant, ZERO, DT_lim
    implicit none
    integer,intent(in) :: ivstep 
    real*8,intent(in)  :: cl1, cl2, cl3, cl4, cl5
    real*8 :: vmaxx

! courant < cl1 即加速
! cl1 < average_100(courant) < cl2　加速 
! cl2 < average_100(courant) < cl3  減速
! cl2 < average_5(courant) < cl4 減速 
! cl5 < courant 即減速
! cl2 と cl3　の間を目指す。

    if (courant.gt.cl5) then
      vmaxx = courant
      do !Iterate to find new dt
        dt = dt / (1.1d0)
        vmaxx = vmaxx / 1.1d0
        if (vmaxx.lt.cl5) exit
      enddo
      dtm = dt
      write(6,*) 'courant=',courant,' > =',cl5,vmaxx,'dt=',dt,vmaxx
      if (dt.lt.DT_lim) then
        icong = 1000
        return
      endif
      icong = 0
      ivmax = 0
      svmax = ZERO
    elseif (courant.lt.cl1.and.dt.lt.dtmax) then
      if (courant*dtm/dt.lt.cl3) then
        dt = min(dtm * 1.1d0,dtmax) ; dtm = dt
        write(6,*) 'ave-courant kasoku=',courant,'dt=',dt
      endif  
      ivmax = 0
      svmax = ZERO
    endif
    
    if (ivmax.lt.ivstep) then
      svmax = svmax + abs(courant)
      ivmax = ivmax + 1
    endif
    
    if (ivmax.gt.ivstep/20.and.svmax/float(ivmax).ge.cl4) then
      dt = dtm / 1.1d0 ; dtm = dt
      write(6,*) 'ave-courant gensoku=',svmax/float(ivmax),'dt=',dt
      ivmax = 0
      svmax = ZERO
    elseif (ivmax.eq.ivstep) then
      if (svmax/float(ivmax).ge.cl3) then
        vmaxx = svmax/float(ivmax)
        do !Iterate to find new dt
            dt = dtm / 1.1d0
            vmaxx = vmaxx / 1.1d0
            if (vmaxx.lt.cl3) exit
        enddo
        dtm = dt
        write(6,*) 'ave-courant gensoku=',svmax/float(ivmax),'dt=',dt
!      elseif(svmax/float(ivmax).lt.0.01d0) then
!        dt=min(dt*1.1d0*1.1d0,dtmax)
!        write(6,*) 'ave-courant kasoku=',svmax/float(ivmax),'dt=',dt
      elseif(svmax/float(ivmax).lt.cl2) then
        if (courant*dtm/dt.lt.cl3) then
            vmaxx = svmax/float(ivmax)
            do !Iterate to find new dt
                dt = dtm * 1.1d0
                vmaxx = vmaxx * 1.1d0
                if (vmaxx.gt.cl2) exit
            enddo
            dt  = min(dt,dtmax) ; dtm = dt
            write(6,*) 'ave-courant kasoku=',svmax/float(ivmax),'dt=',dt
        endif
      else
        write(6,*) 'ave-courant keep=',svmax/float(ivmax),'dt=',dt
      endif
      svmax = ZERO
      ivmax = 0
    endif
endsubroutine check_timestep
    
subroutine hyoko2(kdc,cc,zh,zmax)
    use arrays, only: vt, ivt, ipl, ipo, npo, upl
    use variables, only: ZERO, ONE
    implicit none
    doubleprecision,intent(in)::cc(0:1),zmax
    doubleprecision::zh,cc_pc(1:2),aa(1:2,1:2),b_c(1:3),det,xx(1:2) !,zv(0:2),xmin,xmax,ymin,ymax,yy(1:2)
    doubleprecision::bb(1:2,1:2),a_b(1:2),dd(1:2),pc(0:1),pxmin,pxmax,pymin,pymax  !,ff1(1:2),ff2(1:2),ee(1:2,1:2)
!    doubleprecision,dimension(0:1,0:1)::abi,ab0
    integer::ib,icc,ib1,kdc,idet,ii
    doubleprecision,parameter::eps=1d-9
    !if required data is outside polygon return
    pxmin=minval(vt(kdc,1:ivt(kdc),1));pxmax=maxval(vt(kdc,1:ivt(kdc),1))
    pymin=minval(vt(kdc,1:ivt(kdc),2));pymax=maxval(vt(kdc,1:ivt(kdc),2))
    if (pxmin.gt.cc(0).or.pxmax.lt.cc(0).or.pymin.gt.cc(1).or.pymax.lt.cc(1)) return
    zh=10000.0d0
!!!!omp parallel do private(pxmin,pxmax,pymax,pymin,cc_pc,pc,ib,ib1,b_c,a_b,aa,det,dd,icc) schedule(dynamic) 
do ii=1,ipl(kdc)
    if(ipo(kdc,ii).eq.0) cycle
  !  if(upl(kdc,ii,2).eq.ZERO) cycle ! 面の法線に鉛直成分がない場合
    
    pxmin=minval(vt(kdc,npo(kdc,ii,1:ipo(kdc,ii)),1));if(cc(0).lt.pxmin) cycle
    pxmax=maxval(vt(kdc,npo(kdc,ii,1:ipo(kdc,ii)),1));if(cc(0).gt.pxmax) cycle
    pymin=minval(vt(kdc,npo(kdc,ii,1:ipo(kdc,ii)),2));if(cc(1).lt.pymin) cycle
    pymax=maxval(vt(kdc,npo(kdc,ii,1:ipo(kdc,ii)),2));if(cc(1).gt.pymax) cycle

    pc(0)=sum(vt(kdc,npo(kdc,ii,1:ipo(kdc,ii)),1))/dfloat(ipo(kdc,ii))
    pc(1)=sum(vt(kdc,npo(kdc,ii,1:ipo(kdc,ii)),2))/dfloat(ipo(kdc,ii))    

    icc=0;idet=0
    cc_pc=pc-cc
    do ib=1,ipo(kdc,ii)
        ib1=ib+1
    if(ib.eq.ipo(kdc,ii)) ib1=1
    b_c(1:2)=vt(kdc,npo(kdc,ii,ib1),1:2)-vt(kdc,npo(kdc,ii,ib),1:2)
    a_b(1:2)=vt(kdc,npo(kdc,ii,ib),1:2)-cc
    aa(1:2,1)=cc_pc(1:2)
    aa(1:2,2)=(-ONE)*b_c(1:2)
    det=aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)
    if(det.eq.ZERO) then
        idet=idet+1
        cycle
    endif
    
    bb(1,1)=aa(2,2);bb(2,2)=aa(1,1);bb(2,1)=aa(2,1)*(-ONE);bb(1,2)=aa(1,2)*(-ONE)
    !ee=MATMUL(bb,aa)/det
    dd=MATMUL(bb,a_b)/det
    xx=cc+cc_pc*dd(1)
  !  yy=cc+a_b(1:2)+b_c(1:2)*dd(2)
    !ff1(1:2)=cc_pc*dd(1)
    !ff2(1:2)=a_b(1:2)+b_c(1:2)*dd(2)
     if(sqrt(sum((xx-vt(kdc,npo(kdc,ii,ib),1:2))**2)).lt.eps) then !頂点
        if(abs(dd(1)).lt.eps) then
            goto 234
        else
            cycle 
        endif
     else if (sqrt(sum((xx-vt(kdc,npo(kdc,ii,ib1),1:2))**2)).lt.eps) then
        if(abs(dd(1)).lt.eps) then
            goto 234
        else if(dd(1).ge.ZERO) then
           icc=icc+1
           cycle
        endif
     endif
    !else
        if(dd(2).le.ONE.and.dd(2).ge.ZERO.and.abs(dd(1)).lt.eps) then !オンライン
            goto 234
        else if(dd(2).le.ONE.and.dd(2).ge.ZERO.and.dd(1).ge.ZERO) then
            icc=icc+1
        endif
    !endif
    enddo

    if(icc.eq.2.or.icc.eq.4) then
        cycle
    else if(icc.eq.1.or.icc.eq.3.or.icc.eq.5) then
        goto 234

    else if(icc.eq.0) then
        cc_pc=pc-cc
        if(cc_pc(1)**2+cc_pc(2)**2.eq.ZERO) then
         continue
        goto  234
        else
            if(upl(kdc,ii,3).ne.ZERO) then
         write(6,*) 'upl=',upl(kdc,ii,:)
        write(6,*) 'icc=',icc,'kdc=',kdc,'ii (ipl no.)=',ii,'idet=',idet
        write(6,*) 'cc=',cc(0),cc(1),'cc_pc=',cc_pc(1),cc_pc(2)
            endif
        endif
        !stop
    else
        continue
    endif
enddo   
    
    return
234 if (zmax.eq.-10000d0) then
        zh=10000d0
    else
        zh=(-upl(kdc,ii,1)*cc(0)-upl(kdc,ii,2)*cc(1)-upl(kdc,ii,4))/upl(kdc,ii,3)
    endif
    return
!!!!!omp end parallel do    
    endsubroutine

function first_upwind(uc,vc,vp,Vm)
    use variables, only: ZERO
    use type_list, only: velocity
    implicit none    
    doubleprecision::first_upwind,uc
    type(velocity)::Vp,vc,Vm
    first_upwind = MAX(UC,ZERO)*(vc%v-Vm%v)/(vc%x-Vm%x) &
                  -MAX(-UC,ZERO)*(Vp%v-vc%v)/(Vp%x-vc%x)
endfunction

subroutine naigai2(cc,ipp,pp,icc1,surf)
use interface_list,only:dista
use variables, only: ZERO, ONE
implicit none
integer,intent(in)::ipp
doubleprecision,intent(in)::cc(0:2),pp(1:ipp,0:2),surf(0:3)
doubleprecision::pc(0:2),cc_pc(0:2),b_c(0:2),a_b(0:2),aa(1:2,1:2),det,bb(1:2,1:2),dd(1:2),d3,de(0:2)
doubleprecision::r1(0:2,0:2),r2(0:2,0:2),rr(0:2,0:2) !,across(0:2)
integer::icc,ib,ib1
integer,intent(out)::icc1

do ib=1,ipp
    if(dista(cc-pp(ib,:)).lt.1.0d-9) then
        continue
        icc1=2
        return
    endif
enddo

!ポリゴンの中心位置点pcを計算
pc(0)=sum(pp(1:ipp,0))/dfloat(ipp)
pc(1)=sum(pp(1:ipp,1))/dfloat(ipp)
pc(2)=sum(pp(1:ipp,2))/dfloat(ipp)  
d3=sum(surf(0:2)*pc(0:2))+surf(3)
continue

if(surf(0)**2+surf(1)**2.eq.ZERO) then ! ポリゴンが水平のとき
    rr=ZERO
    rr(0,0)=ONE;rr(1,1)=ONE;rr(2,2)=ONE
else 
!  ポリゴンをz軸周りに回転してx成分を０にする。 ( r1(m.n)とすると、mは行列の行（横列）　nは行列の列（縦列）)
r1(0,0)=surf(1)/sqrt(surf(0)**2+surf(1)**2)
r1(0,1)=(-ONE)*surf(0)/sqrt(surf(0)**2+surf(1)**2)
r1(0,2)=ZERO
r1(1,0)=surf(0)/sqrt(surf(0)**2+surf(1)**2)
r1(1,1)=r1(0,0) !surf(1)/sqrt(surf(0)**2+surf(1)**2)
r1(1,2)=ZERO
r1(2,0)=ZERO
r1(2,1)=ZERO
r1(2,2)=ONE
de=MATMUL(r1,surf(0:2))
!  回転したポリゴンをx軸周りに回転して、y成分も0にする。
r2(0,0)=ONE
r2(0,1)=ZERO
r2(0,2)=ZERO
r2(1,0)=ZERO
r2(1,1)=de(2)/sqrt(de(1)**2+de(2)**2)
r2(1,2)=(-ONE)*de(1)/sqrt(de(1)**2+de(2)**2)
r2(2,0)=ZERO
r2(2,1)=de(1)/sqrt(de(1)**2+de(2)**2)
r2(2,2)=r2(1,1)  !de(2)/sqrt(de(1)**2+de(2)**2)
!de=MATMUL(r2,de(0:2))
rr= MATMUL(r2,r1)  !ポリゴンを回転して、ポリゴン上でｚの値を一定にする回転行列。
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
        aa(1:2,2)=(-ONE)*b_c(0:1)
        det=aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)
        if(abs(det).lt.1.0d-10) cycle
        bb(1,1)=aa(2,2);bb(2,2)=aa(1,1);bb(2,1)=aa(2,1)*(-ONE);bb(1,2)=aa(1,2)*(-ONE)
        dd=MATMUL(bb,a_b(0:1))/det  

    if(dd(2).lt.ONE+1.0d-9.and.dd(2).gt.-1.0d-9.and.abs(dd(1)).lt.1.0d-9) then !オンライン
            icc1=2
            return
    else if(abs(dd(2)).lt.1.0d-9.and.dd(1).ge.ZERO) then
            icc=icc+1
    else if(abs(dd(2)-ONE).lt.1.0d-9.and.dd(1).ge.ZERO) then
            continue
    else if(dd(2).gt.ZERO.and.dd(2).lt.ONE.and.dd(1).ge.ZERO) then
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
    !    write(6,*) 'icc=',icc
        continue
    else
        continue
    endif
endsubroutine  
    
doubleprecision function renzoku_n(nn)
    use interface_list, only: renzoku_fv
    use arrays, only: in, kn, jn, fb, a, un, dx, mn
    use variables, only: ZERO
    implicit none
    integer::nn,i,k,j
    i=in(nn);k=kn(nn);j=jn(nn)
    if(fb(nn).ne.ZERO) then 
#ifndef ENTOU      
      renzoku_n=renzoku_fv( &
                            a(0,mn(I+1,k,j)),a(0,nn),a(1,mn(I,k,j+1)),a(1,nn),a(2,mn(I,k+1,j)),a(2,nn) &
                           ,un(0,mn(I+1,k,j)),un(0,nn),un(1,mn(I,k,j+1)),un(1,nn),un(2,mn(I,k+1,j)),un(2,nn) &
                           ,dx(0,i),dx(1,j),dx(2,k),fb(nn) )
#else                       
      renzoku_n=renzoku_fvr( &
                             a(0,mn(I+1,k,j)),a(0,nn),a(1,mn(I,k,j+1)),a(1,nn),a(2,mn(I,k+1,j)),a(2,nn) &
                            ,un(0,mn(I+1,k,j)),un(0,nn),un(1,mn(I,k,j+1)),un(1,nn),un(2,mn(I,k+1,j)),un(2,nn) &
                            ,dx(0,i),x(1,j+1),dx(1,j),dx(2,k),fb(nn) )
#endif                       
    else
      write(6,*) 'fb=',fb(nn),'nn=',nn,'error'
      renzoku_n=-10000.d0
    endif
endfunction
    
doubleprecision function renzoku_p(nn)
    use interface_list,only:renzoku_fv
    use arrays,only:in,kn,jn,fb,a,u,dx,mn
    use variables, only: ONE, ZERO
    implicit none
    integer::nn,i,k,j
    i=in(nn);k=kn(nn);j=jn(nn)
    if(fb(nn).ne.ZERO) then    
      renzoku_p=renzoku_fv( &
                            a(0,mn(I+1,k,j)),a(0,nn),a(1,mn(I,k,j+1)),a(1,nn),a(2,mn(I,k+1,j)),a(2,nn) &
                           ,u(0,mn(I+1,k,j)),u(0,nn),u(1,mn(I,k,j+1)),u(1,nn),u(2,mn(I,k+1,j)),u(2,nn) &
                           ,dx(0,i),dx(1,j),dx(2,k),fb(nn) )                      
    else
      write(6,*) 'fb=',fb(nn),'nn=',nn,'error'
      renzoku_p=-10000.d0
    endif
endfunction
    
doubleprecision function cal_drr(n1)
    use interface_list, only: renzoku_fv
    use arrays, only: a, fb, un, f, in, jn, kn, mn, dx
    use variables,only: dt, ONE, TWO
    implicit none
    integer,intent(in)::n1
    integer::i,j,k
        i=in(n1);j=jn(n1);k=kn(n1)
#ifdef DRI
#ifdef ENTOU
        cal_drr=renzoku_fv3r(ax2(mn(i+1,k,j)),ax2(n1),ay2(mn(i,k,j+1)),ay2(n1),az2(mn(i,k+1,j)),az2(n1),&
        un(0,mn(I+1,k,j)),UN(n1),un(1,mn(I,k,j+1)),un(1,n1),un(2,mn(I,k+1,j)),un(2,n1),&
        dx(0,i),y(j+1),dx(1,j),dx(2,k),fb(n1),fbd(n1),dt)&
#else
        cal_drr=renzoku_fv3(ax2(mn(i+1,k,j)),ax2(n1),ay2(mn(i,k,j+1)),ay2(n1),az2(mn(i,k+1,j)),az2(n1),&
        un(0,mn(I+1,k,j)),UN(n1),un(1,mn(I,k,j+1)),un(1,n1),un(2,mn(I,k+1,j)),un(2,n1),&
        dx(0,i),dx(1,j),dx(2,k),fb(n1),fbd(n1),dt)  &
#endif ENTOU     
     !              +(ONE-f(n1))/dt*2.0d0/(fb(n1)+fbd(n1))*fb(n1)
     !              +max(-1.0d-5,min(ONE-f(n1),1.0d-5))/dt*2.0d0/(fb(n1)+fbd(n1))*fb(n1)
           +max( -0.01d0,min(0.01d0,(ONE-f(n1))*fb(n1)) ) /dt*TWO/(fb(n1)+fbd(n1))
               
#else     !濁水計算はこれ
        cal_drr=renzoku_fv(a(0,mn(i+1,k,j)),a(0,n1),a(1,mn(i,k,j+1)),a(1,n1),a(2,mn(i,k+1,j)),a(2,n1),&
        un(0,mn(I+1,k,j)),un(0,n1),un(1,mn(I,k,j+1)),un(1,n1),un(2,mn(I,k+1,j)),un(2,n1),&
        dx(0,i),dx(1,j),dx(2,k),fb(n1))&
                   +(ONE-f(n1))*0.001d0/dt*fb(n1)*TWO/fb(n1)
#endif DRI
endfunction cal_drr
!------------------------------------------------------------------------------------            
doubleprecision function cal_drr2(n1)
    use interface_list, only: cal_drr 
    use arrays, only: a ,fb, dx, un, in, jn, kn, dx, mn
    use variables, only: dt
    implicit none
#ifdef DRI
    doubleprecision::renzoku_fv3
#else
    doubleprecision::renzoku_fv
#endif
    integer,intent(in)::n1
    integer::i,j,k
    i=in(n1);j=jn(n1);k=kn(n1)
#ifdef DRI
#ifdef ENTOU
    cal_drr2=renzoku_fv3r(ax2(mn(i+1,k,j)),ax2(n1),ay2(mn(i,k,j+1)),ay2(n1),az2(mn(i,k+1,j)),az2(n1),&
    un(0,mn(I+1,k,j)),un(0,n1),un(1,mn(I,k,j+1)),un(1,n1),un(2,mn(I,k+1,j)),un(2,n1),&
    dx(0,i),y(j+1),dx(1,j),dx(2,k),fb(n1),fbd(n1),dt)
#else
    cal_drr2=renzoku_fv3(ax2(mn(i+1,k,j)),ax2(n1),ay2(mn(i,k,j+1)),ay2(n1),az2(mn(i,k+1,j)),az2(n1),&
    un(0,mn(I+1,k,j)),un(0,n1),un(1,1,mn(I,k,j+1)),un(1,n1),un(2,mn(I,k+1,j)),un(2,n1),&
    dx(0,i),dx(1,j),dx(2,k),fb(n1),fbd(n1),dt)  
#endif ENTOU     
               
#else   !濁水計算はこれ
    cal_drr2=renzoku_fv(a(0,mn(i+1,k,j)),a(0,n1),a(1,mn(i,k,j+1)),a(1,n1),a(2,mn(i,k+1,j)),a(2,n1),&
    un(0,mn(I+1,k,j)),un(0,n1),un(1,mn(I,k,j+1)),un(1,n1),un(2,mn(I,k+1,j)),un(2,n1),&
    dx(0,i),dx(1,j),dx(2,k),fb(n1))
#endif DRI
end function cal_drr2           
