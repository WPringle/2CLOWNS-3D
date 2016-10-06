      subroutine ss_main_loop(flag,nn,mm0_p,mmd_p,mmd2_p,ncpl2,ncpl5,ncount0,pp_p,info_p,ipp_p ,axyz &
#ifdef LARGEDATA
                              ,nn_ia_n,nn_ia &
#endif
                              )
      use variables,only:n,inns,inne,ndri,npln,mm0,ncpl00max,ncpl01max,ncpl02max,pi,minfopath, HALF
      use arrays,only:nff,x,dx,cpl,vtmx,vtmn,in,jn,kn,mn,rhod,a0,vt,ipo,ipl,npo,xxmin,xxmax,upl,pp_0,ncpl0,info_0 &
                      ,mmd,mmd2,idri_kdc,idri_fc_n,idri_fc_all_n,ivt,men,path0,ncplc,drift,objno,upl0,drino2
      use  interface_list,only:crossl,crossv,dista,nodeofLineandPlain
      implicit none

      interface

        function side_num(ic1,ic2)
        implicit none
        integer::side_num
        integer::ic1,ic2
        end function

        function cal_near_point(A,B,face1,length)
        implicit none
        real*8::cal_near_point(0:2)
        real*8,intent(in)::A(0:2),B(0:2),face1(0:3)
        real*8,intent(in)::length
        end function

      end interface
!#include "interface.inc"

      ! main_loop_def
      integer::ipathu_bk,flag,nn,ncount0,icc,ic1,ip,ip1,icpl,ica,icb,icd,loop,loop2,ii1,ip0=8,icc1,icyc,icc2,icc3,icc4,icc00,kdc,objno_nn
      integer::iia,ipath,ipatho,ipathu,spath,incpl2,loop1,ireturn,ib2,ib3,ib3p,nnp,nnp1,i1,ipoia !,ic2p,vt_grnp,nn_ia_n
      integer::ic2,ic22,ic3,ic4,i_small,i_large,idd4(1:ncpl00max),idd5(10),ip2,ip3,ip4,ia1,ips !,icpl_bk,vt_grn
      integer::ncpl20(-2:-1),i_in,mmc,ia,ib,ib1,inear,i,ii,ic,ncpl2s(1:ncpl00max),icme1(5,5),ipcc(20),ipccn(20),ML(20),inear_c !,kk
      integer::mm0_p,ipp_p,mmd_p(1:ndri),mmd2_p(1:ndri)  !,mm0mxp,mmd2mxp
      doubleprecision,dimension(npln,0:3)::cpl2
      doubleprecision::distmax,distmin,dd1,dd2,dd3,dpcc(20),tt1,tt2,tt3,tt4,VV,vv1,tmpp(1:6),length=1.0d-7 !,dist1,ddd,tmp1,tmpQ,tmpR
      integer,dimension(0:ncpl00max,-2:ncpl02max)::ncpl2,ncpl3,ncpl4
      integer,dimension(0:100,-3:ncpl02max)::ncpl5
      doubleprecision,dimension(1:ncpl02max,0:2)::vtia
      doubleprecision,dimension(0:3)::surf2,face1   !,cpl_tmp
      doubleprecision,dimension(0:2)::tmpV,tmpV1,tmpV2,tmpV3 !,tmpBA,tmpV4,tmpV5
      integer,dimension(1:minfopath)::icv0,icv,icv2,icv22,icv3,icv4,icv41,sinfo,icv1,icvm,icvn,icvnt,icv4a  !,ixd,nnx,imen,ixe,ipoint,iko
      doubleprecision,allocatable,dimension(:,:)::pp !,ppp
      doubleprecision,dimension(1:ncpl00max)::dd4,s
      integer,allocatable,dimension(:)::info,info1,info2 !,pl_gr,vt_gr,imax
      integer,allocatable,dimension(:,:)::path,info3
      integer,dimension(0:minfopath,-2:4)::path2
      integer,dimension(0:minfopath,-2:3)::paths,pathu,patho,path3
      doubleprecision,dimension(0:2)::x1max,x1min  !,OA,OB,OQ,OP,AP,AB,PQ,OXg21,g22
      integer,dimension(1:50,1:4)::icdr   !icv5,icv6,
      integer,dimension(1:50,1:3,1:4)::icme
      doubleprecision,allocatable,dimension(:)::leng !,vv0
!      doubleprecision,allocatable,dimension(:,:,:,:)::dummyr4
!      doubleprecision,allocatable,dimension(:,:,:)::dummyr3
!      integer,allocatable,dimension(:,:,:)::dummyi3
      doubleprecision,allocatable,dimension(:,:)::dummyr2
!      integer,allocatable,dimension(:,:)::dummyi2
      integer,allocatable,dimension(:)::dummyi1 !,icv_info
      doubleprecision::axyz(0:6)
      doubleprecision::pp_p(0:ipp_p,0:2),pp_min(0:2),pp_max(0:2)
      doubleprecision::clvt2plypl(1:8) ! セル頂点からポリゴン面までの垂線長さ
      integer::info_p(1:2,0:ipp_p),ncplc00
#ifdef LARGEDATA
      integer::ia_1
      integer::nn_ia(1:nn_ia_n)
#endif

ipathu_bk=0
4538 continue  
ncpl5=0

inear_c=0

#ifdef BUGWIN
if(nn.eq.124215) then
    continue
endif
#endif
!if(flag.eq.0) OBJNO(nn)=0       !固定物が存在する場合，固定物の通し番号を初期化.
ncount0=ncount0+1
!漂流物以外の物体セルは計算しない．次のセルへ
if(nff(nn)%f.eq.-1.and.nff(nn)%b.ne.10) return
if(nff(nn)%b.eq.20) return  !special

!if(in(nn).lt.53.or.in(nn).gt.141) cycle
!if(jn(nn).lt.155.or.jn(nn).gt.217) cycle
!if(in(nn).lt.53.or.in(nn).gt.141) cycle
!if(jn(nn).lt.155.or.jn(nn).gt.217) cycle

#ifdef DDDDDD
#ifndef DRI
!一つ下が気体or流体セルでかつ空隙率が１の場合計算しない．（地形作成で有効））
if(nf(in(nn),kn(nn)-1,jn(nn)).ge.0.and.nf(in(nn),kn(nn)-2,jn(nn)).ge.0.and.flag.eq.0) then
if(fb0(mn(in(nn),kn(nn)-1,jn(nn))).eq.1.0d0.and.fb0(mn(in(nn),kn(nn)-2,jn(nn))).eq.1.0d0) then
!#ifndef BUGWIN    
cycle
!#endif
endif
endif
#endif
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!　セル頂点領域が対象物体の頂点領域に含まれるか調べる．!!!!!!
!!!!!!!!!!!!!!!!!  ここから   !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(allocated(objno)) then
if(objno(nn).ne.0) then ! 地形があるとき
    objno_nn=objno(nn)
    ip=ncpl0(objno_nn,0,1)  !セル頂点数を登録
!    write(6,*) ncpl0(objno_nn,0,0),ncpl0(objno_nn,0,1)
    do i=0,2
        pp_min(i)=minval(pp_0(objno_nn,1:ip,i))
        pp_max(i)=maxval(pp_0(objno_nn,1:ip,i))
    enddo
    goto 434
endif
endif

objno_nn=0
    ip=8  !セル頂点数を登録
    pp_min(0:2)=(/x(0,in(nn)),x(1,jn(nn)),x(2,kn(nn))/)     ! xmin
    pp_max(0:2)=(/x(0,in(nn)+1),x(1,jn(nn)+1),x(2,kn(nn)+1)/)   ! xmax    

434    ip0=ip
 





!セル頂点の一部を設定．

icc=0;ic1=0
kdc_loop_pre: do kdc=1,ndri    !kdc loop
if(rhod(kdc).gt.0.0d0.and.flag.eq.0) cycle   !地形作成で，漂流物は計算しない．
if(rhod(kdc).lt.0.0d0.and.flag.eq.1.and.n.gt.0) cycle   !漂流物計算で，不動物は計算しない．

ic1=ic1+1  !調査した物体数をカウント
if(pp_min(0).gt.vtmx(kdc,0).or.pp_min(1).gt.vtmx(kdc,1).or.pp_min(2).gt.vtmx(kdc,2).or.&
   pp_max(0).lt.vtmn(kdc,0).or.pp_max(1).lt.vtmn(kdc,1).or.pp_max(2).lt.vtmn(kdc,2)) then
!物体頂点の内側にセル頂点が含まれない場合，カウント   
icc=icc+1
endif
enddo kdc_loop_pre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(icc.eq.ic1) return                                 !セル頂点領域が対象物体の頂点領域に含まれない場合終了　-> 次のセルnn+1へ!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  ここまで   !!!!!!!!!!!!!!!!!!!!!!!!
!!!!　セルが対象物体の頂点群に含まれるか調べる．!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(pp(0:ncpl01max,0:2))


if(objno_nn.eq.0) then

!セル頂点を設定，
pp(1,0:2)=pp_min(0:2)  ! xmin
pp(2,0:2)=(/x(0,in(nn)+1),x(1,jn(nn)),x(2,kn(nn))/)   !=pp(1,0:2)+(/dx(0,in(nn)),0.0d0,0.0d0/)
pp(3,0:2)=(/x(0,in(nn)+1),x(1,jn(nn)+1),x(2,kn(nn))/)    !pp(1,0:2)+(/dx(0,in(nn)),dx(1,jn(nn)),0.0d0/)
pp(4,0:2)=(/x(0,in(nn)),x(1,jn(nn)+1),x(2,kn(nn))/)   !=pp(1,0:2)+(/0.0d0,dx(1,jn(nn)),0.0d0/)
pp(5,0:2)=(/x(0,in(nn)),x(1,jn(nn)),x(2,kn(nn)+1)/)  !=pp(1,0:2)+(/0.0d0,0.0d0,dx(2,kn(nn))/)
pp(6,0:2)=(/x(0,in(nn)+1),x(1,jn(nn)),x(2,kn(nn)+1)/)  !=pp(2,0:2)+(/0.0d0,0.0d0,dx(2,kn(nn))/)
pp(7,0:2)=pp_max(0:2)  ! xmax
pp(8,0:2)=(/x(0,in(nn)),x(1,jn(nn)+1),x(2,kn(nn)+1)/)  !=pp(4,0:2)+(/0.0d0,0.0d0,dx(2,kn(nn))/)


!各セル境界面を設定．
cpl2(1:6,:)=cpl(1:6,:) !各セルの境界面情報をcpl2に代入．
cpl2(1,3)=-1.0d0*sum(cpl(1,0:2)*pp(1,0:2))  
cpl2(2,3)=-1.0d0*sum(cpl(2,0:2)*pp(5,0:2))
cpl2(3,3)=-1.0d0*sum(cpl(3,0:2)*pp(1,0:2))
cpl2(4,3)=-1.0d0*sum(cpl(4,0:2)*pp(2,0:2))
cpl2(5,3)=-1.0d0*sum(cpl(5,0:2)*pp(4,0:2))
cpl2(6,3)=-1.0d0*sum(cpl(6,0:2)*pp(1,0:2))
   
ncpl2(0:6,-2:4)=ncplc(0:6,-2:4)  !各セルの境界面情報をncpl2に代入
icpl=6 !ncpl(0,0)
else
    
pp(1:ip,:)=pp_0(objno_nn,1:ip,:)
icpl=ncpl0(objno_nn,0,0)
ncpl2(0:icpl,-2:ncpl0(objno_nn,0,2))=ncpl0(objno_nn,0:icpl,-2:ncpl0(objno_nn,0,2))
do i=1,icpl
    if(ncpl2(i,-1).eq.0) then
        cpl2(i,:)=cpl(ncpl2(i,-2),:)
        cpl2(i,3)=-1.0d0*sum(cpl2(i,0:2)*pp_0(objno_nn,ncpl2(i,1),0:2))  
        continue
    else if(ncpl0(objno_nn,i,-1).eq.1) then
        cpl2(i,:)=upl0(ncpl2(i,-2)-6,:)
        continue
    else
        continue
    endif
enddo
continue

endif

ncplc00=icpl

allocate(info(0:ncpl01max))
ipatho=0
info=0
i_in=0
!iac=1

if(mod(ncount0,50000).eq.1.and.flag.eq.0) then
 !write(6,*) 'nn,mmd=',nn,mmd_p(1),'inne=',inne
#ifdef LARGEDATA
 write(6,*) 'nn,mmd=',nn,mmd(1)  !,'nn_nne=',nn_nne
 !write(6,*) 'nn,mmd=',nn,mmd(1),'nn_1,nn_nne=',nn_1,nn_nne
#endif
endif

!各漂流物について

kdc_loop:  do kdc=1,ndri

icc=0;ic1=0
if(rhod(kdc).gt.0.0d0.and.flag.eq.0) cycle
if(rhod(kdc).lt.0.0d0.and.flag.eq.1.and.n.gt.0) cycle


mmc=0   !セルnn内の物体kdcに関連した面の数を初期化する．

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef LARGEDATA
!ia_loop:do ia_1=1,nn_ia_n  !nn_ia(nn,0)   !!!  対象物体表面iaについて
!    ia=nn_ia(ia_1)
#else
ia_loop:do ia=1,ipl(kdc)   !!!  対象物体表面iaについて
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  


face1=upl(kdc,ia,:)        !!!  対象物体表面iaの面情報をface1に代入.
ipoia=ipo(kdc,ia)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!start start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!物体表面iaを含む領域の中にセル頂点が含まれるか　 !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(pp_min(0)-xxmax(kdc,ia,0).gt.1.0d-9.or.pp_min(1)-xxmax(kdc,ia,1).gt.1.0d-9.or.pp_min(2)-xxmax(kdc,ia,2).gt.1.0d-9.or.&
   pp_max(0)-xxmin(kdc,ia,0).lt.-1.0d-9.or.pp_max(1)-xxmin(kdc,ia,1).lt.-1.0d-9.or.pp_max(2)-xxmin(kdc,ia,2).lt.-1.0d-9) then
    
!物体表面iaを含む領域の中にセル頂点が含まれない場合，（物体面領域とセル領域が交わらない）

    !pp_minとpp_maxの平面との関係を計算   
    distmax=sum(face1(0:2)*pp_min(0:2))+face1(3)
    distmin=sum(face1(0:2)*pp_max(0:2))+face1(3)
    if(distmax.gt.0.0d0.or.distmin.gt.0.0d0) then
        !pp_minとpp_maxのどちらかが平面の外側にある場合，フラグを立てる．
        i_in=1
    endif
#ifdef DRI22
    if(pp_max(0)-xxmin(0).lt.-1.0d-9.and.pp_max(1)-xxmin(1).lt.-1.0d-9) then
        iac=ia

        exit
    endif
#endif
!次のセルに移る．
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cycle                                                               !!!!!!!!!!物体表面iaを含む領域の中にセル頂点が含まれない場合!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if(xxmax(kdc,ia,0).le.pp_max(0)-1.0d-9.and.xxmax(kdc,ia,1).le.pp_max(1)-1.0d-9.and.xxmax(kdc,ia,2).le.pp_max(2)-1.0d-9.and.&
   xxmin(kdc,ia,0).ge.pp_min(0)+1.0d-9.and.xxmin(kdc,ia,1).ge.pp_min(1)+1.0d-9.and.xxmin(kdc,ia,2).ge.pp_min(2)+1.0d-9) then
! セル内に全頂点が含まれる．
    continue
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!物体表面iaを含む領域の中にセル頂点が含まれるか？ !!!!!!!!
!!!!!!!!! end end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
#ifdef BUGWIN
if(nn.eq.124215.and.ia.eq.3) then
continue
endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!start start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!セル頂点群内に対象物体表面が含まれるか？ !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!<<セル頂点と対象物体表面との距離を調べる
!<<また，対象物体表面に含まれるセル頂点数を登録する．
ncpl2s(:)=0
icc=0
do ii=1,ip0                      
    clvt2plypl(ii)=sum(face1(0:2)*pp(ii,0:2))+face1(3)   !セル頂点と対象物体表面iaとの距離を計算
    if(abs(clvt2plypl(ii)).lt.1.0d-9) then
         icc=icc+1          !対象物体表面iaに含まれるセル頂点数iccをカウントして，
         icv0(icc)=ii       !icv0(1:icc)に登録
    endif
enddo
distmax=maxval(clvt2plypl(1:ip0)) ; if(abs(distmax).lt.1.0d-9) distmax=0.0d0    !距離dist1の最大値を調べる
distmin=minval(clvt2plypl(1:ip0)) ; if(abs(distmin).lt.1.0d-9) distmin=0.0d0     !距離dist1の最小値を調べる
!<<セル頂点群に対象物体表面iaが含まれない（対象物体表面iaがセルを切断しない）場合終了　→次の物体表面へ
if(distmax*distmin.gt.0.0d0) then
    if(distmin.gt.0.0d0) then
    !すべてセル頂点が，平面の外側にある場合，フラグを立てる．
    i_in=1
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cycle                    !セル頂点群に対象物体表面iaが含まれない（対象物体表面iaがセルを切断しない）場合終了　→次の物体表面ia+1へ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif

inear=0

!４っつのセル頂点が対象物体表面iaに含まれる場合，face1が何れかのセル境界面に一致するかを調べる．．
if(icc.eq.4) then                           !iccが4の場合， 
    do ib=1,6
        if(sum(abs(cpl2(ib,0:3)-face1(0:3))).lt.1.0d-8) then
            inear=ib   !face1がセル境界面のいずれかに一致する場合inear=0
            exit
        endif
        if(sum(abs(cpl2(ib,0:3)-(-1.0d0)*face1(0:3))).lt.1.0d-8) then
            inear=(-1)*ib  !face1がセル境界面のいずれかに一致する場合inear=0
        endif
    enddo
    if(inear.gt.0) then
        continue  
       ! cycle
    endif
 !   if(distmax.ge.0.0d0.or.distmin.ge.0.0d0) then
    if(distmin.ge.-1.0d-9) then    
    i_in=1
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    if(inear.lt.0) then     
        cycle
   ! else if(inear.gt.0) then
   !     inear_c=inear_c+1
   !     cycle
    endif
        !inear=0のとき，face1を物体表面と見なさない．→次の物体表面ia+1へ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!セル頂点群内に対象物体表面が含まれるか？ !!!!!!!!
!!!!!!!!! end end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!start start !!!!!!!!!!!!!!!!!!!!!!!!!!!
!物体表面によるセルの切断面を求める．!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!セル辺iiごとに物体表面との交点を探して，icv0(1:icc)に登録．
do ii=1,12
    !セル辺iiの両端の頂点，path0(ii,1)，path0(ii,2)と対象物体表面iaの距離を調べ，dd1,dd2に代入
    dd1=clvt2plypl(path0(ii,1))   !sum(face1(0:2)*pp(path0(ii,1),0:2))+face1(3)
    dd2=clvt2plypl(path0(ii,2))   !sum(face1(0:2)*pp(path0(ii,2),0:2))+face1(3)
    if(dd1*dd2.lt.0.0d0.and.abs(dd1).ge.1.0d-9.and.abs(dd2).ge.1.0d-9) then
        !セル辺iiの両端の頂点path0(ii,1)，path0(ii,2)間に対象物体表面との交点がある場合．
        icc=icc+1                    !物体表面とセル辺の交点数iccをカウント
        !交点の候補をtmpVとする．
        !tmpV(0:2)=( pp(path0(ii,1),:)*abs(dd2)+pp(path0(ii,2),:)*abs(dd1) )/(abs(dd1)+abs(dd2))
        tmpV= nodeofLineandPlain( pp(path0(ii,1),:), pp(path0(ii,2),:),face1)
        ip1=0
        !tmpVがすでに登録されている頂点と一致するか調べる．          
        do ib1=1,ip
            dd1=dista(tmpV-pp(ib1,0:2))
            if(dd1.lt.1.0d-9)  then ;ip1=ib1; exit ; endif
        enddo
        if(ip1.eq.0) then !一致しない場合               
            ip=ip+1                 !
            call realloc_pp_info
            pp(ip,0:2)=tmpV(0:2)    !tmpVを新たな頂点(通し番号ip)として登録
            ip1=ip
            endif
        icv0(icc)=ip1            !ipをicc番目の交点として登録
        info(ip1)=ii             !頂点ipが辺番号iiに含まれることを登録    
    endif
enddo

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
if(icc.le.2) cycle                                                     !セル辺を物体面の交点が2つ以下の場合終了　→次の物体表面ia+1 へ  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(inear.gt.0.and.icc.eq.4) then
    continue
endif

!見つけた交点群を含む領域を計算
icb=0
do ic=0,2
    x1max(ic)=maxval(pp(icv0(1:icc),ic))  ! 見つけたセル辺と面の交点の最大値
    x1min(ic)=minval(pp(icv0(1:icc),ic))  ! 見つけたセル辺と面の交点の最小値
enddo

if(x1min(0)-xxmax(kdc,ia,0).gt.1.0d-9.or.x1min(1)-xxmax(kdc,ia,1).gt.1.0d-9.or.x1min(2)-xxmax(kdc,ia,2).gt.1.0d-9.or.&
   x1max(0)-xxmin(kdc,ia,0).lt.-1.0d-9.or.x1max(1)-xxmin(kdc,ia,1).lt.-1.0d-9.or.x1max(2)-xxmin(kdc,ia,2).lt.-1.0d-9) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    continue
    cycle !見つけた交点群を含む領域が，対象物体面を含む領域の外側の場合．次の物体表面ia+1に移る．
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
endif
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!<見つけた交点pp(icv0(1:icc),0:2)を並び替える．
if(allocated(path)) deallocate(path)
allocate(path(0:minfopath,-2:4))
ipath=0
!セル境界面cpl(ib)に二つの交点が存在する場合，二つの交点がパスを形成する．
do ib=1,6      
    ica=0;icb=0
    do ic=1,icc
        !セル境界面cpl(ib)と交点pp(icv0(ic),0:2)の距離を計算
        dd1=sum(cpl2(ib,0:2)*pp(icv0(ic),0:2))+cpl2(ib,3)  
        if(ica.ne.0.and.icb.eq.0.and.abs(dd1).lt.1.0d-9) then
            !別の交点pp(icv0(ic),0:2)がセル境界面cpl(ib)に含まれる場合
            icb=icv0(ic)        !頂点番号をicbに登録
        endif
        if(ica.eq.0.and.abs(dd1).lt.1.0d-9) then
            !交点pp(icv0(ic),0:2)がセル境界面cpl(ib)に含まれる場合
            ica=icv0(ic)        !頂点番号をicaに登録
        endif
    enddo
    if(icb.ne.0) then   
        !セル境界面cpl(ib)に点pp(ica),pp(icb)が含まれる場合
        !path(ica->icb)が登録されていないか調べる．
        icd=0
        do ii=1,ipath
            if(path(ii,1).eq.ica.and.path(ii,2).eq.icb) icd=ii
            if(path(ii,2).eq.ica.and.path(ii,1).eq.icb) icd=ii
        enddo
        if(icd.eq.0) then
            !path(ica->icb)が登録されていない場合．
            ipath=ipath+1         !path(ica->icb)をipath番目のパスに登録．  
            path(ipath,1)=ica     ! 
            path(ipath,2)=icb     !
            path(ipath,3)=1       !未使用パスとして登録．
            path(ipath,-1)=ib     !パスが含まれるセル境界面番号
            path(ipath,0)=0
        else
            if(path(icd,-1).eq.1.and.ib.ge.3.and.ib.le.6) then
                path(icd,0)=98+ib
            else if(path(icd,-1).eq.2.and.ib.ge.3.and.ib.le.6) then
                path(icd,0)=102+ib
            else
                continue
            endif
        endif
    endif
enddo ! ib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!見つけた交点パスがつながるように並べ替える． start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
icv(1)=path(1,1)     !並べ替えの始点を　path(1,1) とする．　　
icd=1;loop=0
237 do ii=1,ipath
    if(path(ii,3).eq.0) cycle          !使用済みパスならば次へ
    if(path(ii,1).eq.icv(icd)) then    
    !ii番目のパスの一端が現在の交点と一致する場合，
    !ii番目のパスの反対側の頂点を次の交点とする．パスiiを使用済みパスとする．
        icd=icd+1
        icv(icd)=path(ii,2);path(ii,3)=0
    else if(path(ii,2).eq.icv(icd)) then
    !ii番目のパスの一端が現在の交点と一致する場合，
    !ii番目のパスの反対側の頂点を次の交点とする．パスiiを使用済みパスとする．
        icd=icd+1
        icv(icd)=path(ii,1);path(ii,3)=0
    endif
enddo
!現在の交点が始点と一致しないと237に戻る．
loop=loop+1
if(loop.gt.100) then
    write(6,*) 'loop 237 nn=',nn,mm0,mmd(1),mmd2(1),'ijk=',in(nn),jn(nn),kn(nn),'ipath=',ipath
    goto 932
endif
if(icv(icd).ne.icv(1)) goto 237

!連続パスicv(1:icc)がface1の法線の向きと一致するように修正．
call surface_dir3(icc,icv(1:icc),face1,icv22(1:icc))
icv(1:icc)=icv22(1:icc)     
icc00=icc !対象物体表面のセルによる切断面の頂点数の並びをicv(1:icc00)とする．
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!見つけた交点パスがつながるように並べ替える． end  (切断面完成）
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!<切断面の頂点の並びをpath2に登録．
do ii=1,icc00
    ii1=ii+1;if(ii1.gt.icc00) ii1=1
    path2(ii,1)=icv(ii)
    path2(ii,2)=icv(ii1)
    path2(ii,3)=0
    path2(ii,4)=side_num(icv(ii),icv(ii1))
    !< 切断面の辺が属する平面情報(-1,0)を登録
    do ib=1,ipath
        if(path2(ii,1).eq.path(ib,1).and.path2(ii,2).eq.path(ib,2))  then
         path2(ii,-1:0)=path(ib,-1:0)
        else if(path2(ii,1).eq.path(ib,2).and.path2(ii,2).eq.path(ib,1)) then
         path2(ii,-1:0)=path(ib,-1:0)
        endif
    enddo
enddo
deallocate(path)

#ifdef BUGWIN
if(nn.eq.2344140.and.ia.eq.1232) then
continue
endif
#endif

do ic=1,ipoia  ! 漂流物の頂点をvtiaにコピー
    vtia(ic,:)=vt(kdc,npo(kdc,ia,ic),:)
enddo

icv3=0  ! 切断面頂点包含情報を初期化
do ii=1,icc00   ! 切断面頂点が漂流物面に含まれるか？　
    call naigai2(pp(icv(ii),:),ipoia,vtia(1:ipoia,0:2),icv3(ii),face1)
enddo


path3(1:ipoia,3)=0    !　漂流物面の辺情報を初期化．
!ikouten=0  !;ic3=0
icme=0;icdr(:,:)=0;;icv4=0;;icv41=0;icv4a=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!対象物体表面iaの頂点(総数ipoia)が切断面icv(1:icc00)に含まれるか？
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ic=1,ipoia
    do ii=1,icc00
    if(dista(vtia(ic,:)-pp(icv(ii),:)).lt.1.0d-9) then ! 物体頂点が切断面頂点と一致するとき
        icv4(ic)=icv(ii)
        icv3(ii)=2
        icv4a(ic)=2
        cycle
    endif
    enddo

 call naigai2(vtia(ic,:),icc00,pp(icv(1:icc00),0:2),icc1,face1)
     if(icc1.ne.0) then
            !物体面頂点が既存の頂点と一致しない場合
            ip1=0
            do ic1=1,ip
              dd3=dista(vtia(ic,:)-pp(ic1,0:2))
              if(dd3.lt.1.0d-9) then
              !交点が既存の頂点と一致する場合
              ip1=ic1;exit
              endif
            enddo
            if(ip1.eq.0) then
            ip=ip+1
            call realloc_pp_info
            pp(ip,:)=vtia(ic,:)
            ip1=ip
      
            endif
            if(info(ip1).eq.0)  info(ip1)=100      !切断面内に存在する物体表面の頂点としてinfo登録     
            icv4a(ic)=icc1
            icv4(ic)=ip1      !icv4にあらたな頂点番号ipを登録
     endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! end end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!物体表面によるセルの切断面を求める．!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef BUGWIN
if(ia.eq.1232.and.nn.eq.2344140) then
    continue 
endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!対象漂流物面の各辺path3と切断面の各辺path2(ii)の交点を求める．!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
path3(:,-1)=0
do ib=1,ipoia !(kdc,ia)
    ib1=ib+1;if(ib1.gt.ipoia) ib1=1
    
ips=ip
    
do ii=1,icc00
!path2(ii)を含むセル境界面cpl(path2(ii,-1),0:2)と漂流物面の辺との交点を調べる．
dd1=sum(cpl2(path2(ii,-1),0:2)*vtia(ib,:))+cpl2(path2(ii,-1),3)    !漂流物面の辺path3(ib)の両端Vt(...)とpath2(ii)を含むセル境界面の距離dd1,dd2
dd2=sum(cpl2(path2(ii,-1),0:2)*vtia(ib1,:))+cpl2(path2(ii,-1),3)

if(abs(dd1).lt.1.0d-9.and.abs(dd2).lt.1.0d-9) then 
    if(path3(ib,-1).ne.0) then
        continue
    else
        path3(ib,-1)=path2(ii,-1)
    endif
    continue
endif

if(dd1*dd2.gt.0.0d0.and.abs(dd1).ge.1.0d-9.and.abs(dd2).ge.1.0d-9) cycle
if(dd1*dd2.lt.0.0d0.and.abs(dd1).ge.1.0d-9.and.abs(dd2).ge.1.0d-9) then       !漂流物面の辺path3(ib)の両端Vt(...)が対象境界面から1.0d-9以上離れ,対象境界面を挟むとき．
    tmpV(0:2)=( vtia(ib,:)*abs(dd2)+vtia(ib1,:)*abs(dd1) )/(abs(dd1)+abs(dd2))   !交点をtmpV
    tmpV1=pp(path2(ii,2),:)-pp(path2(ii,1),:)
    tmpV2=tmpV(0:2)-pp(path2(ii,1),:)
    tmpV3=tmpV(0:2)-pp(path2(ii,2),:)
    if(dista(tmpv1).gt.dista(tmpv2).and.dot_product(tmpv1,tmpv2).gt.0.0d0.and.dista(tmpv2).ge.1.0d-9.and.dista(tmpv3).ge.1.0d-9) then
    !交点の候補とpp(path2(ii,1),:)の距離が，辺の長さdista(tmpv1)より小さく，pp(path2(ii,1),:)から交点候補への向きが，pp(path2(ii,1),:)から
    !pp(path2(ii,2),:)への向きと同じで，かつ，pp(path2(ii,1),:)，pp(path2(ii,2),:)にそんなに近くない．
    
            path2(ii,3)=path2(ii,3)+1    !切断面の辺path2(ii)の交点数path2(ii,3)を準備
            path3(ib,3)=path3(ib,3)+1    !漂流物面の辺path3(ib)の交点数path3(ib,3)を準備
            ip1=0
            do ic1=1,ip
              dd3=dista(tmpV-pp(ic1,0:2))
              if(dd3.lt.1.0d-9) then
              !交点が既存の頂点と一致する場合
              ip1=ic1      
              endif
            enddo
            if(ip1.eq.0) then
            !交点が既存の頂点と一致しない場合
            ip=ip+1
call realloc_pp_info         
            pp(ip,:)=tmpv        
            ip1=ip
            endif  
            ic4=0
            do ic3=1,path2(ii,3)-1   !すでに交点として登録されいるか確認．
                if(icme(ii,1,ic3).eq.ip1) ic4=1
            enddo
            if(ic4.eq.0) then
                icme(ii,1,path2(ii,3))=ip1  !path2(ii)のpath2(ii,3)番目の交点の頂点番号を登録
        !     write(6,*)  'icme(',ii,',1,',path2(ii,3),')=',ip1
                if(dd2.gt.0.0d0) then
                    icme(ii,3,path2(ii,3))=1 ! 中から外の場合
                else
                    icme(ii,3,path2(ii,3))=-1 ! 外から中の場合
                endif
            else
                path2(ii,3)=path2(ii,3)-1 ! すでに登録されている場合，番号を戻す．
            endif
            if(path2(ii,4).ne.0) then
                info(ip1)=path2(ii,4) ! path2(ii)がセル辺と一致するとき．
            endif
            
                icdr(ib,path3(ib,3))=ip1  !path3(ib)のpath3(ib,3)番目の交点の頂点番号をicdrに登録            
         !   write(6,*)' icdr(',ib,',',path3(ib,3),')=',ip1
    else if(dista(tmpv2).lt.1.0d-9) then
    !交点が，辺path2(ii)の一端path2(ii,1)に十分近いとき． 
            ic22=0
 !           if(path3(ib,3).ge.1) then
            do ic2=1,path3(ib,3)
                if(icdr(ib,ic2).eq.path2(ii,1)) then
                 ic22=1
                 exit
                endif 
            enddo
 !           endif
            if(ic22.eq.1) cycle
            path3(ib,3)=path3(ib,3)+1        !漂流物面の辺path3(ib)の交点数path3(ib,3)を登録
            icdr(ib,path3(ib,3))=path2(ii,1)     !path3(ib)のpath3(ib,3)番目の交点の頂点番号を登録
    !        icv3(ii)=2                           !icv(ii)を漂流物に含まれる切断面頂点として登録
    else if(dista(tmpv3).lt.1.0d-9) then
    !交点が，辺path2(ii)の一端path2(ii,2)に十分近いとき
                ic22=0
    !        if(path3(ib,3).ge.1) then
            do ic2=1,path3(ib,3)
                if(icdr(ib,ic2).eq.path2(ii,2)) then
                 ic22=1
                 exit
                endif 
            enddo
    !        endif    
             if(ic22.eq.1) cycle           
            path3(ib,3)=path3(ib,3)+1           !漂流物面の辺path3(ib)の交点数path3(ib,3)を登録
            icdr(ib,path3(ib,3))=path2(ii,2)    !path3(ib)のpath3(ib,3)番目の交点の頂点番号を登録
       !     ic3=ic3+1                           !漂流物に含まれる切断面頂点としてカウント
 !           if(ii.eq.icc00) then                !icv(ii+1)を漂流物に含まれる切断面頂点として登録
 !               icv3(1)=2          
 !           else
 !               icv3(ii+1)=2
 !           endif
        continue
    endif
else if(abs(dd1).lt.1.0d-9) then
!漂流物面の辺path3(ib)の一端Vt(kdc,npo(kdc,ia,ib),:)が対象境界面に含まれるとき．
    tmpV(0:2)=vtia(ib,:)
    tmpV1=pp(path2(ii,2),:)-pp(path2(ii,1),:)
    tmpV2=tmpV(0:2)-pp(path2(ii,1),:)
    tmpV3=tmpV(0:2)-pp(path2(ii,2),:)
#ifdef BUGWIN
!write(6,*) dista(tmpv1),dista(tmpv2),dot_product(tmpv1,tmpv2),dista(tmpv3)
!ddd1=dista(tmpv1);ddd2=dista(tmpv2);ddd3=dista(tmpv3);dddp=dot_product(tmpv1,tmpv2)
#endif
    if(dista(tmpv1).gt.dista(tmpv2).and.dot_product(tmpv1,tmpv2).gt.0.0d0.and.dista(tmpv2).ge.1.0d-9.and.dista(tmpv3).ge.1.0d-9) then
    !流物面の辺path3(ib)の一端Vt(kdc,npo(kdc,ia,ib),:)が，path2(ii)内に存在し，path2(ii)の両端とは離れている場合．
    path2(ii,3)=path2(ii,3)+1     !切断面の辺path2(ii)の交点数path2(ii,3)を登録
    ip1=0
    do ic1=1,ip
        dd3=dista(tmpV-pp(ic1,0:2))
        if(dd3.lt.1.0d-9) then
        !tmpVが既存の頂点と一致する場合
        ip1=ic1       
        endif
    enddo
    if(ip1.eq.0) then
    !tmpVが既存の頂点と一致しない場合
    ip=ip+1
call realloc_pp_info
    pp(ip,:)=tmpv
    
    i_small=min(path2(ii,1),path2(ii,2))
    i_large=max(path2(ii,1),path2(ii,2))
    
    if(i_large.le.8) then
        do ic1=1,12
            if(path0(ic1,1).eq.i_small.and.path0(ic1,2).eq.i_large) then
                info(ip)=ic1
                continue
                exit
            endif
        enddo
    endif    
    
    ip1=ip
    endif   
    icme(ii,1,path2(ii,3))=ip1  !path2(ii)のpath2(ii,3)番目の交点の頂点番号を登録 
    icme(ii,2,path2(ii,3))=npo(kdc,ia,ib)    !漂流物面の頂点であることを登録 
            if(abs(dd2).lt.1.0d-9) then
            icme(ii,3,path2(ii,3))=0            
            else if(dd2.gt.0.0d0) then
            icme(ii,3,path2(ii,3))=1
            else
            icme(ii,3,path2(ii,3))=-1            
            endif
    icv4(ib)=ip1              !Vt(kdc,npo(kdc,ia,ib),:)が切断面頂点ip1に一致することを登録
    icv41(ib)=path2(ii,-1)
    endif
endif

enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ic1=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!切断面の頂点（icv(1:icc00)）が対象物体表面に含まれるか？
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ic1=0;ic2=0;
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!切断面頂点の並びを切断面パスに登録
ipath=0            !セグ表面パスを初期化
allocate(path(0:minfopath,-2:3));path=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(ia.eq.1093) then
    continue
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!<切断面の辺path2について調べる．

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do ii=1,icc00;ii1=ii+1;if(ii1.gt.icc00) ii1=1

if(path2(ii,3).eq.0) then      !交点がない場合
    if(icv3(ii).ge.1.and.icv3(ii1).ge.1) then   !path2(ii)の両端が漂流物面の内部にある場合．
        if(icv3(ii).eq.2.and.icv3(ii1).eq.2) then 
            tmpV=cal_near_point(pp(path2(ii,1),:),pp(path2(ii,2),:),face1,length)
            call naigai2(tmpV,ipoia,vtia(1:ipoia,0:2),icc1,face1)
            if(icc1.eq.1) then
            else if(icc1.eq.0) then
                cycle
            else
            continue
            endif
        endif
        ipath=ipath+1
        path(ipath,-1:2)=path2(ii,-1:2)              !path2(ii)をセグ表面パスに登録
        icd=0
        do ib1=1,12
            if(path0(ib1,1).eq.path(ipath,1).and.path0(ib1,2).eq.path(ipath,2)) icd=ib1
            if(path0(ib1,2).eq.path(ipath,1).and.path0(ib1,1).eq.path(ipath,2)) icd=ib1
        enddo
        if(icd.eq.0) then
        !path(ipath,-1)=ib  
        else
          path(ipath,-1)=icd+100   !格子面の辺番号に100を加えて登録(←面番号と区別するため)
        endif            
    endif 

else if(path2(ii,3).eq.1) then                     !交点が一つある場合
        ip1=icme(ii,1,path2(ii,3))                     !交点の頂点番号をip1に代入
        
    if(icv3(ii).eq.1) then
    else if(icv3(ii).eq.2) then
        tmpV=cal_near_point(pp(path2(ii,1),:),pp(ip1,:),face1,length)
        call naigai2(tmpV,ipoia,vtia(1:ipoia,0:2),icc1,face1)  
        if(icc1.eq.1) then
            continue
        else if(icc1.eq.0) then
            goto 2134
        else
            continue
        endif
    else
       goto 2134
    endif
        
            ipath=ipath+1
            path(ipath,1)=path2(ii,1)  !
            path(ipath,2)=ip1          !
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．             
        
        
2134    if(icv3(ii1).eq.1) then  
                
        else if(icv3(ii1).eq.2) then
            tmpV=cal_near_point(pp(ip1,:),pp(path2(ii,2),:),face1,length)
            call naigai2(tmpV,ipoia,vtia(1:ipoia,0:2),icc1,face1) 
            if(icc1.eq.1) then
                continue
            else if(icc1.eq.0) then
                cycle
            else
                continue
            endif
        else 
            cycle
        endif
            
           ipath=ipath+1
            path(ipath,1)=ip1 
            path(ipath,2)=path2(ii,2)  !
          !
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．      
        
 
#ifdef DDDD        
        if(icc1.ge.1.and.icc2.eq.0) then
        !icv(ii)が物体表面内，icv(ii1)が物体表面外のとき，
            ipath=ipath+1
            path(ipath,1)=path2(ii,1)  !
            path(ipath,2)=ip1          !
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．
        else if(icc1.eq.0.and.icc2.ge.1) then
        !icv(ii)が物体表面外，icv(ii1)が物体表面内のとき，
            ipath=ipath+1
            path(ipath,1)=ip1
            path(ipath,2)=path2(ii,2)
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．            
        else if(icv3(ii).eq.2.and.icv3(ii1).eq.2) then
            
            if(icc1.eq.1) then
            ipath=ipath+1
            path(ipath,1)=path2(ii,1)  !
            path(ipath,2)=ip1          !
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．      
            endif
            if(icc2.eq.1) then
            ipath=ipath+1
            path(ipath,1)=ip1
            path(ipath,2)=path2(ii,2)
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．
            endif
        else
            continue    
        endif
#endif        
        
else if(path2(ii,3).eq.2) then                     !交点が二つある場合
    
            dd1=dista(pp(path2(ii,1),:)-pp(icme(ii,1,1),:))
            dd2=dista(pp(path2(ii,1),:)-pp(icme(ii,1,2),:))
    
            if(dd1.lt.dd2) then
                icme1(1:3,1:2)=icme(ii,1:3,1:2)
            else 
                icme1(1:3,1)=icme(ii,1:3,2);icme1(1:3,2)=icme(ii,1:3,1)
            endif   
            
if(icv3(ii).ge.1) then    
    tmpV=cal_near_point(pp(path2(ii,1),:),pp(icme1(1,1),:),face1,length)
    call naigai2(tmpV,ipoia,vtia(1:ipoia,0:2),icc1,face1) 
    if(icc1.eq.1) then
            ipath=ipath+1
            path(ipath,1)=path2(ii,1)
            path(ipath,2)=icme1(1,1)  
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．  
    endif    
endif
    tmpV=cal_near_point(pp(icme1(1,1),:),pp(icme1(1,2),:),face1,length)
    call naigai2(tmpV,ipoia,vtia(1:ipoia,0:2),icc1,face1) 
if(icc1.eq.1) then
            ipath=ipath+1
            path(ipath,1:2)=icme1(1,1:2)
            path(ipath,-1:0)=path2(ii,-1:0) !path2の面情報をpathに写す． 
endif 
if(icv3(ii1).ge.1) then
    tmpV=cal_near_point(pp(icme1(1,2),:),pp(path2(ii,2),:),face1,length)
    call naigai2(tmpV,ipoia,vtia(1:ipoia,0:2),icc1,face1) 
    if(icc1.eq.1) then
            ipath=ipath+1
            path(ipath,1)=icme1(1,2)  
            path(ipath,2)=path2(ii,2)
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す． 
    endif
endif
    
else if(path2(ii,3).eq.3) then 
            dd1=dista(pp(path2(ii,1),:)-pp(icme(ii,1,1),:))
            dd2=dista(pp(path2(ii,1),:)-pp(icme(ii,1,2),:))
            dd3=dista(pp(path2(ii,1),:)-pp(icme(ii,1,3),:))
            if(dd1.lt.dd2.and.dd2.lt.dd3) then
                icme1(1:3,1:3)=icme(ii,1:3,1:3);
            else if(dd1.lt.dd3.and.dd3.lt.dd2) then
                icme1(1:3,1)=icme(ii,1:3,1);
                icme1(1:3,2)=icme(ii,1:3,3);
                icme1(1:3,3)=icme(ii,1:3,2)
            else if(dd2.lt.dd1.and.dd1.lt.dd3) then    
                icme1(1:3,1)=icme(ii,1:3,2);
                icme1(1:3,2)=icme(ii,1:3,1);
                icme1(1:3,3)=icme(ii,1:3,3)        
            else if(dd2.lt.dd3.and.dd3.lt.dd1) then 
                icme1(1:3,1)=icme(ii,1:3,2);
                icme1(1:3,2)=icme(ii,1:3,3);
                icme1(1:3,3)=icme(ii,1:3,1)
            else if(dd3.lt.dd1.and.dd1.lt.dd2) then  
                icme1(1:3,1)=icme(ii,1:3,3);
                icme1(1:3,2)=icme(ii,1:3,1);
                icme1(1:3,3)=icme(ii,1:3,2)         
            else if(dd3.lt.dd2.and.dd2.lt.dd1) then  
                icme1(1:3,1)=icme(ii,1:3,3);
                icme1(1:3,2)=icme(ii,1:3,2);
                icme1(1:3,3)=icme(ii,1:3,1)         
            else
                continue
            endif             

if(icv3(ii).eq.1) then
            ipath=ipath+1
            path(ipath,1)=path2(ii,1)
            path(ipath,2)=icme1(1,1)  
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．  
endif
if(icme1(3,1).eq.1) then
            ipath=ipath+1
            path(ipath,1:2)=icme1(1,1:2)
            path(ipath,-1:0)=path2(ii,-1:0) !path2の面情報をpathに写す． 
endif 
if(icme1(3,2).eq.1) then
            ipath=ipath+1
            path(ipath,1:2)=icme1(1,2:3)
            path(ipath,-1:0)=path2(ii,-1:0) !path2の面情報をpathに写す． 
endif 

if(icv3(ii1).eq.1) then
            ipath=ipath+1
            path(ipath,1)=icme1(1,3)  
            path(ipath,2)=path2(ii,2)
            path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す． 
endif            
  
else if(path2(ii,3).eq.4) then  
            dd4(1)=dista(pp(path2(ii,1),:)-pp(icme(ii,1,1),:))
            dd4(2)=dista(pp(path2(ii,1),:)-pp(icme(ii,1,2),:))
            dd4(3)=dista(pp(path2(ii,1),:)-pp(icme(ii,1,3),:))  
            dd4(4)=dista(pp(path2(ii,1),:)-pp(icme(ii,1,4),:))      
            idd4 =  minloc(dd4(1:4));idd5(1)=icme(ii,1,idd4(1));dd4(idd4(1))=10000.d0     
            idd4 =  minloc(dd4(1:4));idd5(2)=icme(ii,1,idd4(1));dd4(idd4(1))=10000.d0  
            idd4 =  minloc(dd4(1:4));idd5(3)=icme(ii,1,idd4(1));dd4(idd4(1))=10000.d0  
            idd4 =  minloc(dd4(1:4));idd5(4)=icme(ii,1,idd4(1));dd4(idd4(1))=10000.d0   
             if(icv3(ii).eq.0.and.icv3(ii1).eq.0) then   
                ipath=ipath+1 
                path(ipath,1:2)=idd5(1:2)  !        !
                path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．                
                ipath=ipath+1 
                path(ipath,1:2)=idd5(3:4)  !       !
                path(ipath,-1:0)=path2(ii,-1:0)  !path2の面情報をpathに写す．                                
              continue 
             else
              write(6,*) 'soutei gai sshape Line10110' 
              stop 
             endif           
else
              write(6,*) 'soutei gai sshape Line1010'
              stop
endif
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!<漂流物面上の辺path3について調べる．

do ii=1,ipo(kdc,ia)
ii1=ii+1;if(ii1.gt.ipo(kdc,ia)) ii1=1
if(path3(ii,3).eq.0) then  !交点がない場合
         if(icv4a(ii).ge.1.and.icv4a(ii1).ge.1) then !path3の両端が切断面内の場合，
            tmpv=HALF*(pp(icv4(ii),:)+pp(icv4(ii1),:))
            call naigai2(tmpV,icc00,pp(icv(1:icc00),0:2),icc1,face1)  
            if(icc1.eq.1) then
                ipath=ipath+1
                path(ipath,1)=icv4(ii)
                path(ipath,2)=icv4(ii1)
                path(ipath,-1)=path3(ii,-1)   
            endif
         endif
else if(path3(ii,3).eq.1) then  !交点がひとつの場合
        ip1=icdr(ii,path3(ii,3))
        if(icv4a(ii).ge.1.and.icv4a(ii1).eq.0) then
             tmpv=HALF*(pp(ip1,:)+pp(icv4(ii),:))
             call naigai2(tmpV,icc00,pp(icv(1:icc00),0:2),icc1,face1)  
             if(icc1.eq.1) then
!             icv4(ii)が物体表面外，icv4(ii1)が物体表面内のとき，
                ipath=ipath+1
                path(ipath,1)=icv4(ii)   !!
                path(ipath,2)=ip1        !!
                path(ipath,-1)=path3(ii,-1)             
             else
                 continue
             endif
        !!     
        else if(icv4a(ii).eq.0.and.icv4a(ii1).ge.1) then
        !     tmpv=HALF*(pp(ip1,:)+pp(path2(ii,1),:))
        !      call naigai2(tmpV,ipoia,vtia(1:ipoia,0:2),icc1,face1)
            tmpv=HALF*(pp(ip1,:)+pp(icv4(ii1),:))
            call naigai2(tmpV,icc00,pp(icv(1:icc00),0:2),icc2,face1)   
              if(icc2.eq.1) then            
        !icv4(ii)が物体表面外，icv4(ii1)が物体表面内のとき，
            ipath=ipath+1            !!
            path(ipath,1)=ip1        !!
            path(ipath,2)=icv4(ii1)  !!
            path(ipath,-1)=path3(ii,-1)         !!   
        !    else
        !        continue
              endif
  
        else
            continue    
        endif
    continue
else if(path3(ii,3).eq.2) then !交点が二つの場合
        if(icv4(ii).eq.0.and.icv4(ii1).eq.0) then
        !icv4(ii)，icv4(ii1)がともに物体表面外のとき，
            ipath=ipath+1
            path(ipath,1:2)=icdr(ii,1:2)
 !           path(ipath,2)=icdr(ii,2)
                     !  path(ipath,-1)=0   
        else if(icv4(ii).eq.0.and.icv4(ii1).ge.1) then
        !icv4(ii)が物体表面外，icv4(ii1)が物体表面内のとき，，

            dd1=dista(pp(icdr(ii,1),:)-pp(icv4(ii1),:))
            dd2=dista(pp(icdr(ii,2),:)-pp(icv4(ii1),:))
            if(dd1.lt.dd2) then   !交点icdr(ii,1)がicv4(ii1)に近いとき，これを両端とするパスを登録
                        ipath=ipath+1
                path(ipath,1)=icv4(ii1)
                path(ipath,2)=icdr(ii,1) 
                      !      path(ipath,-1)=0  
            else
                        ipath=ipath+1
                path(ipath,1)=icdr(ii,2)
                path(ipath,2)=icv4(ii1) 
                        !    path(ipath,-1)=0  
            endif 
                       ipath=ipath+1
                path(ipath,1:2)=icdr(ii,1:2)
 !               path(ipath,2)=icdr(ii,2) 
                       !     path(ipath,-1)=0  
        else if(icv4(ii).ge.1.and.icv4(ii1).eq.0) then
        !icv4(ii)が物体表面外，icv4(ii1)が物体表面内のとき，
            dd1=dista(pp(icdr(ii,1),:)-pp(icv4(ii),:))
            dd2=dista(pp(icdr(ii,2),:)-pp(icv4(ii),:))        
            if(dd1.lt.dd2) then !交点icdr(ii,1)がicv4(ii)に近いとき，これを両端とするパスを登録
                        ipath=ipath+1
                path(ipath,1)=icv4(ii)
                path(ipath,2)=icdr(ii,1) 
                      !      path(ipath,-1)=0
            else
                        ipath=ipath+1
                path(ipath,1)=icdr(ii,2)
                path(ipath,2)=icv4(ii) 
                     !       path(ipath,-1)=0  
            endif         
            !交点icdr(ii,1),icdr(ii,2)を両端とするパスを登録
            ipath=ipath+1
            path(ipath,1:2)=icdr(ii,1:2)
!            path(ipath,2)=icdr(ii,2)               
                    !    path(ipath,-1)=0          
        else
            continue    
        endif
else 
    continue
endif
enddo

!if(ikouten.eq.0) then
!if(ic1+ic3.eq.icc.and.in_vt.eq.0) cycle !格子切断面のすべての頂点が漂流物面の外側の場合，終了
!物体表面と切断面表面との交点がない場合
!else
!物体表面と切断面表面との交点がある場合，セグ表面パスを連続させる．

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!pathを並べ替える．
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(ipath.lt.3) then
!pathが2以下の時,面が成立しない．
 continue
 cycle
endif
continue
#ifdef DDDDD
#endif
!pathが重複していないか調べる．
2314 do i=1,ipath
    if(path(i,3).eq.1) cycle
    ib1=0
    do ib=1,ipath
        if(i.eq.ib) cycle
        if(path(i,1).eq.path(ib,1).and.path(i,2).eq.path(ib,2)) then
            ib1=ib;exit
!        path(ib,3)=1   !重複していたら使用済みパスにする．
        else if(path(i,1).eq.path(ib,2).and.path(i,2).eq.path(ib,1)) then
!        path(ib,3)=1   !重複していたら使用済みパスにする．
            ib1=ib;exit
        endif        
    enddo
    if(ib1.ne.0) exit
enddo
if(ib1.ne.0) then
    path(ib1:ipath-1,:)=path(ib1+1:ipath,:)
    ipath=ipath-1
    goto 2314
endif

if(ipath.lt.3) then
!pathが2以下の時,面が成立しない．
 continue
 cycle
endif


!if(nn_1.eq.63510) then
!write(6,*) 'nn_1=',nn_1,ipath
!endif

2371 path(1:ipath,3)=0

icyc=0
loop2=0
236 icv(1)=path(1,1)
ic4=1
loop=0
235 do i=1,ipath
if(path(i,3).eq.1) cycle
if(icv(ic4).eq.path(i,1)) then
    path(i,3)=1
    ic4=ic4+1
    icv(ic4)=path(i,2)
    continue
else if(icv(ic4).eq.path(i,2)) then
    path(i,3)=1
    ic4=ic4+1
    icv(ic4)=path(i,1)
    continue
endif
enddo
            loop=loop+1
        if(loop.gt.100) then  !百回うまくいかない場合．
            path(1:ipath,3)=0  !使用済み判定を初期化
            ic1=0
            !pathが重複していないか調べる．
            do i=1,ipath
                if(path(i,3).eq.1) cycle
            do ib=1,ipath
            if(i.eq.ib) cycle
            if(path(i,1).eq.path(ib,1).and.path(i,2).eq.path(ib,2)) then
            path(ib,3)=1   !重複していたら使用済みパスにする．
        !    path(ib,-1)=max(path(ib,-1),path(i,-1))
        !    path(i,-1)=path(ib,-1)
            ic1=1
            continue
            endif
            if(path(i,1).eq.path(ib,2).and.path(i,2).eq.path(ib,1)) then
            path(ib,3)=1   !重複していたら使用済みパスにする．
            ic1=1    
         !   path(ib,-1)=max(path(ib,-1),path(i,-1))
         !  path(i,-1)=path(ib,-1)    
            continue
            endif        
            enddo
            enddo
            
            if(ic1.eq.1) then
                loop2=loop2+1
                
                if(loop2.gt.100) then
                            write(6,*) 'loop2=100,nn=',nn,'235','i,k,j=',in(nn),kn(nn),jn(nn)
                            goto 932
                endif
                !重複していたらやり直し．
                goto 236
            endif
            
            if(pp(1,0)-xxmax(kdc,ia,0).ge.0.0d0.or.pp(1,1)-xxmax(kdc,ia,1).ge.0.0d0.or.pp(1,2)-xxmax(kdc,ia,2).ge.0.0d0.or.&
   pp(7,0)-xxmin(kdc,ia,0).le.0.0d0.or.pp(7,1)-xxmin(kdc,ia,1).le.0.0d0.or.pp(7,2)-xxmin(kdc,ia,2).le.0.0d0) then
                 continue
                 cycle
   endif
   dd1=0
    do i=2,ipath
            dd1=crossL(pp(path(1,1),:)-pp(path(1,2),:),pp(path(i,1),:)-pp(path(i,2),:))+dd1
    enddo
    if(abs(dd1).lt.1.0d-10) then
        cycle
    else
        write(6,*) 'dd1=',dd1
    endif
     call checkdata2ia
            write(6,*) 'loop=1000,nn=',nn,'235','i,k,j=',in(nn),kn(nn),jn(nn),'ia=',ia
             cycle
#ifdef BUGWIN
        call checkdata
  !      call checkdata2
   !     call checkdata3 
   path(1:ipath,3)=0

icc1=0
icc3=0
do ia1=1,ip
    icc=0
    do ib=1,ipath
  !      if(path(ib,3).eq.1)  cycle       
        if(path(ib,1).eq.ia1) icc=icc+1
        if(path(ib,2).eq.ia1) icc=icc+1
    enddo
    if(icc.gt.2) then
        icc1=icc1+1
        icv1(icc1)=ia1
        icv2(icc1)=icc
    else if(icc.eq.1) then
        icc3=icc3+1
        icv3(icc3)=ia1
    endif
enddo

if(icc3.gt.0) then  !ひげをとる。
    do ib=1,ipath
 !       if(path(ib,3).eq.0)  cycle       
        if(path(ib,1).eq.icv3(1)) then
            path(ib,3)=1
        else if(path(ib,2).eq.icv3(1)) then
            path(ib,3)=1
        endif
    enddo
    goto 236
endif
        continue
              cycle
    !           icv(1)=13;icv(2)=12;icv(3)=17;icc=3;goto 687
#else
             goto 932
         call checkdata2
        call checkdata3  
              stop
#endif
        endif
if(icv(ic4).ne.icv(1).or.ic4.le.3) goto 235
!現在の頂点が始点に一致しない場合，２３５へ


continue
icc=ic4-1 !頂点数をiccに登録
!endif

687 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(inear.ne.0) then
  !  if(maxval(icv(1:icc)).le.8) then
     inear_c=inear_c+1
  !  else
     cycle
  !  endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!メッシュ内漂流物面情報を登録．!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!f!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mmc=mmc+1                             !セルnn内のセグ表面番号mmc
icpl=icpl+1                           !セルnn内の総表面番号icpl
cpl2(icpl,:)=(-1.0d0)*face1           !面icplの面情報を登録
ncpl2(icpl,-1)=1                      !面icplを物体表面として登録
ncpl2(icpl,-2)=idri_fc_all_n(kdc,ia)  !面icplの全体表面番号を登録．
ncpl2(icpl,0)=icc                     !面icplの頂点数を登録．
ncpl2(icpl,1:icc)=icv(1:icc)
call surface_dir3(icc,icv(1:icc),cpl2(icpl,0:2),ncpl2(icpl,1:icc)) !面icplの頂点並びを登録

if(idri_fc_all_n(kdc,ia).eq.317980) then
    continue
endif
!!!!!!　辺「ncpl2(icpl2,ib)-ncpl2(icpl2,ib1)」が属するセル境界面番号を設定
do ib=1,icc
    ib1=ib+1
    if(ib1.gt.icc) ib1=1
    ncpl2(icpl,ib+icc)=0
    do ic1=1,ipath
    if(path(ic1,1).eq.ncpl2(icpl,ib).and.path(ic1,2).eq.ncpl2(icpl,ib1).and.ncpl2(icpl,ib+icc).eq.0) then
        if( path(ic1,0) .gt.100) then
            ncpl2(icpl,ib+icc)=path(ic1,0)          
        else
            ncpl2(icpl,ib+icc)=path(ic1,-1)
        endif
    else if(path(ic1,2).eq.ncpl2(icpl,ib).and.path(ic1,1).eq.ncpl2(icpl,ib1).and.ncpl2(icpl,ib+icc).eq.0) then 
        if( path(ic1,0) .gt.100) then
            ncpl2(icpl,ib+icc)=path(ic1,0)        
        else    
            ncpl2(icpl,ib+icc)=path(ic1,-1)
        endif
    endif
    enddo
enddo
!if(nn_1.eq.63510) then
!write(6,*) 'path',path(1:ipath,1)
!endif

mmd_p(kdc)=mmd_p(kdc)+1
#ifdef TEST
mmd(kdc)=mmd(kdc)+1       ! セグメント表面通し番号
call realloc_IDP_ND_DP(kdc)
IDP(kdc,mmd(kdc),1)=nn    ! セル番号
IDP(kdc,mmd(kdc),2)=idri_fc_n(ncpl2(icpl,-2))     ! 面番号
!         SDP(kdc,mmd(kdc))=S(i)    ! 面積             
! 頂点数，頂点座標を整理する
do kk=1,icc
dp(kdc,mmd(kdc),kk,1:3)=pp(ncpl2(icpl,icc-kk+1),0:2)
enddo
nd(kdc,mmd(kdc))=icc                     
#endif TEST
continue

!デバッグライトstart
#ifdef BUGWIN2
if(cinfn(nn,kdc).ne.0) then
write(6,*) icv(1:icc)
write(6,*) 'nn=',nn
write(6,*) 'face=',face1
do kk=1,icc
write(6,*) pp(icv(kk),:)
enddo
write(6,*) 'upl=',upl(kdc,IDP(kdc,cinf(nn,kdc,mmc),2),:)
do kk=1,nd(kdc,cinf(nn,kdc,mmc))
write(6,*) dp(kdc,cinf(nn,kdc,mmc),kk,1:3)
enddo
endif
#endif
!デバッグライトend

if(ic4.lt.ipath-1.and.icyc.eq.0) then
icyc=1
do i=1,ipath
if(path(i,3).eq.1) cycle
continue
icv(1)=path(i,1)
exit
enddo
ic4=1
goto 235
endif
deallocate(path)

enddo ia_loop ! !物体kdcのすべての物体表面iaについて繰り返す．

 if(inear_c.gt.0.and.icpl.eq.6) then  
    !    nff(nn)%f=-1!    nff(nn)%b=10!    fb0(nn)=0.0d0
        axyz=0.0d0
        drino2(nn,2)=kdc   
        return
endif


!!!!!!!CINFN(nn,kdc)=mmc    ! セル内セグメント数の登録
enddo kdc_loop ! kdc

!f(nn_1.eq.63500.and.flag.eq.0) then
 !write(66,*) 'end _ia nn_1=',nn_1,nn
!endif

!セルnnの表面数を登録
ncpl2(0,0)=icpl

if(objno_nn.ne.0) then
    continue
endif


ipathu=0
pathu=0
!write(6,*) 'icpl=',icpl,nn

#ifdef BUGWIN
if(nn.eq.46434) then
continue
endif
!if(ncpl(0,0).ne.6) then
!continue
!endif
#endif

if(flag.eq.2) then
if(i_in.eq.0.and.icpl.eq.ncplc00) then
    continue
#ifdef DRI  
        mm11=mm11+1
        drino2i(mm11)=nn
        drino2(nn,1)=mm11;drino2(nn,2)=kdc   
#else
        axyz=0.0d0
#endif  
endif

endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(icpl.eq.ncplc00) return  !cycle                                                 !!!!!! セルnn内部に物体表面がなかった時　次のセルnn+1 に移る　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DDDDDD
if(flag.eq.1.and.OBJNO(nn).gt.0) then !固定物体も含まれる場合は、・・・・・・
    ii=ncpl0(OBJNO(nn),0,0)
    do ii=1,ncpl0(OBJNO(nn),0,0)
    if(ncpl0(OBJNO(nn),ii,-1).ne.1) cycle
    icpl=icpl+1
    icc=ncpl0(OBJNO(nn),ii,0)
    ncpl2(icpl,-2:icc+icc)=ncpl0(OBJNO(nn),ii,-2:icc+icc)
    cpl2(icpl,:)=cpl(ncpl2(icpl,-2),:)
    do ia=1,ncpl0(OBJNO(nn),ii,0)
        ic1=0
        do ib=1,ip
            dd1=dista(pp_0(OBJNO(nn),ncpl2(icpl,ia),0:2)-pp(ib,0:2))
            if(dd1.lt.1.0d-9) then
            ncpl2(icpl,ia)=ib;ic1=1
            endif
            enddo 
            if(ic1.eq.0) then
                ip=ip+1
                call realloc_pp_info
                pp(ip,:)=pp_0(OBJNO(nn),ncpl2(icpl,ia),0:2)
                info(ip)=info_0(OBJNO(nn),1,ncpl2(icpl,ia))
                ncpl2(icpl,ia)=ip          
            endif
        enddo          
    enddo
    ncpl2(0,0)=icpl
endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!物体表面が存在する場合
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(allocated(info3)) deallocate(info3)
allocate(info3(1:ip,0:10))
info3=0

do i=ncplc00+1,icpl
    
if(nn.eq.53914) then
   ! write(6,*) i,info2(30)
    continue
endif
    continue
    
!面icplの頂点並びをpathu（空隙部分用）,patho（物体部分用）に登録
do ib=1,ncpl2(i,0)
    if(info(ncpl2(i,ib)).gt.0.and.info(ncpl2(i,ib)).lt.100) then
#ifdef DDDDD
        if(info2(ncpl2(i,ib)).eq.0) then
            info2(ncpl2(i,ib))=i
        else
            if(info2(ncpl2(i,ib)).lt.10000) then
                if(sum((cpl2(info2(ncpl2(i,ib)),:)-cpl2(i,:))**2).ne.0.0d0) then
                info2(ncpl2(i,ib))=info2(ncpl2(i,ib))+10000*i
                else
                continue
                endif
            else if(info2(ncpl2(i,ib)).lt.100000000) then
                info2(ncpl2(i,ib))=info2(ncpl2(i,ib))+i*100000000  
            else if(info2(ncpl2(i,ib)).lt.1000000000000) then
                info2(ncpl2(i,ib))=info2(ncpl2(i,ib))+i*1000000000000                           
            else 
                continue   
            endif
        endif
#endif        
        if(info3(ncpl2(i,ib),0).eq.0) then
            info3(ncpl2(i,ib),1)=i;info3(ncpl2(i,ib),0)=1
        else
            if(info3(ncpl2(i,ib),0).eq.1) then
             !   if(sum((cpl2(info3(ncpl2(i,ib),1),:)-cpl2(i,:))**2).ne.0.0d0) then
                info3(ncpl2(i,ib),2)=i;info3(ncpl2(i,ib),0)=2
              !  else
               ! continue
             !   endif
            else if(info3(ncpl2(i,ib),0).eq.2) then
                info3(ncpl2(i,ib),3)=i;info3(ncpl2(i,ib),0)=3
            else if(info3(ncpl2(i,ib),0).eq.3) then
                info3(ncpl2(i,ib),4)=i;info3(ncpl2(i,ib),0)=4   
            else if(info3(ncpl2(i,ib),0).eq.4) then
                info3(ncpl2(i,ib),5)=i;info3(ncpl2(i,ib),0)=5                                       
            else 
                continue   
            endif
        endif        
    else
        continue
    endif
    
    
    ib1=ib+1;if(ib1.gt.ncpl2(i,0)) ib1=1
    !既登録を確認
    iia=0
    do ii=1,ipathu
    if(pathu(ii,1).eq.ncpl2(i,ib).and.pathu(ii,2).eq.ncpl2(i,ib1)) iia=1
    if(pathu(ii,2).eq.ncpl2(i,ib).and.pathu(ii,1).eq.ncpl2(i,ib1)) iia=1
    enddo

    if(iia.eq.0) then
    !既登録がなければ,pathu（空隙部分用）patho（物体部分用）にパス登録
    ipathu=ipathu+1
    pathu(ipathu,1)=ncpl2(i,ib)
    pathu(ipathu,2)=ncpl2(i,ib1)
    pathu(ipathu,-1)=ncpl2(i,ib+ncpl2(i,0))    !該当パスの属する面番号を登録
    pathu(ipathu,-2)=ncpl2(i,-2)    
continue    
    ipatho=ipatho+1
    patho(ipatho,1)=ncpl2(i,ib)
    patho(ipatho,2)=ncpl2(i,ib1)
!    patho(ipatho,3)=idri_kdc(idri_fc_all_n(kdc,ia))
    patho(ipatho,3)=idri_kdc(ncpl2(i,-2))
    patho(ipatho,-1)=ncpl2(i,ib+ncpl2(i,0))  
    continue
    endif
enddo  

enddo



#ifdef BUGWIN
if(ncpl2(0,0).gt.ncplc00+1) then
    continue
endif    
if(nn.eq.29098) then
continue
endif    
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ncpl2(ib,0)に登場しない頂点は，info=-1とする．
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ia=1,ip
if(info(ia).eq.0) cycle
if(info(ia).eq.-1) cycle
do ib=ncplc00+1,ncpl2(0,0)
    do ib1=1,ncpl2(ib,0)
        if(ia.eq.ncpl2(ib,ib1)) goto 239
    enddo
enddo
info(ia)=-1
239 continue
enddo
    
info(1:8)=0

do ia=1,8
do ib=ncplc00+1,ncpl2(0,0)
    do ib1=1,ncpl2(ib,0)
        if(ia.eq.ncpl2(ib,ib1)) goto 2361
    enddo
enddo
cycle
2361 continue
info(ia)=1
enddo


!do i=1,ipathu
!write(6,*) i,pathu(i,1:2),'pathaaaaa1'
!enddo
if(allocated(info1))  deallocate(info1)
allocate(info1(1:ip))
info1(1:ip)=info(1:ip)
sinfo(1:12)=1
!info(6)=0
Cell_Sides: do i=1,12   !格子の辺について
icc=0;icc2=0;icc3=0;icc4=0
do ia=9,ip  !すべての点について
    if(info(ia).eq.i) then
        icc=icc+1
        ipcc(icc)=ia

    endif
#ifdef DDDDD
    if(info(ia).eq.i.and.icc.ne.0.and.icc2.ne.0.and.icc3.ne.0..and.icc4.eq.0) then !辺番号iと一致した場合，点番号iaをicc4に登録
    continue
    icc4=ia
    endif
    if(info(ia).eq.i.and.icc.ne.0.and.icc2.ne.0.and.icc3.eq.0) then !辺番号iと一致した場合，点番号iaをicc3に登録
    continue
    icc3=ia
    endif
    if(info(ia).eq.i.and.icc.ne.0.and.icc2.eq.0) then !辺番号iと一致した場合，点番号iaをicc2に登録
    continue
    icc2=ia
    endif
    if(info(ia).eq.i.and.icc.eq.0) then !辺番号iと一致した場合，点番号iaをiccに登録
    continue
    icc=ia        
    endif
#endif
enddo

if(icc.eq.0) cycle !辺上に交点なし．
sinfo(i)=0

if(icc.eq.1) then !直線との交点が１つの場合。
!
    ip1=ipcc(icc)
! path0(i,1)と ip1について
    call cal_pathu2(ip1,path0(i,1),path0(i,2))
! path0(i,2)と ip1について
    call cal_pathu2(ip1,path0(i,2),path0(i,1))
    
    cycle
endif

do ii=1,icc
        dpcc(ii)=dista(pp(ipcc(ii),0:2)-pp(path0(i,1),0:2))
enddo
ii=0
do ii=1,icc
!write(6,*) MINLOC(dpcc(1:icc))
ML(1:icc)=MINLOC(dpcc(1:icc))
ipccn(ii)=ipcc(ML(1));dpcc(ML(1))=100000.d0
enddo
ipcc(1:icc)=ipccn(1:icc)

 ! path0(i,1)と ip1について
    call cal_pathu2(ipcc(1),path0(i,1),path0(i,2))
! ip1と ip2について
do ii=1,icc-1 
    call cal_pathu2(ipcc(ii+1),ipcc(ii),path0(i,2))     
enddo
! path0(i,2)と ip4について
    call cal_pathu2(ipcc(icc),path0(i,2),path0(i,1))   

continue
cycle


if(icc4.ne.0) then !直線との交点が４つの場合。
    tt1=dista(pp(icc,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点iccとの距離           
    tt2=dista(pp(icc2,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点icc2との距離 
    tt3=dista(pp(icc3,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点icc3との距離  
    tt4=dista(pp(icc4,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点icc4との距離
    if(tt1.lt.min(tt2,tt3,tt4)) then
        ip1=icc;tt1=10000.d0
    else if(tt2.lt.min(tt1,tt3,tt4)) then
        ip1=icc2;tt2=10000.d0
    else if(tt3.lt.min(tt1,tt2,tt4)) then
        ip1=icc3;tt3=10000.d0
    else if(tt4.lt.min(tt2,tt3,tt1)) then
        ip1=icc4;tt4=10000.d0
    endif   
    
    if(tt1.lt.min(tt2,tt3,tt4)) then
        ip2=icc;tt1=10000.d0
    else if(tt2.lt.min(tt1,tt3,tt4)) then
        ip2=icc2;tt2=10000.d0
    else if(tt3.lt.min(tt1,tt2,tt4)) then
        ip2=icc3;tt3=10000.d0
    else if(tt4.lt.min(tt2,tt3,tt1)) then
        ip2=icc4;tt4=10000.d0
    endif  
    
    if(tt1.lt.min(tt2,tt3,tt4)) then
        ip3=icc;tt1=10000.d0
    else if(tt2.lt.min(tt1,tt3,tt4)) then
        ip3=icc2;tt2=10000.d0
    else if(tt3.lt.min(tt1,tt2,tt4)) then
        ip3=icc3;tt3=10000.d0
    else if(tt4.lt.min(tt2,tt3,tt1)) then
        ip3=icc4;tt4=10000.d0
    endif  
    
     if(tt1.lt.min(tt2,tt3,tt4)) then
        ip4=icc;tt1=10000.d0
    else if(tt2.lt.min(tt1,tt3,tt4)) then
        ip4=icc2;tt2=10000.d0
    else if(tt3.lt.min(tt1,tt2,tt4)) then
        ip4=icc3;tt3=10000.d0
    else if(tt4.lt.min(tt2,tt3,tt1)) then
        ip4=icc4;tt4=10000.d0
    endif         
        
       
    
 ! path0(i,1)と ip1について
    call cal_pathu2(ip1,path0(i,1),path0(i,2))
! ip1と ip2について
    call cal_pathu2(ip2,ip1,path0(i,2))   
! ip2と ip3について
    call cal_pathu2(ip3,ip2,path0(i,2))    
! ip3と ip4について
    call cal_pathu2(ip4,ip3,path0(i,2))       
! path0(i,2)と ip4について
    call cal_pathu2(ip4,path0(i,2),path0(i,1))   
    
        
else if(icc3.ne.0) then !直線との交点が３つの場合。
    tt1=dista(pp(icc,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点iccとの距離           
    tt2=dista(pp(icc2,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点icc2との距離 
    tt3=dista(pp(icc3,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点icc3との距離    
    if(tt1.lt.min(tt2,tt3)) then
        ip1=icc
        if(tt2.lt.tt3) then
        ip2=icc2;ip3=icc3
        else
        ip2=icc3;ip3=icc2       
        endif
    else if(tt2.lt.min(tt1,tt3)) then
        ip1=icc2
        if(tt1.lt.tt3) then
        ip2=icc;ip3=icc3
        else
        ip2=icc3;ip3=icc       
        endif        
    else
        ip1=icc3
        if(tt1.lt.tt2) then
        ip2=icc;ip3=icc2
        else
        ip2=icc2;ip3=icc       
        endif   
    endif
    
! path0(i,1)と ip1について
    call cal_pathu2(ip1,path0(i,1),path0(i,2))
! ip1と ip2について
    call cal_pathu2(ip2,ip1,path0(i,2))   
! ip2と ip3について
    call cal_pathu2(ip3,ip2,path0(i,2))      
! path0(i,2)と ip3について
    call cal_pathu2(ip3,path0(i,2),path0(i,1))
!
    continue
!    
else if(icc2.ne.0) then !直線との交点が２つの場合。
    tt1=dista(pp(icc,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点iccとの距離           
    tt2=dista(pp(icc2,0:2)-pp(path0(i,1),0:2)) !辺の一端path0(i,1)と交点icc2との距離 
    if(tt1.lt.tt2) then
        ip1=icc;ip2=icc2
    else
        ip1=icc2;ip2=icc    
    endif
! path0(i,1)と ip1について
    call cal_pathu2(ip1,path0(i,1),path0(i,2))
! ip1と ip2について
    call cal_pathu2(ip2,ip1,path0(i,2))   
! path0(i,2)と ip2について
    call cal_pathu2(ip2,path0(i,2),path0(i,1))
!       

!
!
endif



continue

enddo Cell_Sides
!do i=1,ipathu
!write(6,*) i,pathu(i,1:2),'pathaaaaa2'
!enddo
!write(6,*) info(7),info(8)
continue


2345 icc=0
do i=1,ip0
if(info(i).ne.-1) cycle
do ia=1,12
if(path0(ia,1).eq.i) then
    if(sinfo(ia).eq.1.and.info(path0(ia,2)).eq.0) then
    info(path0(ia,2))=-1
    icc=icc+1
    endif
    continue
else if(path0(ia,2).eq.i) then
    if(sinfo(ia).eq.1.and.info(path0(ia,1)).eq.0) then
    info(path0(ia,1))=-1
    icc=icc+1
    endif
    continue
endif
enddo
enddo
if(icc.ne.0) goto 2345

do i=1,ip0
if(info(i).ne.1) cycle
do ia=1,12
if(path0(ia,1).eq.i) then
    if(sinfo(ia).eq.1.and.info(path0(ia,2)).eq.0) then
        if(ip.eq.ip0.and.icpl.eq.7) then
            dd1=sum(cpl2(icpl,0:2)*pp(path0(ia,2),:))+cpl2(icpl,3)
            if(dd1.gt.0.0d0) then
             info(path0(ia,2))=-1
            else
             info(path0(ia,2))=1             
            endif
            continue
        else
            info(path0(ia,2))=1
        endif
    endif
else if(path0(ia,2).eq.i) then
    if(sinfo(ia).eq.1.and.info(path0(ia,1)).eq.0) then
        if(ip.eq.ip0.and.icpl.eq.7) then
            dd1=sum(cpl2(icpl,0:2)*pp(path0(ia,1),:))+cpl2(icpl,3)
            if(dd1.gt.0.0d0) then
             info(path0(ia,1))=-1
            else
             info(path0(ia,1))=1              
            endif
            continue
        else    
     info(path0(ia,1))=1
        endif
    endif
endif
enddo
enddo


!write(6,*) pathu(7,1:2),'path'
#ifdef BUGWIN
if(nn.eq.49031) then
continue
endif
#endif

!交点がないセル辺を登録
if(nff(mn(in(nn),kn(nn)-1,jn(nn)))%f.eq.-1) then
where(info(1:4).eq.0) info(1:4)=-1
!info(1:4)=-1
endif
2391 icc1=0
do i=1,12
if(sinfo(i).eq.0) cycle
if(info(path0(i,1)).eq.-1) then
!    if(info(path0(i,2)).ne.-1) then
    if(info(path0(i,2)).eq.0) then    
    info(path0(i,2))=-1
    icc1=icc1+1
    endif
else if(info(path0(i,2)).eq.-1) then
!   if(info(path0(i,1)).ne.-1) then
    if(info(path0(i,1)).eq.0) then    
    info(path0(i,1))=-1
    icc1=icc1+1
    endif   
endif
enddo

if(icc1.ne.0) goto 2391


if(sum(sinfo(1:12)).eq.12) then
if(ipathu_bk.eq.1) goto 2395    
    ipathu_bk=1
 !   icpl=icpl_bk
endif


do i=1,12
if(sinfo(i).eq.0) cycle    !交点がある点は除く


if(info(path0(i,1)).ne.-1.and.info(path0(i,2)).ne.-1) then
icc=0
do ia=1,ipathu
if(pathu(ia,1).eq.path0(i,1).and.pathu(ia,2).eq.path0(i,2)) icc=1 
if(pathu(ia,1).eq.path0(i,2).and.pathu(ia,2).eq.path0(i,1)) icc=1
enddo
if(icc.eq.0) then
    icc1=0
do ia=1,ipathu
    if(pathu(ia,-2).ne.0.and.pathu(ia,-1).eq.100+i) then
        icc1=1
    endif
enddo    
    if(icc1.eq.0) then
        ipathu=ipathu+1        
        pathu(ipathu,1:2)=path0(i,1:2)
                            pathu(ipathu,-1)=100+i
    endif
else
continue        
endif 
else 
        ipatho=ipatho+1        
        patho(ipatho,1:2)=path0(i,1:2)
                            patho(ipatho,-1)=100+i   
endif
enddo

#ifdef BUGWIN
if(nn.eq.59160) then
    continue
endif
#endif

!do i=1,ipathu
!write(6,*) i,pathu(i,1:2),'pathaaaaa'
!enddo
2395 continue
icc1=0
do i=1,ipathu
icc=0
do ib=1,icc1
if(ib.eq.i) cycle
if(pathu(i,1).eq.pathu(ib,1).and.pathu(i,2).eq.pathu(ib,2)) then
icc=ib
else if(pathu(i,1).eq.pathu(ib,2).and.pathu(i,2).eq.pathu(ib,1)) then
icc=ib
endif
enddo
if(icc.eq.0) then
icc1=icc1+1
    pathu(icc1,:)=pathu(i,:)
else
    if(pathu(i,-1).ne.pathu(icc,-1)) then
        if(pathu(i,-1).gt.100.and.pathu(icc,-1).lt.100) then
           pathu(icc,-1)=pathu(i,-1)
        else if(pathu(icc,-1).gt.100.and.pathu(i,-1).lt.100) then
            continue
        else
            continue
        endif
    endif
endif
enddo
ipathu=icc1

!do i=1,ipathu
!write(6,*) i,pathu(i,1:2),'pathbbbb'
!enddo

!
!continue
if(nn.eq.46557) then
    call checkdata3
continue
endif

!write(6,*) pathu(7,1:2),'path'
ncpl3=ncpl2
ncpl2(0,1)=ip
incpl2=ncpl2(0,0)

for_every_face :do i=1,ncpl2(0,0)
if(ncpl2(i,-1).eq.1) cycle

spath=0
do ia=1,ipathu
    if(pathu(ia,-1).eq.0) cycle
    if(pathu(ia,-1).eq.i) then   !pathuが格子面iに含まれる場合
        spath=spath+1
        paths(spath,-2:2)=pathu(ia,-2:2)
        paths(spath,3)=1
    else
        do ib=1,4
            if(pathu(ia,-1).eq.men(i,ib)+100) then   !pathuが格子面iの辺と一致する場合
                spath=spath+1
                paths(spath,-2:2)=pathu(ia,-2:2)
                paths(spath,3)=1
            endif
        enddo
    endif
enddo

!if(spath.le.1) then
!ncpl2(i,0)=0
!goto 444
!endif 

if(spath.le.2) then
    ncpl2(i,0)=0
    goto 444
endif

loop=0
443 icc1=0
icc3=0
do ia=1,ip
    icc=0
    do ib=1,spath
        if(paths(ib,3).eq.0)  cycle       
        if(paths(ib,1).eq.ia) icc=icc+1
        if(paths(ib,2).eq.ia) icc=icc+1
    enddo
    if(icc.gt.2) then
        icc1=icc1+1
        icv1(icc1)=ia
        icv2(icc1)=icc
    else if(icc.eq.1) then
        icc3=icc3+1
        icv3(icc3)=ia
    endif
enddo

if(icc3.gt.0) then  !ひげをとる。
    do ib=1,spath
        if(paths(ib,3).eq.0)  cycle       
        if(paths(ib,1).eq.icv3(1)) then
            paths(ib,3)=0
        else if(paths(ib,2).eq.icv3(1)) then
            paths(ib,3)=0
        endif
    enddo
    goto 443
endif

if(nn.eq.145584.and.i.eq.2) then
    continue
endif

if(icc1.eq.0.or.icc1.eq.1) then !　分岐がない場合 or 八の時
    ii=i;loop=0;loop1=0
icc=1
nnp1=0
ncpl20(-2:-1)=ncpl2(i,-2:-1)  
surf2=cpl2(i,:)
if(icc1.eq.0) then
icv(icc)=paths(1,1)
else
icv(icc)=icv1(icc1)
continue
endif
123 do ib=1,spath
    if(paths(ib,3).eq.0) cycle
    if(paths(ib,1).eq.icv(icc)) then
        icvm(icc)=paths(ib,-1)
        icvn(icc)=paths(ib,-2)
        icc=icc+1
        icv(icc)=paths(ib,2);paths(ib,3)=0;
        exit
    else if(paths(ib,2).eq.icv(icc)) then
        icvm(icc)=paths(ib,-1)
        icvn(icc)=paths(ib,-2)
        icc=icc+1
        icv(icc)=paths(ib,1);paths(ib,3)=0;
        exit
    endif   
    enddo
    if(icv(icc).ne.icv(1)) then
            loop=loop+1
      if(loop.gt.100) then

          if(icc1.eq.1) then
             loop=0
             icc=1
             icv(icc)=icv1(icc1)
             goto 123
          else 
             loop=0
             loop1=loop1+1
             if(loop1.gt.10) then
             write(6,*) 'loop1=',loop1
             continue
             stop
             endif
             icc=1
             icv(icc)=paths(1,1)
             goto 123            
          endif
                tmpV(0:2)=(/minval(pp(1:ncpl2(0,1),0)),minval(pp(1:ncpl2(0,1),1)),minval(pp(1:ncpl2(0,1),2))/)
tmpV2(0:2)=(/maxval(pp(1:ncpl2(0,1),0)),maxval(pp(1:ncpl2(0,1),1)),maxval(pp(1:ncpl2(0,1),2))/)
!write(63,*) 'nn=,',nn,',i,k,j=,',in(nn),',',kn(nn),',',jn(nn)
!write(63,*) vv  !,',',fb(nn)*dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn))
do ib=1,ncpl2(0,1)
write(63,661) ib,(pp(ib,0:2)-tmpV(0:2))/(tmpV2(0:2)-tmpV(0:2))
661 format(1h ,i2,3(',',e14.7))
enddo
do ib=1,ncpl2(0,0)
write(63,662) ib,ncpl2(ib,0:ncpl2(ib,0))
662 format(1h ,i2,8(',',i2))
enddo
        write(6,*) 'SegmentShape No2. icc1=0, loop=',loop,'nn=',nn,'i,k,j=',in(nn),kn(nn),jn(nn)
        continue
#ifdef BUGWIN              
              goto 932
#else
              goto 932
              stop
#endif 
      endif 
        goto 123
    else
        if(icc-1.eq.spath) then
            call surface_dir3(spath,icv(1:spath),cpl2(i,0:2),ncpl2(ii,1:spath))
            ncpl2(ii,0)=spath
        else
            nnp=icc-1
               call surface_dir3(nnp,icv(1:nnp),cpl2(ii,0:2),ncpl2(ii,1:nnp))            
            ib1=0
            do ib=1,nnp
            if(icvn(ib).gt.6) ib1=ib1+1   !  　20110826 ii -> i に変更
            enddo
            
            if(ib1.eq.nnp) then
                 if(icv(1).ne.ncpl2(ii,1)) then
                 icv(1:nnp)=ncpl2(ii,1:nnp)
                 icvnt(1:nnp)=icvn(1:nnp)
                 do ib3=1,nnp-1
                 icvn(ib3)=icvnt(nnp-ib3)
                 enddo
                 endif
                 dd1=10.0d0
                 do ib3=1,nnp 
                     ib3p=ib3+1;if(ib3.eq.nnp) ib3p=1
                     tmpV=crossV( pp(icv(ib3p),0:2)-pp(icv(ib3),0:2) , upl(idri_kdc(icvn(ib3)),idri_fc_n(icvn(ib3)),1:3) )
                     dd1=dot_product(surf2(0:2),tmpV)
               !      if(i.eq.1) dd1=dd1*(-1.0d0)
                     if(dd1.lt.0.0d0) exit
                 enddo
                 if(dd1.lt.0.0d0) then 
                    cpl2(ii,:)=surf2(:)*(-1.0d0)
                    call surface_dir3(nnp,icv(1:nnp),cpl2(ii,:),ncpl2(ii,1:nnp))
                    ncpl2(ii,-1)=1
                    ncpl2(ii,-2)=0   
                    ncpl2s(ii)=i  
                else    
                   cpl2(ii,0:3)=surf2                                   
                endif
            else 
               cpl2(ii,0:3)=surf2  
            endif
 !           call surface_dir(nnp,icv(1:nnp),cpl2(ii,0:2),ncpl2(ii,1:nnp))
            ncpl2(ii,0)=nnp
            nnp1=nnp1+nnp
            if(nnp1.lt.spath.and.spath-nnp1.gt.2.and.sum(paths(1:spath,3)).ne.0) then
                icc=1 
                do ib=1,spath
                    if(paths(ib,3).eq.1) then
                    icv(icc)=paths(ib,1)
                    exit
                    endif
                enddo
                    incpl2=incpl2+1
                    continue
                    cpl2(incpl2,:)=surf2  !cpl2(i,:) 
                    ncpl2(incpl2,-2:-1)=ncpl20(-2:-1)  !ncpl2(i,-2:-1)                   
                    ncpl2(0,0)=incpl2
                    ii=incpl2
                    continue
                    if(sum(paths(1:spath,3)).lt.3) then
                    continue
                   cycle 
                    endif
                    
                    
                 goto 123
            else
                 if(spath-nnp1.le.2.and.spath.ne.nnp1) then
                 continue
                 endif
            endif
        endif

    endif
            goto 444
endif


#ifdef DDDD
if(icc3.gt.0.and.icc1.eq.0) then
    icc4=0
    do ib1=1,icc3
        ia=icv3(ib1)
        icc=0
        do ib=1,spath
            if(paths(ib,1).eq.ia.or.paths(ib,2).eq.ia) then
                icc=ib
                exit
            endif
        enddo
        if(paths(icc,-1).eq.2) then
            paths(icc,3)=0
            icc4=icc4+1
        continue
        endif
    enddo
    if(icc4.ge.1.and.loop.lt.10) then
        loop=loop+1
        goto 443
    endif
    write(6,*) 'arimasen nn=',nn,'i,icc3,icc1=',i,icc3,icc1
    write(6,*) paths(1:spath,1)
    write(6,*) paths(1:spath,2)
    call checkdata3
    ncpl2(i,0)=0
    goto 444
endif
#endif

if(icc1.eq.2) then
    
call path_find_icc1_2(pp(1:ip,0:2),ip,cpl2(i,:),icc1,icv1(1:icc1),icv2(1:icc1),paths(1:spath,-2:3),spath,ncpl4,icpl,ireturn)
    surf2=cpl2(i,:)
    ncpl2(i,0)=ncpl4(1,0)
    ncpl2(i,1:ncpl4(1,0))=ncpl4(1,1:ncpl4(1,0))
    ncpl20(-2:-1)=ncpl2(i,-2:-1) 
    do i1=2,ncpl4(0,0)
        incpl2=incpl2+1
                    cpl2(incpl2,:)=surf2  !cpl2(i,:) 
                    ncpl2(incpl2,-2:-1)=ncpl20(-2:-1)  !ncpl2(i,-2:-1)                   
    ncpl2(incpl2,0)=ncpl4(i1,0)
    ncpl2(incpl2,1:ncpl4(i1,0))=ncpl4(i1,1:ncpl4(i1,0))
    ncpl2(0,0)=incpl2
        continue
    enddo
    continue
    if(sum(paths(1:spath,3)).lt.3) cycle
        icc=1
        do i1=1,spath
        if(paths(i1,3).eq.1) then
            icv(icc)=paths(i1,1)
            exit
        endif
        enddo
223 do ib=1,spath
    if(paths(ib,3).eq.0) cycle
    if(paths(ib,1).eq.icv(icc)) then
        icvm(icc)=paths(ib,-1)
        icvn(icc)=paths(ib,-2)
        icc=icc+1
        icv(icc)=paths(ib,2);paths(ib,3)=0;
        exit
    else if(paths(ib,2).eq.icv(icc)) then
        icvm(icc)=paths(ib,-1)
        icvn(icc)=paths(ib,-2)
        icc=icc+1
        icv(icc)=paths(ib,1);paths(ib,3)=0;
        exit
    endif   
    enddo
    if(icv(icc).ne.icv(1)) then
        goto 223
    else
            nnp=icc-1
            incpl2=incpl2+1
            cpl2(incpl2,:)=surf2  !cpl2(i,:) 
            ncpl2(incpl2,-2:-1)=ncpl20(-2:-1)  !ncpl2(i,-2:-1)  
               call surface_dir3(nnp,icv(1:nnp),cpl2(i,0:2),ncpl2(incpl2,1:nnp))            
            ib1=0
            do ib=1,nnp
            if(icvn(ib).gt.6) ib1=ib1+1   !  　20110826 ii -> i に変更
            enddo
            if(ib1.eq.nnp) then
                continue
            else
                ncpl2(incpl2,0)=nnp
                ncpl2(0,0)=incpl2
            endif
    endif
    !cpl2(ii,:)=surf2  !cpl2(i,:) 
!    ncpl2(ii,-2:-1)=ncpl20(-2:-1)  !ncpl2(i,-2:-1)                   
!    ncpl2(0,0)=incpl2

#ifdef DDDD
if(icv2(1).eq.3.and.icv2(2).eq.3) then
icc3=0
do ib=1,spath
    if(paths(ib,1).eq.icv1(1).and.paths(ib,2).eq.icv1(2)) icc3=1
    if(paths(ib,1).eq.icv1(2).and.paths(ib,2).eq.icv1(1)) icc3=1
enddo
if(icc3.eq.0) then
    iins=0
    do ib=1,spath
    if(paths(ib,1).eq.icv1(1)) then
        do ib1=1,spath
            if(ib.eq.ib1) cycle
            if(paths(ib1,1).eq.paths(ib,2).or.paths(ib1,2).eq.paths(ib,2)) then
                if(paths(ib1,1).eq.icv1(2).or.paths(ib1,2).eq.icv1(2)) then
                    continue
                endif
            endif
        enddo
    else if(paths(ib,2).eq.icv1(1)) then
        do ib1=1,spath
            if(ib.eq.ib1) cycle
            if(paths(ib1,1).eq.paths(ib,1).or.paths(ib1,2).eq.paths(ib,1)) then
                if(paths(ib1,1).eq.icv1(2)) then
                    iins=paths(ib1,2)  
                    exit          
                else if(paths(ib1,2).eq.icv1(2)) then
                    iins=paths(ib1,1)
                    exit
                    continue
                endif
            endif
        enddo
    endif
    enddo

    if(iins.eq.0) then
    write(6,*) 'iins=',iins
    goto 932
    stop
    endif

endif

surf2=cpl2(i,:)
ncpl20(-2:-1)=ncpl2(i,-2:-1)  
icc2=0
    ii=i;loop=0
125 icv(1)=icv1(1)
              !icv1(1)からicv1(2)への輪を探す．
129 icc=1
    icc2=icc2+1
124 do ib=1,spath
    if(paths(ib,3).eq.0) cycle
    if(paths(ib,1).eq.0) cycle    
    if(paths(ib,1).eq.icv1(1).and.paths(ib,2).eq.icv1(2)) then
    paths(ib,3)=0
    cycle
    endif
    if(paths(ib,1).eq.icv1(2).and.paths(ib,2).eq.icv1(1)) then
    paths(ib,3)=0    
    cycle    
    endif
    if(paths(ib,1).eq.icv(icc)) then    
        icc=icc+1
        icv(icc)=paths(ib,2);paths(ib,3)=0
        exit
    else if(paths(ib,2).eq.icv(icc)) then
        icc=icc+1
        icv(icc)=paths(ib,1);paths(ib,3)=0
        exit
    endif
enddo
    if(icc3.eq.0) then
        if(icc2.eq.1) then
            if(icv(icc).ne.icv1(1)) then
                    loop=loop+1
                if(loop.gt.100) then
                write(6,*) 'loop=1000a,nn=',nn
             goto 932 
                endif
                continue
                goto 124
             endif   
         else if(icc2.eq.2) then
             if(icv(icc).ne.icv1(2)) then
                    loop=loop+1
                if(loop.gt.100) then
                write(6,*) 'loop=1000b,nn=',nn
             goto 932 
                endif
                continue
                goto 124
             endif     
         endif      
    else
        if(icv(icc).eq.icv1(1)) then
            continue
        else if(icv(icc).ne.icv1(2)) then
            loop=loop+1
            if(loop.gt.100) then
                write(6,*) 'loop=1000c,nn=',nn
                call checkdata
                call checkdata3        
                goto 932 
            endif
            continue
            goto 124
        endif
     endif 
        if(icv(1).eq.icv(icc)) then
         nnp=icc-1       
        else
        nnp=icc
        endif
        ncpl2(ii,0)=nnp
        
        if(icc3.eq.0.and.icc2.eq.2) then
        nnp=nnp+1
        icv(nnp)=iins
        ncpl2(ii,0)=nnp
        continue
        endif
                
        call surface_dir3(nnp,icv(1:nnp),cpl2(ii,0:2),ncpl2(ii,1:nnp))  
        continue
        if(icc2.eq.1) then  !
        nnp1=nnp
        incpl2=incpl2+1
                    cpl2(incpl2,:)=surf2  !cpl2(i,:) 
                    ncpl2(incpl2,-2:-1)=ncpl20(-2:-1)  !ncpl2(i,-2:-1)                   
                    ncpl2(0,0)=incpl2
                    ii=incpl2        
                    do ib=1,spath
                     if(paths(ib,3).eq.1) then
                        icv(1)=paths(ib,1)
                        icv1(2)=paths(ib,1)
                        exit
                     endif
                    enddo 
        goto 129    
        
        else if(icc2.eq.2.and.sum(paths(1:spath,3)).ge.3) then
        incpl2=incpl2+1
                    cpl2(incpl2,:)=surf2  !cpl2(i,:) 
                    ncpl2(incpl2,-2:-1)=ncpl20(-2:-1)  !ncpl2(i,-2:-1)                   
                    ncpl2(0,0)=incpl2
                    ii=incpl2
                    do ib=1,spath
                     if(paths(ib,3).eq.1) then
                        icv(1)=paths(ib,1)
                        icv1(2)=paths(ib,1)
                        exit
                     endif
                    enddo 
        goto 129                    
        continue                
        endif 
!   endif 
   
   
else
continue
endif
#endif

else if(icc1.eq.4) then
    
call path_find_icc1_2(pp(1:ip,0:2),ip,cpl2(i,:),icc1,icv1(1:icc1),icv2(1:icc1),paths(1:spath,-2:3),spath,ncpl4,icpl,ireturn)
    surf2=cpl2(i,:)
    ncpl2(i,0)=ncpl4(1,0)
    ncpl2(i,1:ncpl4(1,0))=ncpl4(1,1:ncpl4(1,0))
    ncpl20(-2:-1)=ncpl2(i,-2:-1) 
    do i1=2,ncpl4(0,0)
        incpl2=incpl2+1
                    cpl2(incpl2,:)=surf2  !cpl2(i,:) 
                    ncpl2(incpl2,-2:-1)=ncpl20(-2:-1)  !ncpl2(i,-2:-1)                   
    ncpl2(incpl2,0)=ncpl4(i1,0)
    ncpl2(incpl2,1:ncpl4(i1,0))=ncpl4(i1,1:ncpl4(i1,0))
    ncpl2(0,0)=incpl2
        continue
    enddo
    
else 
    call checkdata3
    call checkdata
    write(6,*) 'icc1=',icc1,'nn=',nn
    goto 932 
continue
endif


!444 if(i.ge.1.and.i.le.6) then
444    if(i.eq.1.or.i.eq.3.or.i.eq.6) then
    call for_2d_vis  ! 二次元可視化用データ作成
endif


enddo for_every_face
ic3=0

icc1=0
do ib=1,ipathu
icc=0;icd=0
do ib1=1,ncpl2(0,0)
do ib2=1,ncpl2(ib1,0)
ib3=ib2+1;if(ib2.eq.ncpl2(ib1,0)) ib3=1
if(ncpl2(ib1,ib2).eq.pathu(ib,1).and.ncpl2(ib1,ib3).eq.pathu(ib,2)) then
    icd=icd+1
endif
if(ncpl2(ib1,ib2).eq.pathu(ib,2).and.ncpl2(ib1,ib3).eq.pathu(ib,1)) then
    icc=icc+1
endif
enddo
enddo
if(icc.eq.1.and.icd.eq.1) then
    continue
else if(icc+icd.le.2) then
write(6,224) ib,pathu(ib,1),pathu(ib,2),pathu(ib,-1),info(pathu(ib,1)),info(pathu(ib,2)),pathu(ib,3)
if(info(pathu(ib,1)).eq.-1.or.info(pathu(ib,2)).eq.-1)  ic3=1
continue
icc1=icc1+1
else if(icc+icd.gt.2) then
continue
else
continue
endif
enddo
224 format(1h ,'ib=',i3,',(p1,p2)=',i3,',',i3,',pu(-1)=',i3,',inf(p1,p2)=',i3,',',i3,',pu3=',i3)
if(icc1.eq.1) then
    write(6,*) 'ic3=1,icc1=1'
else if(icc1.ne.0) then    
  !  call checkdata3
  !  call checkdata
    write(6,*) 'pathu fu^{icchi} nn,icc1=',nn,icc1,'i,k,j=',in(nn),kn(nn),jn(nn)
    write(6,*) 'ubound(info,1)=',ubound(info,1)
 !   stop
    goto 932
endif

#ifdef DDDDDD
!!!!!!!!!!!!!!!!!!!　二分割の解消

allocate(pl_gr(1:ncpl2(0,0)))
pl_gr=0
i1=0
allocate(vt_gr(0:ncpl2(0,1)))
398 do ic3=1,ncpl2(0,0)
if(pl_gr(ic3).gt.0) cycle
if(ncpl2(ic3,0).ne.0) then
goto 397
endif
enddo !はじめの面をic3とする。

!はじめの面がない　終了。
if(i1.eq.0) then
deallocate(pl_gr,vt_gr)
goto 295
else
goto 396
endif
    
397 i1=i1+1 
vt_gr(1:ncpl2(ic3,0))=ncpl2(ic3,1:ncpl2(ic3,0)); !はじめの面をic3、two_grにはじめの面の頂点を登録v
vt_grn=ncpl2(ic3,0);   !現在のtwo_grnの頂点数
pl_gr(ic3)=i1  !　はじめの面 ic3をグループi1とする。
ic2=0;vt_grnp=vt_grn
222 do ic=ic3+1,ncpl2(0,0) 
    if(pl_gr(ic).eq.1) cycle  !すでに１に属しているものをのぞく
    icc=0
    do ia=1,vt_grn
        do ib=1,ncpl2(ic,0) 
            if(vt_gr(ia).eq.ncpl2(ic,ib)) then
                continue
                icc=icc+1
            endif
        enddo
    enddo
    if(icc.gt.1) then   !面ic の頂点が２点以上　two_gr(1:two_grn)　に含まれる場合　面icも　グループ１とする。
        pl_gr(ic)=i1
        icc=vt_grn
        do ib=1,ncpl2(ic,0)
        ic1=0
            do ia=1,vt_grn
            if(vt_gr(ia).eq.ncpl2(ic,ib)) ic1=1
            enddo
            if(ic1.eq.0) then
                icc=icc+1
                vt_gr(icc)=ncpl2(ic,ib)  !重複しないように　面icの点を　two_grに追加する。
            endif    
        enddo
        vt_grn=icc
    else
#ifdef BUGWIN2
        call checkdata3
        call checkdata
#endif        
        continue
    endif
    continue
enddo

if(vt_grnp.ne.vt_grn) then
vt_grnp=vt_grn
goto 222
endif

!ic2=ic2+1
!if(ic2.eq.1) goto 222
ic2=sum(pl_gr(1:ncpl2(0,0)))
if(i1.eq.1.and.ic2.eq.ncpl2(0,0)) then
deallocate(pl_gr,vt_gr)
goto 295
else 
    goto 398
endif

396 allocate(vv0(1:i1));VV0=0.0d0
!vv1=0.0d0;s=0.0d0  !VV01=0.0d0;VV02=0.0d0;VV03=0.0d0
Hyomen_menseki_loop2: do i=1,ncpl2(0,0)
    nnp=ncpl2(i,0)
    if(nnp.eq.0) cycle
    tmpV=0.0d0
    do ib=2,nnp-1
        ib1=ncpl2(i,ib)
        ib2=ncpl2(i,ib+1)
        tmpV=crossV(pp(ib1,0:2)-pp(ncpl2(i,1),0:2),pp(ib2,0:2)-pp(ncpl2(i,1),0:2))+tmpV
    enddo
    if(pl_gr(i).ne.0) VV0(pl_gr(i))=VV0(pl_gr(i))+sum(tmpV*pp(ncpl2(i,1),0:2))/6.0d0
enddo Hyomen_menseki_loop2

if(sum(vv0(1:i1)).gt.dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn)).and.ipathu_bk.ne.0) then
   ! ipathu=ipathu_bk
#ifdef BUGWIN
   call checkdata
#endif
   write(6,*) 'sum(vv0).gt.VV, yarinaoshi, nn=',nn,sum(vv0(1:i1))/(dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn)))
   deallocate(pl_gr,vt_gr,vv0)
    goto 4538
endif
    allocate(imax(1:i1))
        imax(1:i1)=MAXLOC(VV0(1:i1))
    Hyomen_menseki_loop31: do i=1,ncpl2(0,0)
        if(pl_gr(i).ne.imax(1)) then
        ncpl2(i,0)=0;s(i)=0.0d0
        endif
    enddo Hyomen_menseki_loop31
    VV=vv0(imax(1))

continue
   deallocate(pl_gr,vt_gr,imax,vv0)
!endif
!!!!!!!!!!!!!!!!!!!　二分割の解消
!endif
#endif


#ifdef DDDDD
do i=7,ncpl2(0,0)
    if(ncpl2(i,0).ne.ncpl3(i,0)) then
        continue
    else
      if(ncpl2(i,1).eq.ncpl3(i,1)) then
        do ia=2,ncpl2(i,0)
        if(ncpl2(i,ia).ne.ncpl3(i,ia)) then
        continue
        endif
        enddo
      else
         icc=0
        do ib=1,ncpl2(i,0)
      if(ncpl2(i,1).eq.ncpl3(i,ib)) then
         icc=ib
      endif     
        enddo 
        if(icc.eq.0) then
        continue
        else
        
        do ib=1,ncpl2(i,0)
        ib1=icc+ib-1
        if(ib1.gt.ncpl2(i,0)) ib1=ib1-ncpl2(i,0)
      if(ncpl2(i,ib).ne.ncpl3(i,ib1)) then
         continue
      endif     
        enddo         
        continue
        endif
      endif    
    endif
enddo
#endif
!
295 continue

do i=1,ip
do ia=1,ncpl2(0,0)
do ib=1,ncpl2(ia,0)
if(ncpl2(ia,ib).eq.i) goto 5551
enddo
enddo
do ia=1,ncpl5(0,0)
do ib=1,ncpl5(ia,0)
if(ncpl5(ia,ib).eq.i) goto 5551
enddo
enddo
pp(i,:)=-10000.d0
5551 continue
enddo

icc=0
if(allocated(info2)) deallocate(info2)
allocate(info2(1:ip))
do i=1,ip
    if(pp(i,1).eq.-10000.0d0) cycle 
    icc=icc+1
    pp(icc,:)=pp(i,:)
    info(icc)=info(i)
    info1(icc)=info1(i)
    info2(icc)=info2(i)
    do ia=1,ncpl2(0,0)
        do ib=1,ncpl2(ia,0)
        if(ncpl2(ia,ib).eq.i) ncpl2(ia,ib)=icc
        enddo
    enddo
    do ia=1,ncpl5(0,0)
        do ib=1,ncpl5(ia,0)
        if(ncpl5(ia,ib).eq.i) ncpl5(ia,ib)=icc
        enddo
    enddo
enddo

ncpl2(0,1)=icc



icc=0
do i=1,ncpl2(0,0)
if(ncpl2(i,0).eq.0) cycle
icc=icc+1
        ncpl2s(icc)=ncpl2s(i)
        cpl2(icc,:)=cpl2(i,:)
        nnp=ncpl2(i,0)
        ncpl2(icc,-2:nnp+nnp)=ncpl2(i,-2:nnp+nnp)
enddo

ncpl2(0,0)=icc

#ifdef BUGWIN
allocate(leng(1:ncpl2(0,0)))
#endif
!体積のチェック start
vv1=0.0d0;VV=0.0d0;s=0.0d0
Hyomen_menseki_loop: do i=1,ncpl2(0,0)
nnp=ncpl2(i,0)
if(nnp.eq.0) cycle
tmpV=0.0d0
do ib=2,nnp-1
ib1=ncpl2(i,ib)
ib2=ncpl2(i,ib+1)
tmpV=crossV(pp(ib1,0:2)-pp(ncpl2(i,1),0:2),pp(ib2,0:2)-pp(ncpl2(i,1),0:2))+tmpV
continue
enddo
s(i)=sqrt(sum(tmpV**2))*HALF
continue
#ifdef BUGWIN    
    leng(i)=sum(tmpV*pp(ncpl2(i,1),0:2))/6.0d0
#endif   
VV=VV+sum(tmpV*pp(ncpl2(i,1),0:2))/6.0d0
enddo Hyomen_menseki_loop

if(abs(vv1-vv).gt.1.0d-8) then
continue
endif

#ifdef BUGWIN
deallocate(leng)
#endif


if(flag.ne.0) then
    if(vv/(dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn))).lt.1.0d-5) then
#ifdef BUGWIN
    write(6,*) nn,vv/(dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn))),in(nn),kn(nn),jn(nn)
#endif
#ifdef DRI
        mm11=mm11+1
        drino2i(mm11)=nn
        drino2(nn,1)=mm11;drino2(nn,2)=kdc    
    continue
#endif
        return
    endif
endif
!ic2p=ubound(pp_0,2)
!mm0mxp=ubound(OBJNOI,1)
!if(ncpl2(0,1).gt.ic2p) then 
!    write(6,*) 'ncpl2(0,1)=',ncpl2(0,1)
!endif

mm0_p=mm0_p+1


!OBJNO(nn)=mm0   
ncpl2(0,2)=2*maxval(ncpl2(1:ncpl2(0,0),0))  !一面あたりの最大頂点数
pp_p(1:ncpl2(0,1),0:2)=pp(1:ncpl2(0,1),0:2)
info_p(1,1:ncpl2(0,1))=info1(1:ncpl2(0,1))
info_p(2,1:ncpl2(0,1))=info2(1:ncpl2(0,1))
#ifdef TEST
OBJNOI(mm0)=nn
ncpl0(mm0,:,:)=0
ncpl0(mm0,0:ncpl2(0,0),-2:ncpl2(0,2))=ncpl2(0:ncpl2(0,0),-2:ncpl2(0,2))
pp_0(mm0,1:ncpl2(0,1),0:2)=pp(1:ncpl2(0,1),0:2)
info_0(mm0,1,1:ncpl2(0,1))=info1(1:ncpl2(0,1))
info_0(mm0,2,1:ncpl2(0,1))=info2(1:ncpl2(0,1))
#endif TEST

if(flag.eq.0) then
!if(mm0.eq.1) then
!write(6,*) ncpl0(mm0,1:2,-2:-1)
!call checkdata
!stop
!endif


#ifdef DRI
else
!if(nff(nn)%f.eq.-1) then
! nfd(in(nn),kn(nn),jn(nn))=1
! nfbd(in(nn),kn(nn),jn(nn))=0
!endif
DRINO(nn)=mm0
DRINOI(mm0)=nn
ncpl1(mm0,:,:)=ncpl2(:,:)
pp_1(mm0,1:ncpl2(0,1),0:2)=pp(1:ncpl2(0,1),0:2)


ib1=0
idvd=1
info(1:ncpl2(0,1))=1
info2(1:ncpl2(0,0))=0
do i=1,ncpl2(0,0)
if(ncpl2(i,0).eq.0) cycle
info2(i)=idvd;icc=1
    do ib=1,ncpl2(i,0)
    info(ncpl2(i,ib))=0
    enddo
exit
enddo
598 do i=1,ncpl2(0,0)
    if(ncpl2(i,0).eq.0) cycle
    if(info2(i).ne.0) cycle
    do ib=1,ncpl2(i,0)
        if(info(ncpl2(i,ib)).eq.0) goto 599
    enddo
    cycle
    599 info2(i)=idvd
        icc=icc+1
    do ib=1,ncpl2(i,0)
        info(ncpl2(i,ib))=0
    enddo
enddo
if(icc.lt.ncpl2(0,0)) then
ib1=ib1+1
if(ib1.gt.5) then
    continue
    idvd=idvd+1
     do i=1,ncpl2(0,0)
        if(info2(i).ne.0) cycle
        if(ncpl2(i,0).eq.0) cycle
        info2(i)=idvd;icc=icc+1
        do ib=1,ncpl2(i,0)
        info(ncpl2(i,ib))=0
        enddo
        exit
    enddo
    goto 598
endif
goto 598
endif

#ifdef DRI2
if(idvd.ge.2) then
mm2p=mm2p+1
DRI2PNO(mm2p,1)=nn
DRI2PNO(mm2p,2)=idvd
DRI2PNOI(nn)=mm2p

do ib3=1,idvd
axyz(mm2p,ib3,1:6)=0.0d0
axyzd(mm2p,ib3,1:6)=0.0d0
vv2=0.0d0

Hyomen_menseki_loop2: do i=1,ncpl2(0,0)
    nnp=ncpl2(i,0)
    if(info2(i).ne.ib3) cycle
    if(nnp.eq.0) cycle
    tmpV=0.0d0
    do ib=2,nnp-1
    ib1=ncpl2(i,ib)
    ib2=ncpl2(i,ib+1)
    tmpV=crossV(pp(ib1,0:2)-pp(ncpl2(i,1),0:2),pp(ib2,0:2)-pp(ncpl2(i,1),0:2))+tmpV
    continue
    enddo
    s2(i)=sqrt(sum(tmpV**2))*HALF

    if(ncpl2(i,-1).eq.0) then
    if(ncpl2(i,-2).eq.1) then
        if(mn(in(nn),kn(nn)-1,jn(nn)).ge.inns.and.mn(in(nn),kn(nn)-1,jn(nn)).le.inne) then
            axyzd(mm2p,ib3,1)=axyz(mm2p,ib3,1)+s2(i)/(dx(0,in(nn))*dx(1,jn(nn)))
        endif
    else if(ncpl2(i,-2).eq.2) then
        if(mn(in(nn),kn(nn)+1,jn(nn)).ge.inns.and.mn(in(nn),kn(nn)+1,jn(nn)).le.inne) then
            axyzd(mm2p,ib3,2)=axyz(mm2p,ib3,2)+s2(i)/(dx(0,in(nn))*dx(1,jn(nn)))
        endif
    else if(ncpl2(i,-2).eq.3) then
        if(mn(in(nn),kn(nn),jn(nn)-1).ge.inns.and.mn(in(nn),kn(nn),jn(nn)-1).le.inne) then    
            axyzd(mm2p,ib3,3)=axyz(mm2p,ib3,3)+s2(i)/(dx(2,kn(nn))*dx(0,in(nn)))
        endif
    else if(ncpl2(i,-2).eq.4) then
        if(mn(in(nn)+1,kn(nn),jn(nn)).ge.inns.and.mn(in(nn)+1,kn(nn),jn(nn)).le.inne) then      
            axyzd(mm2p,ib3,4)=axyz(mm2p,ib3,4)+s2(i)/(dx(2,kn(nn))*dx(1,jn(nn)))
        endif
    else if(ncpl2(i,-2).eq.5) then
        if(mn(in(nn),kn(nn),jn(nn)+1).ge.inns.and.mn(in(nn),kn(nn),jn(nn)+1).le.inne) then    
            axyzd(mm2p,ib3,5)=axyz(mm2p,ib3,5)+s2(i)/(dx(2,kn(nn))*dx(0,in(nn)))
        endif
    else if(ncpl2(i,-2).eq.6) then
        if(mn(in(nn)-1,kn(nn),jn(nn)).ge.inns.and.mn(in(nn)-1,kn(nn),jn(nn)).le.inne) then 
            axyzd(mm2p,ib3,6)=axyz(mm2p,ib3,6)+s2(i)/(dx(2,kn(nn))*dx(1,jn(nn)))
        endif
    endif
    endif

    continue
 
    VV2=VV2+sum(tmpV*pp(ncpl2(i,1),0:2))/6.0d0
    continue
enddo Hyomen_menseki_loop2

FB2PD(mm2p,ib3)=VV2/(dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn)))
continue
enddo

if(DRI2PNOIP(nn).ne.0) then
continue


else
if(abs(FB2PD(mm2p,1)-fb(nn)).lt.abs(FB2PD(mm2p,2)-fb(nn)) ) then
F2P(mm2p,1)=fp(nn);F2P(mm2p,2)=0.0d0
FB2P(mm2p,1)=fb(nn);FB2P(mm2p,2)=0.0d0
P2P(mm2p,1)=p(nn);P2P(mm2p,2)=0.0d0;
nf2P(mm2p,1)=1;nf2P(mm2p,2)=0;
    axyz(mm2p,1,1)=a(2,nn)
    axyz(mm2p,1,2)=a(2,mn(in(nn),kn(nn)+1,jn(nn)))
    axyz(mm2p,1,3)=a(1,nn)
    axyz(mm2p,1,4)=a(0,mn(in(nn)+1,kn(nn),jn(nn)))
    axyz(mm2p,1,5)=a(1,mn(in(nn),kn(nn),jn(nn)+1))
    axyz(mm2p,1,6)=a(0,nn)
    axyz(mm2p,2,1:6)=0.0d0    
    
    Vxyz(mm2p,1,1)=w(nn)
    Vxyz(mm2p,1,2)=w(mn(in(nn),kn(nn)+1,jn(nn)))
    Vxyz(mm2p,1,3)=v(nn)
    Vxyz(mm2p,1,4)=u(mn(in(nn)+1,kn(nn),jn(nn)))
    Vxyz(mm2p,1,5)=v(mn(in(nn),kn(nn),jn(nn)+1))
    Vxyz(mm2p,1,6)=u(nn)
    Vxyz(mm2p,2,1:6)=0.0d0    

    axyz2=(axyzd+axyz)*HALF        

continue
else
F2P(mm2p,2)=fp(nn);F2P(mm2p,1)=0.0d0
FB2P(mm2p,2)=fp(nn);FB2P(mm2p,1)=0.0d0
P2P(mm2p,2)=p(nn);P2P(mm2p,1)=0.0d0;
    axyz(mm2p,2,1)=a(2,nn)
    axyz(mm2p,2,2)=a(2,mn(in(nn),kn(nn)+1,jn(nn)))
    axyz(mm2p,2,3)=a(1,nn)
    axyz(mm2p,2,4)=a(0,mn(in(nn)+1,kn(nn),jn(nn)))
    axyz(mm2p,2,5)=a(1,mn(in(nn),kn(nn),jn(nn)+1))
    axyz(mm2p,2,6)=a(0,nn)
    axyz(mm2p,1,1:6)=0.0d0    
    
    Vxyz(mm2p,2,1)=w(nn)
    Vxyz(mm2p,2,2)=w(mn(in(nn),kn(nn)+1,jn(nn)))
    Vxyz(mm2p,2,3)=v(nn)
    Vxyz(mm2p,2,4)=u(mn(in(nn)+1,kn(nn),jn(nn)))
    Vxyz(mm2p,2,5)=v(mn(in(nn),kn(nn),jn(nn)+1))
    Vxyz(mm2p,2,6)=u(nn)
    Vxyz(mm2p,1,1:6)=0.0d0    

    axyz2=(axyzd+axyz)*HALF
endif

endif

else
DRI2PNOI(nn)=0
endif
#endif
#endif
endif

#ifdef BUGWIN
!if(nn.eq.17415) then
!call checkdata
!call checkdata3
!continue
!endif
#endif
if(flag.eq.0) then
    axyz(0)=vv/(dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn)))

    if(axyz(0).gt.1.0d0+1.0d-10.or.axyz(0).lt.0.0d0) then
    if(ipathu_bk.eq.1) then
           write(6,*) 'fb0(nn).gt.1.0d0, yarinaoshi, nn=',nn
       !     icpl=icpl_bk
        goto 4538
    endif
    write(6,*) 'nn=',nn,'fb0=',axyz(0),in(nn),jn(nn),kn(nn)

    do i=1,ncpl2(0,0)
    write(62,664) i,s(i),ncpl2(i,-2:ncpl2(i,0))
    !write(6,664) i,s(i),ncpl2(i,-2:3)
    664 format(1h ,i2,(',',e14.7),20(',',i6))
    enddo
    call checkdata
    !stop
    continue
    endif

else
    !fbd(nn)=vv/(dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn)))
    axyz(0)=vv/(dx(0,in(nn))*dx(1,jn(nn))*dx(2,kn(nn)))
endif
tmpp(1:6)=0.0d0;icc=0
kaiko_rate_hyoumen_info_loop0:do i=1,ncpl2(0,0)
if(ncpl2(i,-1).eq.0) then    !物体表面の場合
tmpp(ncpl2(i,-2))=tmpp(ncpl2(i,-2))+s(i)
continue
else if(ncpl2(i,-1).eq.1.and.ncpl2(i,-2).eq.0.and.ncpl2s(i).ne.0) then !漂流物表面の場合，
tmpp(ncpl2s(i))=tmpp(ncpl2s(i))-s(i)
endif
enddo kaiko_rate_hyoumen_info_loop0
#ifdef DRIOLD
if(axyzfb1(mm0,0).gt.1.0d-3) then
axyzfb1(mm0,1)=tmpp(1)/dx(1,jn(nn))/dx(0,in(nn))
axyzfb1(mm0,2)=tmpp(2)/dx(1,jn(nn))/dx(0,in(nn))
axyzfb1(mm0,3)=tmpp(3)/dx(0,in(nn))/dx(2,kn(nn))
axyzfb1(mm0,4)=tmpp(4)/dx(1,jn(nn))/dx(2,kn(nn))
axyzfb1(mm0,5)=tmpp(5)/dx(0,in(nn))/dx(2,kn(nn))
axyzfb1(mm0,6)=tmpp(6)/dx(1,jn(nn))/dx(2,kn(nn))
else
        mm11=mm11+1
        drino2i(mm11)=nn
        drino2(nn,1)=mm11;drino2(nn,2)=kdc   
        drino(nn)=0
        mm0=mm0-1
        cycle
endif 
#endif
if(flag.ge.0) then
axyz(1)=tmpp(1)/dx(1,jn(nn))/dx(0,in(nn))
axyz(2)=tmpp(2)/dx(1,jn(nn))/dx(0,in(nn))
axyz(3)=tmpp(3)/dx(0,in(nn))/dx(2,kn(nn))
axyz(4)=tmpp(4)/dx(1,jn(nn))/dx(2,kn(nn))
axyz(5)=tmpp(5)/dx(0,in(nn))/dx(2,kn(nn))
axyz(6)=tmpp(6)/dx(1,jn(nn))/dx(2,kn(nn))
endif

#ifdef BUGWIN2
if(nn.eq.4550543) then
call checkdata
call checkdata3
stop
endif
#endif
return

932 continue
#ifdef BUGWIN
call checkdata
call checkdata3
! write(6,*) 'nn=',nn,fb0(nn),vv,objno(nn),mm0
stop
#endif

    contains

subroutine checkdata3
integer::ib
write(612,*) nn,ip   !ncpl2(0,1)
do ib=1, ip   !ncpl2(0,1)
write(612,13) pp(ib,:)
13 format(1h ,e12.5,10(' ',e12.5))
continue
enddo
write(612,*) ipathu
do ib=1,ipathu
write(612,14) pathu(ib,1),pathu(ib,2)
14 format(1h ,i5,' ',i5)
enddo
write(612,*) ipatho
do ib=1,ipatho
write(612,14) patho(ib,1),patho(ib,2)
enddo
write(612,*) 'nn,i,k,j=',nn,in(nn),kn(nn),jn(nn)
end subroutine
subroutine checkdata
implicit none
integer::i,ib
write(611,*) ncpl2(0,1),ncpl2(0,1) !,objno(nn)
do ib=1, ncpl2(0,1)
write(611,11) pp(ib,:)
11 format(1h ,e12.5,10(' ',e12.5))
continue
enddo
write(611,*) ncpl2(0,0)
do ib=1, ncpl2(0,0)
write(611,12) ncpl2(ib,-2:ncpl2(ib,0))
12 format(1h ,i7,12(' ',i3))
enddo
!#ifdef BUGWIN
do ib=1, ncpl2(0,1)
do i=1,ivt(1)
if(dista(vt(1,i,:)-pp(ib,:)).lt.1.0d-10) then
write(611,*) ib,i
exit
endif
enddo
enddo
!#endif
write(611,*) 'nn=',nn

end subroutine
subroutine checkdata2
write(611,*) nn,ip
do ib=1, ip
write(611,11) pp(ib,:)
11 format(1h ,e12.5,10(' ',e12.5))
continue
enddo
write(611,*) icpl
do ib=1, icpl
write(611,12) ncpl2(ib,-2:ncpl2(ib,0))
12 format(1h ,i6,12(' ',i2))
enddo
end subroutine
subroutine checkdata2ia
integer::icvv(1:ipoia)
write(619,*) nn,ip+ipoia
do ib=1, ip
write(619,11) pp(ib,:)
enddo
do ib=1, ipoia
write(619,11) vtia(ib,:)
enddo

11 format(1h ,e12.5,10(' ',e12.5))
continue
ib=2+6
write(619,*) ib
do ib=1, 6
write(619,12) ncpl2(ib,-2:ncpl2(ib,0))
enddo
ib=1
write(619,12) ib,ib,icc00,(icv22(i),i=1,icc00)
do ib=1, ipoia
    if(icv4(ib).ne.0) then
        icvv(ib)=icv4(ib)
    else
        icvv(ib)=ib+ip
    endif
enddo
ib=2
write(619,12) ib,ib,ipoia,(icvv(i),i=1,ipoia)
!write(611,12) ib,ib,
12 format(1h ,i6,12(' ',i2))
!enddo
write(619,*) 'nn=',nn,'ia=',ia
end subroutine

subroutine cal_pathu(icc,ic2)
implicit none
integer,intent(in)::icc
integer,intent(out)::ic2
integer::path(minfopath,1:3)
integer::ic1,ipath,ib,ic
ipath=0
do ib=1,info3(icc,0)
do ic=1,ncpl2(info3(icc,ib),0)
ic1=ic+1; if(ic1.gt.ncpl2(info3(icc,ib),0)) ic1=1
if(ncpl2(info3(icc,ib),ic).eq.icc) then
ipath=ipath+1
path(ipath,1)=ncpl2(info3(icc,ib),ic1)
path(ipath,3)=info3(icc,ib)
continue
else if(ncpl2(info3(icc,ib),ic1).eq.icc) then
ipath=ipath+1
path(ipath,1)=ncpl2(info3(icc,ib),ic)
path(ipath,3)=info3(icc,ib)
continue
endif
enddo
enddo
ic2=0
do ic=1,ipath
ib=0
do ic1=1,ipath
if(ic.eq.ic1) cycle
if(path(ic,1).eq.path(ic1,1))  ib=1
enddo
if(ib.eq.0) then
ic2=ic2+1
path2(ic2,1:3)=path(ic,1:3)
!info3(icc,ic2)=path(ic,3)  !面番号
continue
endif
enddo
end subroutine cal_pathu

subroutine cal_pathu4(icc,p02)
implicit none
integer,intent(in)::p02,icc
            ipathu=ipathu+1
            pathu(ipathu,1)=p02  !path0(i,2)
            pathu(ipathu,-1)=100+i
            if(p02.le.8.and.info(p02).eq.0) info(p02)=1;
            pathu(ipathu,2)=icc
end subroutine cal_pathu4

subroutine cal_pathu3(icc,imen,p02)
implicit none
integer,intent(in)::p02,imen,icc
doubleprecision::dd1
        dd1=sum(cpl2(imen,0:2)*pp(p02,0:2))+cpl2(imen,3) 
        if(dd1.lt.0.0d0) then
            ipathu=ipathu+1
            pathu(ipathu,1)=p02  !path0(i,2)
            pathu(ipathu,-1)=100+i
            if(p02.le.8.and.info(p02).eq.0) info(p02)=1;
            pathu(ipathu,2)=icc
        else
            ipatho=ipatho+1
            patho(ipatho,1)=p02  !path0(i,2)  
            if(p02.le.8.and.info(p02).eq.0) info(p02)=-1;
            patho(ipatho,2)=icc         
            patho(ipatho,-1)=100+i  
        endif
end subroutine cal_pathu3        

subroutine cal_pathu2(icc,pc,pc2)   !対象点、
integer,intent(in)::icc,pc,pc2
integer::ic2,imen !,it1,it2,ii13,ii2,ii23,ii1
doubleprecision::tt1,tt2  !,tt3,tt4,dd1
!write(6,*) info3(icc,0)
    if(info3(icc,0).eq.1) then !交点iccが一つの平面と交わってるとき、端点pcとの関係を調べる。
        imen=info3(icc,1)
        call cal_pathu3(icc,imen,pc)
    else if(info3(icc,0).ge.2) then
        call cal_pathu(icc,ic2) !交点iccが二つの平面と交わってるとき、両端の点をpath2(1,1:ic2)に入力
        if(ic2.eq.2) then
        else if(ic2.eq.4) then
     !       write(6,*) path2(1:4,1)
      !      write(6,*) path2(1:4,3)            
      !      tt1=sum( (pp(path2(1,1),0:2)-pp(pc,0:2))**2 ) 
      !      tt2=sum( (pp(path2(2,1),0:2)-pp(pc,0:2))**2 )
      !      tt3=sum( (pp(path2(3,1),0:2)-pp(pc,0:2))**2 ) 
      !      tt4=sum( (pp(path2(4,1),0:2)-pp(pc,0:2))**2 )            
        !write(6,*) 'ic2 gt.4 ic2=',ic2,' nn=',nn
       !    if(tt1.le.tt2.and.tt1.le.tt3.and.tt1.le.tt4) then
       !        ii1=path2(1,1);ii13=path2(1,3);tt1=10000.
       !    else if(tt2.le.tt1.and.tt2.le.tt3.and.tt2.le.tt4) then
       !        ii1=path2(2,1);ii13=path2(2,3);tt2=10000.            
       !    else if(tt3.le.tt1.and.tt3.le.tt2.and.tt3.le.tt4) then
       !        ii1=path2(3,1);ii13=path2(3,3);tt3=10000.  
       !    else if(tt4.le.tt1.and.tt4.le.tt2.and.tt4.le.tt3) then
       !        ii1=path2(4,1);ii13=path2(4,3);tt4=10000.   
       !    endif
       !    if(tt1.le.tt2.and.tt1.le.tt3.and.tt1.le.tt4) then
      !         ii2=path2(1,1);ii23=path2(1,3);
      !     else if(tt2.le.tt1.and.tt2.le.tt3.and.tt2.le.tt4) then
      !         ii2=path2(2,1);ii23=path2(2,3);        
      !     else if(tt3.le.tt1.and.tt3.le.tt2.and.tt3.le.tt4) then
       !        ii2=path2(3,1);ii23=path2(3,3);
     !      else if(tt4.le.tt1.and.tt4.le.tt2.and.tt4.le.tt3) then
     !          ii2=path2(4,1);ii23=path2(4,3); 
     !      endif
      !     if(ii13.eq.ii23) then
      !     imen=ii13
      !      call cal_pathu3(icc,imen,pc)
      !  else
       !         path2(1,1)=ii1
         !       path2(1,3)=ii13
        !        path2(2,1)=ii2
         !       path2(2,3)=ii23
                path2(2,1)=path2(4,1)
                path2(2,3)=path2(4,3)
  !      endif
           

!       info3(icc,2)=info3(icc,4)
        continue
        else
            write(6,*) 'ic2 gt.4 ic2=',ic2,' nn=',nn
            continue
        endif        
            if(path2(2,1).eq.pc.or.path2(1,1).eq.pc) then
                call cal_pathu4(icc,pc)
            else
            tt1=sum( (pp(path2(1,1),0:2)-pp(pc,0:2))*(pp(pc2,0:2)-pp(pc,0:2)) ) 
            tt2=sum( (pp(path2(2,1),0:2)-pp(pc,0:2))*(pp(pc2,0:2)-pp(pc,0:2)) )
            if(tt1.lt.tt2) then            
            imen=path2(1,3) !info3(icc,1)  
            else
            imen=path2(2,3) !info3(icc,2)              
            endif    
                call cal_pathu3(icc,imen,pc)
            endif
    else    
        continue
    endif
end subroutine cal_pathu2    

subroutine pathu_patuo(pcc,ip1,ip2)   !(path0(i,1),ip1)
implicit none
integer,intent(in)::pcc,ip1,ip2
if(info3(ip1,0).eq.1) then
!    call cal_pathu2(ip1,pcc)
    continue
else if(info3(ip1,0).eq.2) then
    continue
else
    continue
endif
    continue
end subroutine

subroutine for_2d_vis
implicit none
integer::kdc,tpath,icc3,icc,ia,ib,tpathm
integer,dimension(1:minfopath)::icv
integer,dimension(0:minfopath,-2:4)::patht

kdc=0;tpath=0
do ia=1,ipatho

if(patho(ia,-1).eq.0) cycle
if(patho(ia,-1).eq.i) then
    tpath=tpath+1
    patht(tpath,-2:3)=patho(ia,-2:3)
    patht(tpath,4)=1
        if(patho(ia,3).ne.0.and.kdc.eq.0) then
        kdc=patho(ia,3)
        endif  
else
    do ib=1,4
    if(patho(ia,-1).eq.men(i,ib)+100) then
        tpath=tpath+1
        patht(tpath,-2:3)=patho(ia,-2:3)
        patht(tpath,4)=1
            if(patho(ia,3).ne.0.and.kdc.eq.0) then
            kdc=patho(ia,3)
            endif  
    endif
    enddo
endif
enddo

if(kdc.eq.0) then
return !cycle
endif
;tpathm=0
4439 icc1=0
icc3=0
do ia=1,ip
    icc=0
    do ib=1,tpath
        if(patht(ib,4).eq.0)  cycle       
        if(patht(ib,1).eq.ia) icc=icc+1
        if(patht(ib,2).eq.ia) icc=icc+1
    enddo
    if(icc.eq.1) then
        icc3=icc3+1
        icv3(icc3)=ia
        exit
    endif
enddo

if(icc3.gt.0) then  !ひげをとる。
    
    do ib=1,tpath
        if(patht(ib,4).eq.0)  cycle       
        if(patht(ib,1).eq.icv3(1)) then
            patht(ib,4)=0
            tpathm=tpathm+1
        else if(patht(ib,2).eq.icv3(1)) then
            patht(ib,4)=0
            tpathm=tpathm+1            
        endif
    enddo
    goto 4439
endif

if(tpath-tpathm.le.1) then
ncpl3(i,0)=0
return !cycle
else if(tpath-tpathm.le.2) then
return !cycle
endif

do ia=2,tpath
    if(patht(ia,-1).ne.patht(1,-1)) then
        goto 9124
    endif
enddo

return
loop2=0
9124 do ia=1,tpath
    do ib=1,2
        ic1=0
        do ib1=1,tpath
            if(patht(ib1,4).eq.0) cycle        
            do ib2=1,2
                if(ia-ib1.eq.0.and.ib-ib2.eq.0) cycle

                if(patht(ia,ib).eq.patht(ib1,ib2)) ic1=1
                continue
            enddo
        enddo
        if(ic1.eq.0) then
            patht(ia,4)=0
        endif
    enddo
enddo

continue

do kdc=1,ndri
if(rhod(kdc).lt.0.0d0.and.flag.eq.1.and.n.gt.0) cycle
ii=i;loop=0
icc=1
nnp1=0
icv(icc)=0
do i1=1,tpath
if(patht(i1,3).eq.kdc) then
icv(icc)=patht(i1,1)
continue
exit
endif
enddo
if(icv(icc).eq.0) then 
continue
cycle
endif

9123 do ib=1,tpath
    if(patht(ib,4).eq.0) cycle
    if(patht(ib,1).eq.icv(icc)) then
        icvm(icc)=patht(ib,-1);        icvn(icc)=paths(ib,-2)
        icc=icc+1
        icv(icc)=patht(ib,2);patht(ib,4)=0;
        exit
    else if(patht(ib,2).eq.icv(icc)) then
        icvm(icc)=patht(ib,-1);        icvn(icc)=paths(ib,-2)
        icc=icc+1
        icv(icc)=patht(ib,1);patht(ib,4)=0;
        exit
    endif
    enddo
    
    
    if(icv(icc).ne.icv(1)) then
      loop=loop+1
      if(loop.gt.100) then
         loop2=loop2+1
         !     call checkdata3
              continue
              ib1=0;ib2=0
              patht(1:tpath,4)=1
              do ic=1,ip
              icc=0
              do ib=1,tpath
              if(patht(ib,1).eq.ic.or.patht(ib,2).eq.ic) icc=icc+1
              enddo 
              
              if(icc.eq.3) then
              if(ib1.eq.0) then
              ib1=ic
              else
              ib2=ic
              endif
              endif
              enddo 
              if(ib2.eq.0) then
              write(6,*) 'SegmentShape patht. icc1=0, loop=',loop,'nn=',nn,'i=',i,'tpath=',tpath
              call checkdata3
              continue
              cycle
              endif
              do ib=1,tpath
              if(patht(ib,1).eq.ib1.and.patht(ib,2).eq.ib2) patht(ib,4)=0
              if(patht(ib,1).eq.ib2.and.patht(ib,2).eq.ib1) patht(ib,4)=0               
              enddo 
              continue
        if(loop2.gt.100) then
              write(6,*) 'SegmentShape patht. icc1=0, loop=',loop,'nn=',nn,'i=',i,'tpath=',tpath  
             call checkdata3
             call checkdata
             return
        else
            
              goto 9124             
        endif
#ifdef BUGWIN              
  !            goto 932
#else
       
  !            goto 932
  !            stop
#endif 
              cycle        
      endif
     goto 9123
    endif

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!start!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!二次元可視化用データ取得!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(icc-1.eq.4) then
if(icv(1).gt.8) goto 230
if(icv(2).gt.8) goto 230
if(icv(3).gt.8) goto 230
if(icv(4).gt.8) goto 230
cycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif       

230 if(icc.eq.1) cycle

!230 mmd2(kdc)=mmd2(kdc)+1
    mmd2_p(kdc)=mmd2_p(kdc)+1    
#ifdef TEST                 
if(mmd2(kdc).gt.mmd2mx) then
    mmd2mxp=mmd2mx
    mmd2mx=mmd2mx+1000
    allocate(dummyi3(1:ndri,1:mmd2mxp,1:2))
    dummyi3=IDP2
    deallocate(IDP2)
    allocate(IDP2(1:ndri,1:mmd2mx,1:2))
    IDP2(1:ndri,1:mmd2mxp,1:2)=dummyi3(1:ndri,1:mmd2mxp,1:2)
    allocate(dummyi2(1:ndri,1:mmd2mxp))
    dummyi2=nd2
    deallocate(nd2)
    allocate(nd2(1:ndri,1:mmd2mx))
    nd2(1:ndri,1:mmd2mxp)=dummyi2(1:ndri,1:mmd2mxp)
    allocate(dummyr4(1:ndri,1:mmd2mxp,1:mpx2,1:3))
    dummyr4=dp2
    deallocate(dp2)
    allocate(dp2(1:ndri,1:mmd2mx,1:mpx2,1:3))
    dp2(1:ndri,1:mmd2mxp,1:mpx2,1:3)=dummyr4(1:ndri,1:mmd2mxp,1:mpx2,1:3)
    deallocate(dummyi2,dummyi3,dummyr4)
endif              
              select case(i)
                case(1);IDP2(kdc,mmd2(kdc),1)=3
                case(3);IDP2(kdc,mmd2(kdc),1)=2 
                case(6);IDP2(kdc,mmd2(kdc),1)=1  
                case(2);IDP2(kdc,mmd2(kdc),1)=3
                case(5);IDP2(kdc,mmd2(kdc),1)=2 
                case(4);IDP2(kdc,mmd2(kdc),1)=1 
              end select
              select case(i)
                case(3);IDP2(kdc,mmd2(kdc),2)=jn(nn)
                case(6);IDP2(kdc,mmd2(kdc),2)=in(nn)
                case(1);IDP2(kdc,mmd2(kdc),2)=kn(nn)
                case(5);IDP2(kdc,mmd2(kdc),2)=jn(nn)+1
                case(4);IDP2(kdc,mmd2(kdc),2)=in(nn)+1
                case(2);IDP2(kdc,mmd2(kdc),2)=kn(nn)+1
              end select
!!!!!              CINFN2(nn,kdc,IDP2(kdc,mmd2(kdc),1))=mmd2(kdc)  
              nd2(kdc,mmd2(kdc))=icc-1       
    if(nd2(kdc,mmd2(kdc)).gt.mpx2) then
    allocate(dummyr4(1:ndri,1:mmd2mx,1:mpx2,1:3))
    dummyr4=dp2
    deallocate(dp2)
    allocate(dp2(1:ndri,1:mmd2mx,1:nd2(kdc,mmd2(kdc))+5,1:3))
    dp2(1:ndri,1:mmd2mx,1:mpx2,1:3)=dummyr4(1:ndri,1:mmd2mx,1:mpx2,1:3)
    mpx2=nd2(kdc,mmd2(kdc))+5
    deallocate(dummyr4)              
    endif                            
              do kk=1,nd2(kdc,mmd2(kdc))
                dp2(kdc,mmd2(kdc),nd2(kdc,mmd2(kdc))-kk+1,1:3)=pp(icv(kk),0:2)
              enddo
#endif TEST
    ncpl5(0,0)=ncpl5(0,0)+1
                ncpl5(ncpl5(0,0),-3)=kdc
              select case(i)
                case(1);ncpl5(ncpl5(0,0),-2)=3
                case(3);ncpl5(ncpl5(0,0),-2)=2 
                case(6);ncpl5(ncpl5(0,0),-2)=1               
              end select
              select case(i)
                case(3);ncpl5(ncpl5(0,0),-1)=jn(nn)
                case(6);ncpl5(ncpl5(0,0),-1)=in(nn)
                case(1);ncpl5(ncpl5(0,0),-1)=kn(nn)
              end select
                ncpl5(ncpl5(0,0),0)=icc-1 
                ncpl5(ncpl5(0,0),1:icc-1)=icv(1:icc-1)
    continue
enddo

end subroutine for_2d_vis

subroutine realloc_pp_info
integer::iip1
if(ip.gt.ubound(pp,1)) then
    iip1=ubound(pp,1)
    allocate(dummyr2(0:iip1,0:2))
    allocate(dummyi1(0:iip1))    
    dummyr2=pp;dummyi1=info
    deallocate(pp,info)
    allocate(pp(0:ip+20,0:2),info(0:ip+20))
    pp(0:iip1,0:2)=dummyr2(0:iip1,0:2)
    info(0:iip1)=dummyi1(0:iip1);info(iip1+1:ip+20)=0
    deallocate(dummyr2,dummyi1)
endif
end subroutine

subroutine surface_dir3(nnp,icv,cpls,ncp)
doubleprecision,dimension(0:2)::tmpV
doubleprecision,dimension(0:2),intent(in)::cpls
integer,dimension(1:nnp),intent(out)::ncp
integer,dimension(1:nnp),intent(in)::icv
doubleprecision::dd2 !,dis1,dista
integer,intent(in)::nnp
integer::ib,ic

tmpV(0:2)=0.0d0
do ib=2,nnp-1
tmpV=crossV(pp(icv(ib),0:2)-pp(icv(1),0:2),pp(icv(ib+1),0:2)-pp(icv(1),0:2))+tmpV
enddo

dd2=dot_product(tmpV,cpls(0:2))

 if(dd2.gt.0.0d0) then
        ncp(1:nnp)=icv(1:nnp)
 else
        do ic=1,nnp
        ncp(ic)=icv(nnp-ic+1)
        enddo
 endif 

 
end subroutine surface_dir3

end subroutine ss_main_loop

function cal_near_point(A,B,face1,length)
use  interface_list,only:crossv,dista
use variables, only: HALF
implicit none
doubleprecision::cal_near_point(0:2)
doubleprecision,intent(in)::A(0:2),B(0:2),face1(0:3)
doubleprecision,intent(in)::length
doubleprecision,dimension(0:2)::tmpv1,tmpv2,tmpv5

        tmpv1=HALF*(A+B)  ! P
        tmpv2=(B-A)/dista(B-A)  ! AB(s,t,u)
        tmpv5=crossv(face1(0:2),tmpv2)
        cal_near_point=tmpv1+crossv(face1(0:2),tmpv2)*length

end function
