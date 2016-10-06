subroutine set_drift
use variables,only:ndri,mm0mx,mmdmx,mmd2mx,ncpl02max,inne2,inne,title,testa,mm0
use arrays,only:alloc_drift,alloc_for_segment,alloc_adfbd,a0,fb0,a,fb,ad,fbd,nff
implicit none
integer::nn !,n1,i
!!a,a0,fb0,objno,objnoi,pp_0,ncpl0,info_0,,nd,dp,mmd,in,jn,kn,mn,fb,x,alloc_afb,nff
    call alloc_adfbd(inne2,inne)
    call alloc_drift(ndri)
!    call set_cpl_upl

    call alloc_for_segment(ndri,mm0mx,mmdmx,mmd2mx,ncpl02max,inne2)
    ad(:,1:inne)=a0(:,1:inne)
    fbd(1:inne)=fb0(1:inne)
    do nn=1,inne
        if(nff(nn)%b.eq.10) then
                nff(nn)%f=0
                nff(nn)%b=0
        endif
    enddo
    call Segment_shape(1)
    call find_close_area_d
    a(:,1:inne)=ad(:,1:inne)
    fb(1:inne)=fbd(1:inne)
    call set_nflag
end subroutine
    
subroutine find_close_area_d
use variables,only:inns,inne,is,ie,ks,ke,js,je,ZERO
use arrays,only:nff,mn,in,jn,kn,ad,fbd
implicit none
integer,dimension(0:inne)::nflag
integer,dimension(0:inne)::nsearch
integer::nnow,i,j,k,icc,nt,icon,icon1,icon2,nn
! 流入した水が届かない領域（閉鎖領域）を探す
nflag=0;nsearch=0
! 最初の点を決める
do i=is,ie-1
do j=js,je-1
do k=ks,ke-1
if(nff(mn(i,k,j))%f.ge.0) then
    nnow=mn(i,k,j);goto 124
endif    
enddo
enddo
enddo
 
nflag(0)=1
124 nflag(nnow)=1
icc=0
i=in(nnow);j=jn(nnow);k=kn(nnow)
nt=mn(i-1,k,j);if(nff(nt)%f.ne.-1) then
if(ad(0,nnow).gt.ZERO.and.nff(nt)%f.ge.0.and.nflag(nt).ne.1) then
nflag(nt)=1;icc=icc+1;nsearch(icc)=nt;endif;endif
nt=mn(i,k,j-1);if(nff(nt)%f.ne.-1) then
if(ad(1,nnow).gt.ZERO.and.nff(nt)%f.ge.0.and.nflag(nt).ne.1) then
nflag(nt)=1;icc=icc+1;nsearch(icc)=nt;endif;endif
nt=mn(i,k-1,j);if(nff(nt)%f.ne.-1) then
if(ad(2,nnow).gt.ZERO.and.nff(nt)%f.ge.0.and.nflag(nt).ne.1) then
nflag(nt)=1;icc=icc+1;nsearch(icc)=nt;endif;endif
nt=mn(i+1,k,j);if(nff(nt)%f.ne.-1) then
if(ad(0,nt).gt.ZERO.and.nff(nt)%f.ge.0.and.nflag(nt).ne.1) then
nflag(nt)=1;icc=icc+1;nsearch(icc)=nt;endif;endif
nt=mn(i,k,j+1);if(nff(nt)%f.ne.-1) then
if(ad(1,nt).gt.ZERO.and.nff(nt)%f.ge.0.and.nflag(nt).ne.1) then
nflag(nt)=1;icc=icc+1;nsearch(icc)=nt;endif;endif
nt=mn(i,k+1,j);if(nff(nt)%f.ne.-1) then
if(ad(2,nt).gt.ZERO.and.nff(nt)%f.ge.0.and.nflag(nt).ne.1) then
nflag(nt)=1;icc=icc+1;nsearch(icc)=nt;endif;endif
continue

icon2=0
123 icon=0
    icon1=0    
do nnow=1,inne
if(nff(nnow)%f.eq.-1) cycle   !物体は関係ない
icon1=icon1+1
if(nflag(nnow).ne.1) cycle   !nflag=1 の点の隣接点を調べる
icon=icon+1
call findnext_d(nnow,icc)   !nflag=1 の点の隣接点を調べ，連続した隣接点についてnflag=1とする．
if(icc.gt.0) then
    nflag(nnow)=2;icon2=icon2+1  !nflag=1 の点に隣接点があればnflag=2として，その点の検索は終了
endif
enddo
if(icon.ne.0) then
  !  write(6,*) 'icon2=',icon2,icon1   !nflag=1 の点がなくなるまで繰り返す．
    goto 123
endif


continue

do nn=inns,inne
if( nflag(nn).ne.0) cycle
if(nff(nn)%f.eq.-1) cycle
fbd(nn)=ZERO
ad(0,nn)=ZERO;ad(0,mn(in(nn)+1,kn(nn),jn(nn)))=ZERO
ad(1,nn)=ZERO;ad(1,mn(in(nn),kn(nn),jn(nn)+1))=ZERO
ad(2,nn)=ZERO;ad(2,mn(in(nn),kn(nn)+1,jn(nn)))=ZERO
nff(nn)%f=-1
nff(nn)%b=10
enddo

contains



subroutine findnext_d(nnow,icc)
integer::nnow,icc,i,k,j,nt
icc=0
i=in(nnow);j=jn(nnow);k=kn(nnow)
nt=mn(i-1,k,j);if(ad(0,nnow).gt.ZERO.and.nff(nt)%f.ge.0) then;if(nflag(nt).ne.2) nflag(nt)=1;icc=icc+1;endif
nt=mn(i,k,j-1);if(ad(1,nnow).gt.ZERO.and.nff(nt)%f.ge.0) then;if(nflag(nt).ne.2) nflag(nt)=1;icc=icc+1;endif
nt=mn(i,k-1,j);if(ad(2,nnow).gt.ZERO.and.nff(nt)%f.ge.0) then;if(nflag(nt).ne.2)nflag(nt)=1;icc=icc+1;endif
nt=mn(i+1,k,j);if(ad(0,nt).gt.ZERO.and.nff(nt)%f.ge.0) then;if(nflag(nt).ne.2) nflag(nt)=1;icc=icc+1;endif
nt=mn(i,k,j+1);if(ad(1,nt).gt.ZERO.and.nff(nt)%f.ge.0) then;if(nflag(nt).ne.2) nflag(nt)=1;icc=icc+1;endif
nt=mn(i,k+1,j);if(ad(2,nt).gt.ZERO.and.nff(nt)%f.ge.0) then;if(nflag(nt).ne.2) nflag(nt)=1;icc=icc+1;endif
end subroutine findnext_d

end subroutine  
       
subroutine LT(A,R)
        use variables, only: ZERO
        implicit none
        doubleprecision,dimension(3),intent(inout)::A
        doubleprecision,dimension(3)::B
        doubleprecision,dimension(3,3),intent(in)::R
        integer::i,j
        
        B(:)=ZERO
        do i=1,3
          do j=1,3
            B(i)=B(i)+R(i,j)*A(j)
          enddo
        enddo
        A(:)=B(:)
        
end subroutine
        
subroutine RV(DV,theta,R)        
        implicit none
        doubleprecision,dimension(3),intent(in)::DV
        doubleprecision,dimension(3,3),intent(out)::R
        doubleprecision,intent(in)::theta
        doubleprecision::ct,st
        
        ct=cos(theta)
        st=sin(theta)
        
        R(1,1)=DV(1)**2*(1.0d0-ct)+ct
        R(1,2)=DV(1)*DV(2)*(1.0d0-ct)-DV(3)*st
        R(1,3)=DV(3)*DV(1)*(1.0d0-ct)+DV(2)*st
        
        R(2,1)=DV(1)*DV(2)*(1.0d0-ct)+DV(3)*st
        R(2,2)=DV(2)**2*(1.0d0-ct)+ct
        R(2,3)=DV(2)*DV(3)*(1.0d0-ct)-DV(1)*st
        
        R(3,1)=DV(3)*DV(1)*(1.0d0-ct)-DV(2)*st
        R(3,2)=DV(2)*DV(3)*(1.0d0-ct)+DV(1)*st
        R(3,3)=DV(3)**2*(1.0d0-ct)+ct
end subroutine   
        
          
#ifdef DDDDD        
subroutine rotation_and_move(dGD1,theta1)
use arrays,only:ivt,vt,GD,DVEC
implicit none
doubleprecision,dimension(1:3)::theta1,dGD1
doubleprecision,dimension(20,3)::DL
doubleprecision,dimension(3,3)::R
doubleprecision::DVECL
integer::ii,jj
        
        DL(:,:)=ZERO
        do ii=1,ivt(kdc)  ! ii:頂点番号,ivt:頂点数
           DL(ii,:)=vt(kdc,ii,:)-GD(kdc,1:3)
        enddo
        
        do jj=1,3  ! jj:方向ベクトルの番号
          call RV(DVEC(kdc,jj,:),theta1(jj),R)  ! 回転行列Rを計算する
          
          do ii=1,ivt(kdc)  ! ii:頂点番号,ivt:頂点数
            call LT(DL(ii,:),R)
          enddo
          do ii=1,3
            if(ii.ne.jj) then
              call LT(DVEC(kdc,ii,:),R)  ! 方向ベクトルを回転
              DVECL=sqrt(sum(DVEC(kdc,ii,:)**2))
#ifdef BUGWIN
              if(abs(DVECL-1.0d0).gt.1.0d-10) then
                continue
              endif
#endif
              if(DVECL.ne.ZERO) then
                DVEC(kdc,ii,:)=DVEC(kdc,ii,:)/DVECL
              endif
            endif
          enddo
        enddo
        
        ! 移動後の重心と頂点の位置の計算
        GD(kdc,1:3)=GD(kdc,1:3)+dGD1(1:3)
        do ii=1,ivt(kdc)
          vt(kdc,ii,:)=GD(kdc,1:3)+DL(ii,:)
        enddo        

end subroutine rotation_and_move
#endif

 
        