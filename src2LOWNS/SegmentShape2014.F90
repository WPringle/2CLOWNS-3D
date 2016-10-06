subroutine Segment_shape(flag)
!info(ip)=ia>0:ipÇÕï”iaè„ÇÃåì_,ia=-1ÅFipÇÕï®ëÃÇ…ä‹Ç‹ÇÍÇ»Ç¢êÿífñ ÇÃí∏ì_,ia=100ÅFëŒè€ïYó¨ï®ñ è„ÇÃì_
!output (flag = 1 ) nfd,nfbd,axd,ayd,azd,fbd,pp_1,ncpl1,drino,drinoi
!output (flag = 0 ) nf,nfb,ax0,ay0,az0,fb0,pp_0,ncpl0,objno,objnoi
use variables,only:ncpl00max,ncpl01max,ncpl02max,is,ie,js,je,mm0,ndri,inns,inne,n,mm1,mm1f,mm1fmx
use arrays
!#ifdef DRI
!use axayazd
!#endif
!use norm
!use f_p
!use OJT
#ifdef DRI
!use drift2
#ifdef DRI2
!use uvw
#endif
#endif
use omp_lib
implicit none
#ifdef LARGEDATA
integer::nn_1,idiv(0:20),idive,xl1,xl2,nn_nne,iis,iie,jjs,jje,kks,kke,cell_ia_n_p,k,ia_num_max,nnn,nnne,ncount3,cell_ia_2
integer,allocatable,dimension(:,:)::cell_ia,cell_ia_t,cell_ia_n_t
integer,allocatable,dimension(:,:)::nn_nn,idiv_ia
integer,allocatable,dimension(:)::idiv_ia_n,ncount,ncount2,cell_ia_p,cell_ia_n,nn_ia
#endif

integer::j,ia,ib,i,ii,nthr,mythr,nn,ncount0,ipp_p,flag,kk,icc,kdc,k !,iplmax,ic
doubleprecision,allocatable,dimension(:,:,:,:)::dummyr4,dp_p,dp2_p
doubleprecision,allocatable,dimension(:,:,:)::dummyr3,pp_0p
integer,allocatable,dimension(:,:,:)::dummyi3,ncpl0p,info_0p,idp2_p
!doubleprecision,allocatable,dimension(:,:)::dummyr2
integer,allocatable,dimension(:,:)::dummyi2,info_p,mmdp,mmd2p,nd_p,nd2_p
integer,allocatable,dimension(:)::dummyi1,mm0p,objnoip
integer::mm0_p,mmd_p(1:ndri),mmd2_p(1:ndri),mm0mx_p,mmdmx_p,mmd2mx_p,mm_s,mm_e
integer::ncpl2_0,ncpl2_1,ncpl2_2,ncpl2_ii0,ncpl5_ii0
integer,dimension(0:ncpl00max,-2:ncpl02max)::ncpl2
integer,dimension(0:100,-3:ncpl02max)::ncpl5
doubleprecision,allocatable,dimension(:,:)::pp_p
doubleprecision,dimension(inns:inne,0:6)::axyzn

call ss_initial_setting(flag)
      
ncount0=0

!$omp parallel private(mythr,OBJNOIp,mm0_p,mm0mx_p,mmd_p,mmdmx_p,mmd2_p,mmd2mx_p, &
!$omp dummyi1,dummyi2,dummyi3,dummyr3,dummyr4,mm_s,mm_e, &
!$omp ncpl2,ncpl5,ncpl0p,pp_0p,pp_p,ipp_p,info_p,info_0p,ii,kdc,icc, &
#ifdef LARGEDATA
!$omp nn, &
#endif
!$omp dp_p,nd_p,dp2_p,nd2_p,idp2_p)
ipp_p=ncpl01max*2
nthr=omp_get_num_threads();mythr=omp_get_thread_num()
mm0mx_p=int((ie-is)*(je-js)/nthr)
mmdmx_p=int((ie-is)*(je-js)/nthr)
mmd2mx_p=int((ie-is)*(je-js)/nthr)
if(mythr.eq.0) then
    allocate( mm0p(0:nthr-1) );mm0p=0
    allocate( mmdp(1:ndri,0:nthr-1) );mmdp=0
    allocate( mmd2p(1:ndri,0:nthr-1) );mmd2p=0
endif
!$omp barrier

allocate( OBJNOIp(1:mm0mx_p) )  
allocate( ncpl0p(1:mm0mx_p,0:ncpl00max,-2:ncpl02max) )  
allocate( dp_p(1:ndri,1:mmdmx_p,1:ncpl02max,1:3))
allocate( nd_p(1:ndri,1:mmdmx_p))
allocate( idp2_p(1:ndri,1:mmd2mx_p,1:2))
allocate( dp2_p(1:ndri,1:mmd2mx_p,1:ncpl02max,1:3))
allocate( nd2_p(1:ndri,1:mmd2mx_p))
allocate( pp_p(0:ipp_p,0:2))
allocate( pp_0p(1:mm0mx_p,0:ipp_p,0:2))
allocate(info_p(1:2,0:ipp_p))
allocate(info_0p(1:mm0mx_p,1:2,0:ipp_p))

write(6,*) 'mythr=',mythr,nthr,mm0mx_p,ipp_p
!$omp barrier

!!if(mythr.eq.0) then
!deallocate( ncpl0 )
!endif
ncpl2_ii0=0;ncpl5_ii0=0;ncpl2_0=0;ncpl2_1=0;ncpl2_2=0
!$omp do reduction(max:ncpl2_0,ncpl2_1,ncpl2_2,ncpl2_ii0,ncpl5_ii0) !schedule(dynamic) 
nn_loop:do nn=1,inne
#ifdef BUGWIN2
 if(nn.eq.121975) then
     write(6,*) 'nn=',nn,'mythr=',mythr,mm0p(0),mmdp(1:2,0)
     continue
 endif
#endif 
!!!!!!!!!!!!!!!!!!!!!!!!!
if(mod(nn,50000).eq.1.and.flag.eq.0) then
    write(6,*) 'nn=',nn,'mythr=',mythr,mm0p(mythr),mmdp(1,mythr)
endif

mm0_p=0;mmd_p=0;mmd2_p=0 

axyzn(nn,1)=a0(2,nn)
axyzn(nn,2)=a0(2,mn(in(nn),kn(nn)+1,jn(nn)))
axyzn(nn,3)=a0(1,nn)
axyzn(nn,4)=a0(0,mn(in(nn)+1,kn(nn),jn(nn)))
axyzn(nn,5)=a0(1,mn(in(nn),kn(nn),jn(nn)+1))
axyzn(nn,6)=a0(0,nn)
axyzn(nn,0)=fb0(nn)  

call ss_main_loop(flag,nn,mm0_p,mmd_p(1:ndri),mmd2_p(1:ndri),ncpl2,ncpl5,ncount0,pp_p,info_p,ipp_p,axyzn(nn,:) &
#ifdef LARGEDATA
,cell_ia_n_t(nn_1,2),cell_ia_t(nn_1,1:cell_ia_n_t(nn_1,2)) &
#endif
)

if(sum(mmd_p(1:ndri)).ne.0) then
if(mm0_p.eq.0) then
    continue
endif
    do ii=1,ncpl2(0,0)
    if(ncpl2(ii,-1).ne.1) cycle
    if(ncpl2(ii,-2).eq.0) then
        cycle
    endif
    kdc=idri_kdc(ncpl2(ii,-2))
    mmdp(kdc,mythr)=mmdp(kdc,mythr)+1

    if(mmdp(kdc,mythr).gt.mmdmx_p) then
        allocate( dummyr4(1:ndri,1:mmdmx_p,1:ncpl02max,1:3))
        dummyr4=dp_p;deallocate(dp_p)
        allocate( dp_p(1:ndri,1:mmdmx_p*2,1:ncpl02max,1:3))
        dp_p(1:ndri,1:mmdmx_p,1:ncpl02max,1:3)=dummyr4(1:ndri,1:mmdmx_p,1:ncpl02max,1:3)
        deallocate(dummyr4)
        allocate( dummyi2(1:ndri,1:mmdmx_p))
        dummyi2=nd_p;deallocate(nd_p)
        allocate( nd_p(1:ndri,1:mmdmx_p*2))  
        nd_p(1:ndri,1:mmdmx_p)=dummyi2(1:ndri,1:mmdmx_p)
        deallocate(dummyi2)    
        mmdmx_p=mmdmx_p*2
  !      write(6,*) 'mmdmx_p=',mmdmx_p,'mythr=',mythr
    endif
    
    do kk=1,ncpl2(ii,0)
    dp_p(kdc,mmdp(kdc,mythr),kk,1:3)=pp_p(ncpl2(ii,ncpl2(ii,0)-kk+1),0:2)
    if(sum(pp_p(ncpl2(ii,ncpl2(ii,0)-kk+1),0:2)).eq.0.0d0) then
        continue
    endif
    enddo
#ifdef BUGWIN
    if(mmdp(kdc,mythr).eq.36680) then
        continue
    endif
#endif
    nd_p(kdc,mmdp(kdc,mythr))=ncpl2(ii,0) 
    ncpl2_ii0=max(ncpl2_ii0,ncpl2(ii,0))
    enddo
    continue
endif

if(sum(mmd2_p(1:ndri)).ne.0) then
    
         do ii=1,ncpl5(0,0)   
    kdc=ncpl5(ii,-3) 
    mmd2p(kdc,mythr)=mmd2p(kdc,mythr)+1   
    if(mmd2p(kdc,mythr).gt.mmd2mx_p) then

        allocate( dummyi3(1:ndri,1:mmd2mx_p,1:2))
        dummyi3=idp2_p;deallocate(idp2_p)
        allocate( idp2_p(1:ndri,1:mmd2mx_p*2,1:2))
        idp2_p(1:ndri,1:mmd2mx_p,1:2)=dummyi3(1:ndri,1:mmd2mx_p,1:2)
        deallocate(dummyi3)
        
        allocate( dummyr4(1:ndri,1:mmd2mx_p,1:ncpl02max,1:3))
        dummyr4=dp2_p;deallocate(dp2_p)
        allocate( dp2_p(1:ndri,1:mmd2mx_p*2,1:ncpl02max,1:3))
        dp2_p(1:ndri,1:mmd2mx_p,1:ncpl02max,1:3)=dummyr4(1:ndri,1:mmd2mx_p,1:ncpl02max,1:3)
        deallocate(dummyr4)
        
        allocate( dummyi2(1:ndri,1:mmd2mx_p)) 
        dummyi2=nd2_p;deallocate(nd2_p)
        allocate( nd2_p(1:ndri,1:mmd2mx_p*2))
        nd2_p(1:ndri,1:mmd2mx_p)=dummyi2(1:ndri,1:mmd2mx_p)
        deallocate(dummyi2) 
        
        mmd2mx_p=mmd2mx_p*2
  !      write(6,*) 'mmd2mx_p=',mmd2mx_p,'mythr=',mythr        
    endif
    
    IDP2_p(kdc,mmd2p(kdc,mythr),1:2)=ncpl5(ii,-2:-1)
    nd2_p(kdc,mmd2p(kdc,mythr))=ncpl5(ii,0)
    do kk=1,ncpl5(ii,0)
    dp2_p(kdc,mmd2p(kdc,mythr),ncpl5(ii,0)-kk+1,1:3)=pp_p(ncpl5(ii,kk),0:2)
    enddo
    ncpl5_ii0=max(ncpl5_ii0,ncpl5(ii,0))
    continue
        enddo
endif

    
if(mm0_p.eq.0) cycle
mm0p(mythr)=mm0p(mythr)+1
if(mm0p(mythr).gt.mm0mx_p) then

    allocate(dummyi1(1:mm0p(mythr)))
    dummyi1=OBJNOIp
    deallocate(OBJNOIp)
    allocate(OBJNOIp(1:mm0p(mythr)*2))    
    OBJNOIp(1:mm0p(mythr))=dummyi1(1:mm0p(mythr))
    deallocate(dummyi1)  

    allocate(dummyi3(1:mm0mx_p,0:ncpl00max,-2:ncpl02max))
    dummyi3=ncpl0p
    deallocate(ncpl0p)    
    allocate(ncpl0p(1:mm0p(mythr)*2,0:ncpl00max,-2:ncpl02max))
    ncpl0p(1:mm0mx_p,:,:)=dummyi3(1:mm0mx_p,:,:)
    deallocate(dummyi3)

    allocate(dummyr3(1:mm0mx_p,0:ipp_p,0:2))
    dummyr3=pp_0p
    deallocate(pp_0p)    
    allocate(pp_0p(1:mm0p(mythr)*2,0:ipp_p,0:2))
    pp_0p(1:mm0mx_p,:,:)=dummyr3(1:mm0mx_p,:,:)
    deallocate(dummyr3)  
    
    allocate(dummyi3(1:mm0mx_p,1:2,0:ipp_p))
    dummyi3=info_0p
    deallocate(info_0p)    
    allocate(info_0p(1:mm0p(mythr)*2,1:2,0:ipp_p))
    info_0p(1:mm0mx_p,:,:)=dummyi3(1:mm0mx_p,:,:)
    deallocate(dummyi3)  
    
    mm0mx_p=mm0p(mythr)*2
    write(6,*) 'mm0mx_p=',mm0mx_p,'mythr=',mythr
    continue
endif
OBJNOIp(mm0p(mythr))=nn
#ifdef BUGWIN            
            if(nn.eq.2488) then
                write(6,*) 'mythr=',mythr,nn,mm0p(mythr)
             continue
            endif 
#endif 
do i=1,ncpl2(0,1)   
do ia=1,ncpl2(0,0)
do ib=1,ncpl2(ia,0)
if(ncpl2(ia,ib).eq.i) goto 5551
enddo
enddo
pp_p(i,:)=-10000.d0
5551 continue
enddo

icc=0
do i=1,ncpl2(0,1)
    if(pp_p(i,1).eq.-10000.0d0) cycle 
    icc=icc+1
    pp_p(icc,:)=pp_p(i,:)
    info_p(1:2,icc)=info_p(1:2,i)
    do ia=1,ncpl2(0,0)
        do ib=1,ncpl2(ia,0)
        if(ncpl2(ia,ib).eq.i) ncpl2(ia,ib)=icc
        enddo
    enddo
enddo
ncpl2(0,1)=icc
ncpl0p(mm0p(mythr),0:ncpl2(0,0),-2:ncpl2(0,2))=ncpl2(0:ncpl2(0,0),-2:ncpl2(0,2))
pp_0p(mm0p(mythr),1:ncpl2(0,1),0:2)=pp_p(1:ncpl2(0,1),0:2)
info_0p(mm0p(mythr),1:2,1:ncpl2(0,1))=info_p(1:2,1:ncpl2(0,1))
ncpl2_0=max(ncpl2_0,ncpl2(0,0))
ncpl2_1=max(ncpl2_1,ncpl2(0,1))
ncpl2_2=max(ncpl2_2,ncpl2(0,2))
enddo  nn_loop
!$omp end do

if(mythr.eq.0) then
    write(6,*) ncpl2_0,ncpl2_1,ncpl2_2,ncpl2_ii0,ncpl5_ii0
    ncpl00max=max(ncpl2_0,ncpl00max)
    ncpl01max=max(ncpl2_1,ncpl01max)
    ncpl02max=max(ncpl2_2,ncpl02max)
if(.not.drift) then
    mm0=sum(mm0p(0:nthr-1))
    write(6,*) 'end loop mm0=',mm0
    allocate(OBJNOI(1:mm0))
    allocate(ncpl0(1:mm0,0:ncpl2_0,-2:ncpl2_2))
    allocate(pp_0(1:mm0,1:ncpl2_1,0:2))
    allocate(info_0(1:mm0,1:2,1:ncpl2_1))
else
    mm1=sum(mm0p(0:nthr-1))
    write(6,*) 'end loop mm1=',mm1
    allocate(DRINOI(1:mm1))
    allocate(ncpl1(1:mm1,0:ncpl2_0,-2:ncpl2_2))
    allocate(pp_1(1:mm1,1:ncpl2_1,0:2))
    allocate(info_1(1:mm1,1:2,1:ncpl2_1))  
endif
endif

!$omp barrier

mm_s=sum(mm0p(0:mythr-1))+1
mm_e=sum(mm0p(0:mythr))
write(6,*) 'mm0=',mm_s,mm_e,mythr
if(.not.drift) then
    OBJNOI(mm_s:mm_e)=OBJNOIp(1:mm0p(mythr))
    ncpl0(mm_s:mm_e,0:ncpl2_0,-2:ncpl2_2)=ncpl0p(1:mm0p(mythr),0:ncpl2_0,-2:ncpl2_2)
    pp_0(mm_s:mm_e,1:ncpl2_1,0:2)=pp_0p(1:mm0p(mythr),1:ncpl2_1,0:2)
    info_0(mm_s:mm_e,1:2,1:ncpl2_1)=info_0p(1:mm0p(mythr),1:2,1:ncpl2_1)
else
    DRINOI(mm_s:mm_e)=OBJNOIp(1:mm0p(mythr))
    ncpl1(mm_s:mm_e,0:ncpl2_0,-2:ncpl2_2)=ncpl0p(1:mm0p(mythr),0:ncpl2_0,-2:ncpl2_2)
    pp_1(mm_s:mm_e,1:ncpl2_1,0:2)=pp_0p(1:mm0p(mythr),1:ncpl2_1,0:2)
    info_1(mm_s:mm_e,1:2,1:ncpl2_1)=info_0p(1:mm0p(mythr),1:2,1:ncpl2_1)
endif
    deallocate(OBJNOIp,ncpl0p,pp_0p,info_0p)
    
if(mythr.eq.0) then
    do kdc=1,ndri
        mmd(kdc)=sum(mmdp(kdc,0:nthr-1))
    enddo
    write(6,*) 'mmd(1)=',mmd(1:ndri),ncpl2_ii0
    allocate( dp(1:ndri,1:maxval(mmd(1:ndri)),1:ncpl2_ii0,1:3) )
    allocate( nd(1:ndri,1:maxval(mmd(1:ndri)))) 
endif
!$omp barrier    
    

    do kdc=1,ndri
    mm_s=sum(mmdp(kdc,0:mythr-1))+1
    mm_e=sum(mmdp(kdc,0:mythr))        
    dp(kdc,mm_s:mm_e,1:ncpl2_ii0,1:3)=dp_p(kdc,1:mmdp(kdc,mythr),1:ncpl2_ii0,1:3)
    nd(kdc,mm_s:mm_e)=nd_p(kdc,1:mmdp(kdc,mythr))
    enddo
    deallocate(dp_p,nd_p)
    

if(mythr.eq.0) then
    do kdc=1,ndri
    mmd2(kdc)=sum(mmd2p(kdc,0:nthr-1))
    enddo
    allocate( idp2(1:ndri,1:mmd2(1),1:2) )
    allocate( dp2(1:ndri,1:mmd2(1),1:ncpl5_ii0,1:3) )
    allocate( nd2(1:ndri,1:mmd2(1)) )
    continue
endif
!$omp barrier 
 do kdc=1,ndri 
mm_s=sum(mmd2p(kdc,0:mythr-1))+1
mm_e=sum(mmd2p(kdc,0:mythr))
!
    dp2(kdc,mm_s:mm_e,1:ncpl5_ii0,1:3)=dp2_p(kdc,1:mmd2p(kdc,mythr),1:ncpl5_ii0,1:3)
    nd2(kdc,mm_s:mm_e)=nd2_p(kdc,1:mmd2p(kdc,mythr))
    idp2(kdc,mm_s:mm_e,1:2)=idp2_p(kdc,1:mmd2p(kdc,mythr),1:2)
enddo
    deallocate(dp2_p,nd2_p,idp2_p)
    
!if(mythr.eq.0) then
!    continue
!endif
continue
!$omp end parallel
if(.not.drift) then
    allocate(objno(1:inne));objno=0       
    do i=1,mm0
    objno(objnoi(i))=i
    enddo
    do nn=inns,inne
        i=in(nn);j=jn(nn);k=kn(nn)
        if(axyzn(nn,0).eq.0.0d0) then
            continue
        endif
        if(mn(i,k-1,j).ge.inns.and.mn(i,k-1,j).le.inne) then
            a0(2,nn)=min(axyzn(nn,1),axyzn(mn(i,k-1,j),2))
        endif
        if(mn(i,k,j-1).ge.inns.and.mn(i,k,j-1).le.inne) then
            a0(1,nn)=min(axyzn(nn,3),axyzn(mn(i,k,j-1),5))
        endif
        if(mn(i-1,k,j).ge.inns.and.mn(i-1,k,j).le.inne) then
            a0(0,nn)=min(axyzn(nn,6),axyzn(mn(i-1,k,j),4))
        endif
            fb0(nn)=axyzn(nn,0)
    enddo
else
    allocate(drino(1:inne));drino=0         
    do i=1,mm1
    drino(drinoi(i))=i
    enddo
    allocate(drino2i(1:mm1fmx));drino2i=0;mm1f=0     
    do nn=inns,inne
        if(axyzn(nn,0).eq.0.0d0.and.nff(nn)%f.ne.-1) then
            mm1f=mm1f+1
            drino2i(mm1f)=nn
            drino2(nn,1)=mm1f
            nff(nn)%f=-1;nff(nn)%b=10
            continue
        endif
        i=in(nn);j=jn(nn);k=kn(nn)
        if(mn(i,k-1,j).ge.inns.and.mn(i,k-1,j).le.inne) then
            ad(2,nn)=min(axyzn(nn,1),axyzn(mn(i,k-1,j),2))
        endif
        if(mn(i,k,j-1).ge.inns.and.mn(i,k,j-1).le.inne) then
            ad(1,nn)=min(axyzn(nn,3),axyzn(mn(i,k,j-1),5))
        endif
        if(mn(i-1,k,j).ge.inns.and.mn(i-1,k,j).le.inne) then
            ad(0,nn)=min(axyzn(nn,6),axyzn(mn(i-1,k,j),4))
        endif
            fbd(nn)=axyzn(nn,0)
    enddo
endif


#ifdef DRI
if(flag.eq.1) then 
    mm1=mm0
endif    
continue
#endif


contains



subroutine naiten_gaiten(cc,vt1,nvt1,ica)
    implicit none
    integer,intent(in)::nvt1
    doubleprecision,intent(in)::vt1(1:nvt1,0:2),cc(0:2)
    doubleprecision::cc_pc(0:2),aa(1:2,1:2),b_c(0:2),det,a_bb(0:1) !,zh,zv(0:2),xmin,xmax,ymin,ymax
    doubleprecision::bb(1:2,1:2),a_b(0:2),dd(1:2),pc(0:2) !,pmin,pmax,pymin,pymax  !,ff1(1:2),ff2(1:2),ee(1:2,1:2)
!    doubleprecision,dimension(0:1,0:1)::abi,ab0
    integer::ib,icc,ib1,ica
    
!    if(nvt1.eq.0) 
    ica=-1
    icc=0

pc(0)=sum(vt1(1:nvt1,0) )/dfloat(nvt1)
pc(1)=sum(vt1(1:nvt1,1) )/dfloat(nvt1)
pc(2)=sum(vt1(1:nvt1,2) )/dfloat(nvt1)
 

    icc=0
    cc_pc=pc-cc
continue  
    do ib=1,nvt1
        ib1=ib+1
    if(ib.eq.nvt1) ib1=1
    b_c=vt1(ib1,:)-vt1(ib,:)
    a_b=vt1(ib,:)-cc
    aa(1,1)=cc_pc(0);aa(2,1)=cc_pc(1)-cc_pc(2)
    aa(1,2)=(-1.0d0)*b_c(0);aa(2,2)=(-1.0d0)*b_c(1)-(-1.0d0)*b_c(2)
    det=aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)
    if(det.eq.0.0d0) then
            aa(1,1)=cc_pc(1);aa(2,1)=cc_pc(0)-cc_pc(2)
            aa(1,2)=(-1.0d0)*b_c(1);aa(2,2)=(-1.0d0)*b_c(0)-(-1.0d0)*b_c(2)
            det=aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)
            if(det.eq.0.0d0) then
             continue
            else
                bb(1,1)=aa(2,2);bb(2,2)=aa(1,1);bb(2,1)=aa(2,1)*(-1.0d0);bb(1,2)=aa(1,2)*(-1.0d0)
                a_bb(0)=a_b(1);a_bb(1)=a_b(0)-a_b(2)
                dd=MATMUL(bb,a_bb(0:1))/det
            endif
    else
        bb(1,1)=aa(2,2);bb(2,2)=aa(1,1);bb(2,1)=aa(2,1)*(-1.0d0);bb(1,2)=aa(1,2)*(-1.0d0)
        a_bb(0)=a_b(0);a_bb(1)=a_b(1)-a_b(2)
        dd=MATMUL(bb,a_bb(0:1))/det
    endif
        if(dd(2).lt.1.0d0.and.dd(2).ge.0.0d0.and.abs(dd(1)).lt.1.0d-9) then !ÉIÉìÉâÉCÉì
            ica=0
            return
        else if(dd(2).lt.1.0d0.and.dd(2).ge.0.0d0.and.dd(1).ge.0.0d0) then
            icc=icc+1
        endif
    enddo

    if(icc.eq.2.or.icc.eq.4) then
        ica=-1
        return
    else if(icc.eq.1.or.icc.eq.3.or.icc.eq.5) then
        ica=1
        return
    else if(icc.eq.0) then
        write(6,*) 'icc=',icc
        continue
    else
        continue
    endif
end subroutine
end subroutine Segment_shape

subroutine ss_initial_setting(flag)
use arrays,only:rhod,mmd,mmd2,vtmx,vtmn,vt,ivt,xxmax,xxmin,ipl,ipo,npo,drift,drino2
use variables,only:ndri,n,iplmax,mm1f,inne
implicit none
integer::kdc,flag,ia,ib,ic

if(drift) then 
    allocate(drino2(1:inne,1:2));drino2=0 
    mm1f=0                 !ïYó¨ï®ì‡ïî    
endif

#ifdef DRI

#endif
#ifdef DRI2  
mm2p=0                 !ïYó¨ï®Ç…ÇÊÇËï™ífÇ≥ÇÍÇƒÇ¢ÇÈÉZÉãÇÃêî
#endif
#ifndef HENDO
do kdc=1,ndri
if(rhod(kdc).gt.0.0d0.and.flag.eq.0) cycle
if(rhod(kdc).lt.0.0d0.and.flag.eq.1.and.n.gt.0) cycle
mmd(kdc)=0
mmd2(kdc)=0
vtmx(kdc,0)=maxval(vt(kdc,1:ivt(kdc),1))
vtmn(kdc,0)=minval(vt(kdc,1:ivt(kdc),1))
vtmx(kdc,1)=maxval(vt(kdc,1:ivt(kdc),2))
vtmn(kdc,1)=minval(vt(kdc,1:ivt(kdc),2))
vtmx(kdc,2)=maxval(vt(kdc,1:ivt(kdc),3))
vtmn(kdc,2)=minval(vt(kdc,1:ivt(kdc),3))
enddo
#endif

#ifdef DRI
drino2=0
drino=0
#endif

iplmax=maxval(ipl(1:ndri))
allocate(xxmax(1:ndri,1:iplmax,0:2),xxmin(1:ndri,1:iplmax,0:2));xxmax=-10000.d0;xxmin=10000.d0
do kdc=1,ndri
    do ia=1,ipl(kdc)
    !ëŒè€ï®ëÃï\ñ iaÇä‹ÇﬁóÃàÊÇåvéZ
        do ib=1,ipo(kdc,ia)
            do ic=0,2
                xxmax(kdc,ia,ic)=max(xxmax(kdc,ia,ic),vt(kdc,npo(kdc,ia,ib),ic+1))
                xxmin(kdc,ia,ic)=min(xxmin(kdc,ia,ic),vt(kdc,npo(kdc,ia,ib),ic+1))
            enddo
        enddo
    enddo
enddo
!deallocate(OBJNOI,pp_0,info_0,dp,nd,dp2,nd2,idp2)

end subroutine ss_initial_setting