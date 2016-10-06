subroutine drift_pressure  !(tw,dtw,t0)   !(FFD,rot)
use variables,only:ie,je,ke,mm1,mm1f,rho0,g,inns,inne,ndri,ZERO,HALF,t,n,testa
use arrays,only:ffd,rot,in,jn,kn,drino,drinoi,mn,p,nff,drino2,drino2i,x,dx,cpl&
,SFNO,ncpls,pp_s,idri_kdc,rhod,nor,ncpl1,pp_1,md,u
use interface_list,only:crossV,centroid2,crossL
        implicit none
        integer::imkj,ipkj,ikmj,ikpj,ikjm,ikjp
        doubleprecision,dimension(3)::PL !,FV,DL
        integer::kk !,ii,IC,jj,NFD1B
        doubleprecision::S,pd !,d,height
!        doubleprecision,dimension(3)::UV,Vw,Fvis,Vd,Vr
        doubleprecision::dd1,rot22 !,UVL
        integer::ia,kdc  !,kk2
#ifdef BUGWIN
        doubleprecision::DLL1,DLL2
        doubleprecision,dimension(1:10,3)::FFD2
#endif
        doubleprecision,dimension(0:2)::tmpV,pp_11
        doubleprecision,dimension(1:4)::face1,face0        
        doubleprecision,dimension(30,0:2)::cin        
!        doubleprecision,dimension(1:mdri,1:3)::rot  !FFD     
        INTEGER::i,k,j,nn,ib,ib1,ib2,SF1  !,nnp
        
        !!!!!!!!!!!!!!!!!
        integer,allocatable::counter(:,:,:)
        integer:: cell !,n_loop,m,mm,ipd,jpd,kpd
!        doubleprecision::Uave  !ffdd1
 !       doubleprecision,intent(in)::tw,dtw,t0
!        character(6)::times
!        character(7)::dtime
        doubleprecision,dimension(inns:inne)::fbb,fbbp
        doubleprecision,dimension(inns:inne)::pdp,pdpk
        doubleprecision,dimension(3)::pdr,pdl
        doubleprecision::rollp,rollm,pitchp,pitchm,yawp,yawm
!        doubleprecision::Uave_limit
        
!        integer::Scount,sn
!        integer,allocatable::Snum(:,:)
!        doubleprecision,allocatable::vertex(:,:)
!        doubleprecision::face_free(4)
!        doubleprecision::phy,pdy,seg_mass,PLSZ
!        doubleprecision,dimension(3)::PLS,Gdri,Fmass
        doubleprecision,dimension(3)::DL1,DL2,FFD1
        doubleprecision,dimension(1:ndri,1:3)::FFD0
        integer :: kdc1,ncpl1kkia0
        integer :: kdc11, kdc12,nss,nee,nslash1
        real(8) :: Fvis_all(ndri, 1:3)
        !------------------------------!
!        integer :: nnn, kk2, num_vt
!        real(8) :: diff
        !------------------------------!
        
        if(n.eq.0) then
            call owari(testa,nee); call hajime(testa,nss); call sura(testa,nslash1)
            write(6,*) testa(nss:nee)//"/DAT/"//testa(nslash1+1:nee)//".force.csv",nslash1
            open(unit=666,file=testa(nss:nee)//"/DAT/"//testa(nslash1+1:nee)//".force.csv",status='unknown')
        endif
        allocate(counter(ie+2,ke+2,je+2))
        !counter(:,:,:)=0
        !n_loop=int(t/dt)
        !!!!!!!!!!!!!!!!!
        !FFD(kdc,:)=ZERO
        !rot(kdc,:)=ZERO
        FFD(:,:)=ZERO
        rot(:,:)=ZERO
        Fvis_all(:, :) = ZERO
#ifdef BUGWIN2
        FFD2(kdc,:)=ZERO
#endif
        
        ! face1 ñ ÇÃäOå¸Ç´ñ@ê¸

fbb(:)=1.0d0
pdpk(:)=ZERO
pdr(:)=ZERO;pdl(:)=ZERO

do kk=1,mm1f
nn=drino2i(kk)
    !if(nn == 93982) then
    !    continue
    !endif
i=in(nn);k=kn(nn);j=jn(nn)
imkj=mn(i-1,k,j);FFD0=0.0
if(nff(imkj)%f.ge.1.and.nff(imkj)%b.ne.-1) then
if(drino(imkj).eq.0) then
    !if(counter(i-1,k,j).eq.1) goto 1027
    !kdc=1;
    cell=1
    PL(1)=x(0,i);PL(2)=HALF*(x(1,j)+x(1,j+1));PL(3)=HALF*(x(2,k)+x(2,k+1))
    face1=cpl(6,:);face1(4)=-sum(face1(1:3)*PL(1:3))
    s=dx(2,k)*dx(1,j);pd=p(imkj) !phy=p(mn(i-1,k,j))
    !call find_surface(PL, kdc1)
    kdc1 = drino2(nn, 2)
    call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
    pdp(nn)=pdpk(nn)
    cell=0
    !counter(i-1,k,j)=1
endif
endif
ipkj=mn(i+1,k,j)
if(nff(ipkj)%f.ge.1.and.nff(ipkj)%b.ne.1) then
if(drino(ipkj).eq.0) then
    !if(counter(i+1,k,j).eq.1) goto 1027
    !kdc=1;
    cell=2
    PL(1)=x(0,i+1);PL(2)=HALF*(x(1,j)+x(1,j+1));PL(3)=HALF*(x(2,k)+x(2,k+1))
    face1=cpl(4,:);face1(4)=-sum(face1(1:3)*PL(1:3))
    s=dx(2,k)*dx(1,j);pd=p(ipkj) !phy=p(mn(i+1,k,j))
    !call find_surface(PL, kdc1)
    kdc1 = drino2(nn, 2)
    call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
    pdp(nn)=pdpk(nn)
    cell=0
    !counter(i+1,k,j)=1
    continue
endif
endif
ikmj=mn(i,k-1,j)
if(nff(ikmj)%f.ge.1.and.nff(ikmj)%b.ne.-3) then
if(drino(ikmj).eq.0) then
    !kdc=1;
    cell=8
    PL(1)=0.5*(x(0,i)+x(0,i+1));PL(2)=HALF*(x(1,j)+x(1,j+1));PL(3)=x(2,k)
    face1=cpl(1,:);face1(4)=-sum(face1(1:3)*PL(1:3))
    s=dx(0,i)*dx(1,j);pd=p(ikmj)-g*HALF*dx(2,k-1)  !phy=p(mn(i,k-1,j))-g*HALF*dx(2,k-1) 
    !call find_surface(PL, kdc1)
    kdc1 = drino2(nn, 2)
    call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
    !pdp(mn(i,k,j))=pdpk(mn(i,k,j))
    !cell=0
    continue
endif
endif
ikpj=mn(i,k+1,j)
if(nff(ikpj)%f.ge.1.and.nff(ikpj)%b.ne.3) then
if(drino(ikpj).eq.0) then
    !kdc=1;
    cell=9
    PL(1)=0.5*(x(0,i)+x(0,i+1));PL(2)=HALF*(x(1,j)+x(1,j+1));PL(3)=x(2,k+1)
    face1=cpl(2,:);face1(4)=-sum(face1(1:3)*PL(1:3))
    s=dx(0,i)*dx(1,j);pd=p(mn(i,k+1,j))+g*HALF*dx(2,k+1) !phy=p(mn(i,k+1,j))+g*HALF*dx(2,k+1)
    !call find_surface(PL, kdc1)
    kdc1 = drino2(nn, 2)
    call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
    !pdp(mn(i,k,j))=pdpk(mn(i,k,j))
    !cell=0
    continue
endif
endif
ikjp=mn(i,k,j+1)
if(nff(ikjp)%f.ge.1.and.nff(ikjp)%b.ne.2) then
if(drino(ikjp).eq.0) then
    !kdc=1;
    cell=6
    PL(1)=0.5*(x(0,i)+x(0,i+1));PL(2)=x(1,j+1);PL(3)=0.5*(x(2,k)+x(2,k+1))
    face1=cpl(5,:);face1(4)=-sum(face1(1:3)*PL(1:3))
    s=dx(0,i)*dx(2,k);pd=p(ikjp) !phy=p(mn(i,k,j+1))
    !call find_surface(PL, kdc1)
    kdc1 = drino2(nn, 2)
    call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
    !pdp(mn(i,k,j))=pdpk(mn(i,k,j))
    !cell=0
endif
endif
ikjm=mn(i,k,j-1)
if(nff(ikjm)%f.ge.1.and.nff(ikjm)%b.ne.-2) then
if(drino(ikjm).eq.0) then
    !kdc=1;
    cell=7
    PL(1)=0.5*(x(0,i)+x(0,i+1));PL(2)=x(1,j);PL(3)=0.5*(x(2,k)+x(2,k+1))
    face1=cpl(3,:);face1(4)=-sum(face1(1:3)*PL(1:3))
    s=dx(0,i)*dx(2,k);pd=p(ikjm) !phy=p(mn(i,k,j-1))
    !call find_surface(PL, kdc1)
    kdc1 = drino2(nn, 2)
    call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
    !pdp(mn(i,k,j))=pdpk(mn(i,k,j))
    !cell=0
endif
endif

!1027 continue
         FFD(:,:)=FFD(:,:)+FFD0(:,:)
     
enddo        

kdc11 = 0; kdc12 = 0
        
mm1_loop:  do kk=1,mm1
          nn=DRINOI(kk)         
          i=in(nn)
          j=jn(nn)
          k=kn(nn)
          if(nff(nn)%f.eq.0) cycle
          
          rot22=rot(1,2)   
       
#ifdef BUGWIN          
if(nn.eq.164100) then
continue
endif 
#endif     
FFD0=0.0
!êÖñ Çä‹ÇﬁèÍçá
if(SFNO(nn).gt.0.and.nff(nn)%f.ne.1) then
    SF1=SFNO(nn)

    do ia=1,ncpls(SF1,0,0)
        if(ncpls(SF1,ia,-1).ne.1) cycle    
        if(ncpls(SF1,ia,-1).eq.1.and.ncpls(SF1,ia,-2).eq.0) cycle
        if(ncpls(SF1,ia,0).lt.3) cycle
        kdc1=idri_kdc(ncpls(SF1,ia,-2))
        !=============================================================!
#ifdef WM
        if(kdc1 == 0) cycle
#endif
        !=============================================================!
        if(rhod(kdc1).le.ZERO) cycle

            tmpV=ZERO
            do ib=2,ncpls(SF1,ia,0)-1
            ib1=ncpls(sf1,ia,ib)
            ib2=ncpls(sf1,ia,ib+1)
            tmpV=crossV(pp_s(sf1,ib1,0:2)-pp_s(sf1,ncpls(sf1,ia,1),0:2),pp_s(sf1,ib2,0:2)-pp_s(sf1,ncpls(sf1,ia,1),0:2))+tmpV
            continue
            enddo
            s=sqrt(sum(tmpV**2))*HALF       
            do ib=1,ncpls(SF1,ia,0)
            cin(ib,0:2)=pp_s(sf1,ncpls(sf1,ia,ib),0:2)
            enddo
            PL(:)=centroid2(cin(1:ncpls(SF1,ia,0),0:2),ncpls(SF1,ia,0))
            dd1=sum(nor(nn,0:2)*PL(:))+nor(nn,3)
            !!!!!!!!!!!!!!!!!!!!!!
            !if(dot_product(tmpG(0:2),PL(1:3)).gt.ZERO) cycle
            !pd=p(mn(i,k,j))+g*(-1.0d0)*dd1
            !!!!!!!!!!!!!!!!!!!!!!
            pd=g*(-1.0d0)*dd1
            !phy=g*(-1.0d0)*dd1
#ifdef BUGWIN
      if(pd.lt.ZERO) then
      continue
      endif      
#endif            
            face0=cpl(ncpls(SF1,ia,-2),:)*(-1.0d0)
            DL1(1:3)=pp_s(sf1,ncpls(sf1,ia,2),0:2)-pp_s(sf1,ncpls(sf1,ia,1),0:2)
            DL2(1:3)=pp_s(sf1,ncpls(sf1,ia,ncpls(sf1,ia,0)),0:2)-pp_s(sf1,ncpls(sf1,ia,1),0:2)
            call cal_face1(DL1,DL2,face1)
            face1(4) = -sum(PL(1:3) * face1(1:3))
!#ifdef DRI
            cell=3
            call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)  !(pd,s)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
            pdp(nn)=pdpk(nn)
            cell=0
!#else            
!            FFD(kdc1,:)=FFD(kdc1,:)+pd*rho0*S*(-face1(1:3))
!#endif
            continue
                
    enddo
    
else
!#ifdef DRI
    do ia=1,ncpl1(kk,0,0)
        if(ncpl1(kk,ia,-1).ne.1) cycle   
        if(ncpl1(kk,ia,-1).eq.1.and.ncpl1(kk,ia,-2).eq.0) cycle             
        kdc1=idri_kdc(ncpl1(kk,ia,-2))
        !=============================================================!
#ifdef WM
        if(kdc1 == 0) cycle
#endif
        !=============================================================!
        if(rhod(kdc1).le.ZERO) cycle
            ncpl1kkia0=ncpl1(kk,ia,0)
            if(ncpl1kkia0.lt.3) cycle
            tmpV=ZERO
            pp_11=pp_1(kk,ncpl1(kk,ia,1),0:2)
            do ib=2,ncpl1kkia0-1
            ib1=ncpl1(kk,ia,ib)
            ib2=ncpl1(kk,ia,ib+1)
            tmpV=crossV(pp_1(kk,ib1,0:2)-pp_11,pp_1(kk,ib2,0:2)-pp_11)+tmpV
            continue
            enddo
            s=sqrt(sum(tmpV**2))*HALF       
            do ib=1,ncpl1kkia0
            cin(ib,0:2)=pp_1(kk,ncpl1(kk,ia,ib),0:2)
            enddo
            PL(:)=centroid2(cin(1:ncpl1kkia0,0:2),ncpl1kkia0)           
            pd=p(nn)+g*((x(2,k)+x(2,k+1))/2.0d0-PL(3)) 
            
            face0=cpl(ncpl1(kk,ia,-2),:)*(-1.0d0)
            DL1(1:3)=pp_1(kk,ncpl1(kk,ia,2),0:2)-pp_11
            DL2(1:3)=pp_1(kk,ncpl1(kk,ia,ncpl1(kk,ia,0)),0:2)-pp_11
            call cal_face1(DL1,DL2,face1)
            face1(4) = -sum(PL(1:3) * face1(1:3))
!#ifdef DRI
            cell=4
            call nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)  !(pd,s)
    FFD0(kdc1,:)=FFD0(kdc1,:)+FFD1
            write(621,*) p(nn),pl(3),s
            pdp(nn)=pdpk(nn)
            cell=0
!#else            
!            FFD(kdc1,:)=FFD(kdc1,:)+pd*rho0*S*(-face1(1:3))
!#endif
    
    enddo
endif  
!write(622,*) i,k,ffd0(1),ffd0(3)
continue

    FFD(:,:)=FFD(:,:)+FFD0(:,:)
!#ifdef BUGWIN
!if(t.ge.3.d0)  write(65,431) i,k,j,rot(2,2)-rot22,pd,t,pl(:),face1(:)
431 format(1h ,i3,2(',',i3),10(',',e12.5))
!#endif   
enddo  mm1_loop

!write(6, *) kdc11, kdc12
!write(6, *) 'FFD'
!write(6, *) (FFD(kk, 1:3), kk = 1, ndri)
!write(6, *) 'Fvis'
!write(6, *) (Fvis_all(kk, 1:3), kk = 1, ndri)
!write(6, *) 'Fvis / FFD'
!write(6, *) ((Fvis_all(kk, 1:3) / FFD(kk, 1:3)), kk = 1, ndri)

#ifdef DDD
#endif
do kdc = 1, ndri
    FFD(kdc,3)=FFD(kdc,3)-md(kdc)*g
enddo


#ifdef DRI 
#ifdef FLAPGATE
rot(kdc,1)=ZERO;  rot(kdc,3)=ZERO
write(63,*) t,rot(kdc,2),md(kdc)*(GD(kdc,1)-fixp(kdc,1))*g
rot(kdc,2)=rot(kdc,2)+md(kdc)*(GD(kdc,1)-fixp(kdc,1))*g
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
666 format(3(e12.5))
#ifdef DDDD

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!2014/04/24!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
where(abs(FFD(:,:)) < 1.0d-10) FFD=ZERO
where(abs(rot(:,:)) < 1.0d-10) rot=ZERO
!write(*,*) n
if(ffd(1,1).ne.0.0d0) then
    write(*,*) 'FFD1=',t,FFD(1,:)
    write(666,662) t,FFD(1,:)
endif
662 format(1h ,e12.4,3(',',e12.4))
!write(*,*) 'FFD2=',FFD(2,:)
!write(*,*) 'ROT=',rot
!close(27)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains      

!===================== 2014/05/20 ===========================================
subroutine cal_Vw2(Pin,Vout)
!use UVW
!use F_P
!use UNVNWN
implicit none
doubleprecision,dimension(3),intent(in)::Pin
doubleprecision,dimension(3),intent(out)::Vout
doubleprecision,dimension(64,3)::Vcell,rcell
doubleprecision,dimension(3)::rv
doubleprecision::rr,rsum
integer,dimension(3)::dir
integer::i2,j2,k2,is2,ie2,js2,je2,ks2,ke2
integer::counter,nn1 !,nn2

i2=in(nn);j2=jn(nn);k2=kn(nn)
call OD(dir)

if(dir(1).eq.-1) then
    is2=i2-2;ie2=i2
elseif(dir(1).eq.1) then
    is2=i2;ie2=i2+2
else
    is2=i2-1;ie2=i2+1
endif

if(dir(2).eq.-1) then
    js2=j2-2;je2=j2
elseif(dir(2).eq.1) then
    js2=j2;je2=j2+2
else
    js2=j2-1;je2=j2+1
endif

if(dir(3).eq.-1) then
    ks2=k2-2;ke2=k2
elseif(dir(3).eq.1) then
    ks2=k2;ke2=k2+2
else
    ks2=k2-1;ke2=k2+1
endif

counter=0
do j2=js2,je2
    do k2=ks2,ke2
        do i2=is2,ie2
            nn1=mn(i2,k2,j2)
            !if(drino2(nn1,1).ne.0) cycle
            if(nff(nn1)%f.le.0) cycle
            counter=counter+1
            Vcell(counter,1)=HALF*(u(0,nn1)+u(0,mn(i2+1,k2,j2)))
            Vcell(counter,2)=HALF*(u(1,nn1)+u(1,mn(i2,k2,j2+1)))
            Vcell(counter,3)=HALF*(u(2,nn1)+u(2,mn(i2,k2+1,j2)))
            rcell(counter,1:3)=(/PL(1)-(x(0,i2)+HALF*dx(0,i2)),PL(2)-(x(1,j2)+HALF*dx(1,j2)),PL(3)-(x(2,k2)+HALF*dx(2,k2))/)
        enddo
    enddo
enddo

rsum=ZERO;rv(:)=ZERO
do i2=1,counter
    rr=sqrt(sum(rcell(i2,:)**2))
    rsum=rsum+1.0d0/rr
    rv(:)=rv(:)+Vcell(i2,:)/rr
enddo

if(rsum.lt.1.0d-10) then
    write(*,*) nn
    write(*,*) 'rsum is nearly ZERO!'
endif
Vout(:)=rv(:)/rsum
!write(*,*) 'Vout from cal_Vw2'
!write(*,*) Vout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!open(20,file='uvw2.txt')
!do i2=inns,inne
!    write(20,*) u(i2),v(i2),w(i2)
!enddo
!close(20)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine cal_Vw2

subroutine OD(obj_d)
implicit none
integer::i2,j2,k2,mm
integer,dimension(6)::obj !1:k-1,2:k+1,3:j-1,4:i+1,5:j+1,6:i-1
integer,dimension(3),intent(out)::obj_d

obj(:)=0;obj_d(:)=0
i2=in(nn);j2=jn(nn);k2=kn(nn)


mm=mn(i2,k2-1,j2);if(drino2(mm,1).ne.0) obj(1)=1
mm=mn(i2,k2+1,j2);if(drino2(mm,1).ne.0) obj(2)=1
mm=mn(i2,k2,j2-1);if(drino2(mm,1).ne.0) obj(3)=1
mm=mn(i2+1,k2,j2);if(drino2(mm,1).ne.0) obj(4)=1
mm=mn(i2,k2,j2+1);if(drino2(mm,1).ne.0) obj(5)=1
mm=mn(i2-1,k2,j2);if(drino2(mm,1).ne.0) obj(6)=1

if(obj(1).gt.obj(2)) obj_d(3)=1
if(obj(1).lt.obj(2)) obj_d(3)=-1
if(obj(3).gt.obj(5)) obj_d(2)=1
if(obj(3).lt.obj(5)) obj_d(2)=-1
if(obj(6).gt.obj(4)) obj_d(1)=1
if(obj(6).lt.obj(4)) obj_d(1)=-1

end subroutine OD

!===========================================================================
!===================== 2014/07/02 ==========================================
subroutine cal_face1(DL1,DL2,face1)
implicit none
doubleprecision,dimension(3),intent(in)::DL1,DL2
doubleprecision,dimension(4),intent(out)::face1
doubleprecision::DLL

face1=ZERO
call cross(DL1,DL2,face1(1:3))
DLL=crossL(DL1,DL2)
face1(1:3)=-face1(1:3)/DLL

end subroutine cal_face1
!==========================================================================

!!#ifdef DRI

subroutine ffd1_fix(Uave)
use variables,only:js,ks
use arrays,only:nff,mn,u,drinoi,drino2i
implicit none
integer i,ii,jj,kk,min,c_da !,j,k
doubleprecision,intent(out)::Uave
!doubleprecision,intent(out)::ffdd1

min=ie
do ii=1,mm1
    jj=drinoi(ii)
    kk=drino2i(ii)
    if(jj.eq.kk) cycle
    i=in(jj)
    if(i.le.min)then
        min=i
    endif
enddo

Uave=ZERO;c_da=0
do jj=js,je
    do kk=ks,ke
        if(nff(mn(min-5,kk,jj))%f.le.0) cycle
        Uave=Uave+u(0,mn(min-5,kk,jj))
        c_da=c_da+1
    enddo
enddo

if(c_da.eq.0) goto 88
Uave=Uave/dble(c_da)
!write(555,*) min
88 continue
   
!if(abs(Uave).lt.1.0d-3)then
!    ffdd1=ZERO
!else
!    ffdd1=-1.0d0
!endif

end subroutine ffd1_fix
!!#endif  

subroutine cross_points(mm,sn,Snum,Scount)
implicit none
integer ia,ib,ic
integer isurf !,idri
integer,intent(in)::mm,sn
integer,intent(inout)::Snum(sn,0:1)
integer,intent(inout)::Scount

!allocate(Snum(sn,0:1))
!#ifdef DDD
do ia=1,ncpls(mm,0,0)
    if(ncpls(mm,ia,-1).eq.2)then
        isurf=ia
    endif
enddo
do ia=1,ncpls(mm,0,0)
    if(ncpls(mm,ia,-1).eq.1) then
        do ib=1,ncpls(mm,ia,0)
            do ic=1,ncpls(mm,isurf,0)
                if(ncpls(mm,ia,ib).eq.ncpls(mm,isurf,ic)) then
                    Scount=Scount+1
                    Snum(Scount,0)=mm
                    Snum(Scount,1)=ncpls(mm,ia,ib)
                endif
            enddo
        enddo
    endif
enddo
!#endif
#ifdef DDD
do ia=1,ncpls(mm,0,0)
    if(ncpls(mm,ia,-1).eq.2) then
        do ib=1,ncpls(mm,ia,0)
            Scount=Scount+1
            Snum(Scount,0)=mm
            Snum(Scount,1)=ncpls(mm,ia,ib)
        enddo
    endif
enddo
#endif

    end subroutine cross_points
    
    subroutine segg(kk,Gdri,seg_mass)
    implicit none
    integer,intent(in)::kk
    integer::ggn,ia,ib,ic,ic1,ic2,id
    doubleprecision,allocatable::shimen(:,:)
    doubleprecision::dd1,dd2,VV,Vcell
    doubleprecision,dimension(3)::Ss,Gs,Gcell,GV
    doubleprecision,dimension(3),intent(out)::Gdri
    doubleprecision,intent(out)::seg_mass
    
    !! For allocation
    ggn=0
    do ia=1,ncpl1(kk,0,0)
        ib=ncpl1(kk,ia,0)
        do ic=2,ib-1
            ggn=ggn+1
        enddo
    enddo
    !!
    allocate(shimen(ggn,1:4))
    
    VV=ZERO;GV(:)=ZERO;ggn=0
    do ia=1,ncpl1(kk,0,0)
        ib=ncpl1(kk,ia,0)
    if(ib.eq.0) cycle
    Ss=ZERO;Gs(:)=ZERO
    do ic=2,ib-1
        ic1=ncpl1(kk,ia,ic)
        ic2=ncpl1(kk,ia,ic+1)
        
        Ss=crossV(pp_1(kk,ic1,0:2)-pp_1(kk,ncpl1(kk,ia,1),0:2),pp_1(kk,ic2,0:2)-pp_1(kk,ncpl1(kk,ia,1),0:2))
        Gs(1:3)=(pp_1(kk,ncpl1(kk,ia,1),0:2)+pp_1(kk,ic1,0:2)+pp_1(kk,ic2,0:2))/3.0d0
        dd2=dot_product(cpl(ncpl1(kk,ia,-2),0:2),Ss)
        !if(sqrt(sum(Ss**2)).gt.1.0d-8) dd2=dd2/abs(dd2)
        if(dd2.lt.ZERO) dd2=-1.0d0;if(dd2.gt.1.0d0) dd2=1.0d0
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ggn=ggn+1
        dd1=dot_product(cpl(ncpl1(kk,ia,-2),0:2)/sqrt(sum(cpl(ncpl1(kk,ia,-2),0:2)**2)),Gs(1:3))
        
        shimen(ggn,4)=dd2*sqrt(sum(Ss**2))*dd1/6.0d0
        do id=1,3
            shimen(ggn,id)=(pp_1(kk,ic1,id-1)+pp_1(kk,ic2,id-1)+pp_1(kk,ncpl1(kk,ia,1),id-1))/4.0d0
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        continue
    enddo
    continue  
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    VV=ZERO
    do ia=1,ggn
        VV=VV+shimen(ia,4)
        GV(1:3)=GV(1:3)+shimen(ia,4)*shimen(ia,1:3)
    enddo
    GV(1:3)=GV(1:3)/VV
    Vcell=dx(0,i)*dx(1,j)*dx(2,k)
    Gcell(1:3)=(/x(0,i)+HALF*dx(0,i),x(1,j)+HALF*dx(1,j),x(2,k)+HALF*dx(2,k)/)
    !DVP=(Vcell-VV)/Vcell;VVP=VV/Vcell
    do ia=1,3
        !Gdri(jj)=((DVP+VVP)*Gcell(jj)-(VVP)*GV(jj))/DVP
        Gdri(ia)=(Vcell*Gcell(ia)-VV*GV(ia))/(Vcell-VV)
    enddo
   
    seg_mass=rhod(kdc)*(Vcell-VV)
    
    end subroutine segg
    
    !==================== FOR MULTIPLE FLOATING OBJECTS ==========================!
    subroutine find_surface(PL, kdc1)
    use arrays,only:GD
    implicit none
    real(8), intent(in) :: PL(1:3)
    integer, intent(out) :: kdc1
    
    integer :: i !, j, k
    real(8) :: dist
    real(8) :: dis_min
    
    if(ndri == 1) then
        kdc1 = 1
        return
    endif
    
    kdc1 = 0; dis_min = 1.0d8
    do i = 1, ndri
        dist = sqrt(sum((GD(i, 1:3) - PL(1:3)) ** 2))
        if(dist < dis_min) then
            dis_min = dist
            kdc1 = i
            !write(6, *) kdc1
        endif
    enddo
    if(kdc1 == 0) then
        write(6, *) 'Cannot find the surface of floating object!'
        !write(6, *) ndri, dis_min, kdc1
        !write(6, *) (GD(i, 1:3), i = 1, ndri)
        write(6, *) (PL(i), i = 1, 3)
        write(6, *) cell
        stop
    endif
    
    end subroutine find_surface
    !
    !subroutine calc_S_PL()
    !    implicit none
    !end subroutine calc_S_PL
    
        
        
    
    
    
    end subroutine
    
    
subroutine LS(X,A,B)
use variables, only: ZERO
implicit none
doubleprecision,dimension(3),intent(in)::X,A
doubleprecision,dimension(3),intent(out)::B
doubleprecision::s,BL,XL
integer::ii !,i
        
s=sum(X(:)*A(:))
!   write(6,*) 's=',s
!   write(6,*) 'x=',x        
!   write(6,*) 'a=',a          
B(:)=X(:)-s*A(:)
BL=sqrt(sum(B(:)**2))
!        write(6,*) 'bl=',bl   
if(BL.eq.ZERO) return
XL=sqrt(sum(X(:)**2))
      
do ii=1,3
    if(abs(B(ii)/BL-1.0d0).lt.1.0d-8.and.BL/XL.lt.1.0d-8) then
    B(ii)=ZERO
    endif
enddo
        
end subroutine
subroutine RotationalMoment(PL,Fseg,kdcc) !,drot)
use arrays,only:rot,DVECp
implicit none
integer,intent(in)::kdcc
doubleprecision,dimension(3),intent(in)::PL,Fseg
doubleprecision,dimension(3)::Tq !,PLt,Fsegt
!doubleprecision,dimension(1:3),intent(out)::drot
integer::ii
        
!drot(:)=ZERO 
do ii=1,3  ! ii:iié≤
    call cross(PL(:),Fseg(:),Tq(:))
    rot(kdcc,ii)=rot(kdcc,ii)+sum(DVECp(kdcc,ii,:)*Tq(:))
    !drot(ii)=sum(DVECp(kdcc,ii,:)*Tq(:))
enddo

end subroutine   
    
subroutine nenseihoka(nn,kdc1,cell,pdp,pdpk,fbb,fbbp,pdr,pdl,rollp,rollm,pitchp,pitchm,yawp,yawm,PL,FFD1,pd,s,rho0,face1)   !(pd,S)
use arrays,only:DVECp,omegap,GDp,VVDp,DX,in,jn,kn
use variables,only:t,js,ks,inns,inne,ZERO
implicit none
integer::i,j,k
integer,intent(in)::kdc1,cell,nn
doubleprecision::UVL !,rot0
doubleprecision,dimension(3)::Vd,Fvis,Vr,Vw !,Vw0,Vr0
!doubleprecision,dimension(3)::drot
doubleprecision,dimension(3),intent(inout)::pdr,pdl
doubleprecision,intent(inout)::rollp,rollm,pitchp,pitchm,yawp,yawm
doubleprecision::pdk,pdd,pd_c
doubleprecision,dimension(inns:inne),intent(inout)::pdp,pdpk
doubleprecision,dimension(inns:inne),intent(inout)::fbb,fbbp
integer::ii !,kk2,ss,tt,uu
!!!!!!!!!!!!!!!!!!!!!!!!!
!doubleprecision::Cd=1.12d0
doubleprecision,dimension(3)::Fdrag,DLn !,Vs,Vse,cs,Se,evec1,evec2,evec3
!doubleprecision,dimension(2)::dir
doubleprecision,dimension(4)::dxyz
doubleprecision::Re,Vwa
!doubleprecision,dimension(3)::pdy2
!==========================================================================!
doubleprecision,dimension(1:3),intent(inout)::FFD1
doubleprecision,dimension(3),intent(inout)::PL
doubleprecision,dimension(4),intent(in)::face1
doubleprecision,intent(in)::pd,s,rho0
doubleprecision::anu=1.0d-6
!==========================================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef BUGWIN
doubleprecision::pd1,s1
pd1=pd;s1=s
#endif
!        kdc1 = drino2(nn, 2)
          !pdd=pd
          i=in(nn);j=jn(nn);k=kn(nn)  
          ! îSê´óÕÇÃåvéZ
          Fvis(:)=ZERO
          Vd(:)=ZERO
          UVL=min(dx(0,i),dx(1,j),dx(2,k))/1000.0d0     !1.0d-2
          
          if(kdc1 == 0) then
              write(6, *) cell
              write(6, *) i, k, j
          endif
          
          ! ïYó¨ï®ï\ñ ÇÃë¨ìxÇÃåvéZ
          do ii=1,3
#ifdef FLAPGATE
            call cross(DVECp(kdc1,ii,:)*omegap(kdc1,ii),PL(:)-fixp(kdc1,1:3),DLn(:)) ! cross2 Ç©ÇÁcross Ç…ïœçX
#else
            call cross(DVECp(kdc1,ii,:)*omegap(kdc1,ii),PL(:)-GDp(kdc1,1:3),DLn(:)) ! cross2 Ç©ÇÁcross Ç…ïœçX
#endif
 !           
          Vd(:)=Vd(:)+DLn(:)     
          enddo
          Vd(:)=Vd(:)+VVDp(kdc1,:)     !ï®ëÃÇÃâÒì]Ç∆ï¿êië¨ìxÇÃòa
          
       
          Vwa=sqrt(sum(Vw(:)**2))
          dxyz(1:3)=(/dx(0,i),dx(1,j),dx(2,k)/)
          Dxyz(4)=sqrt(sum(dxyz(1:3)**2))
          Re=abs(Vwa*dxyz(4)/anu)
          
          !wriTE(*,*) Remax
          
          Fdrag(:)=ZERO;Fvis(:)=ZERO
          !if(Re.gt.100.0d0)then
#ifdef DDD          
 
#endif          
          !else
          
#ifdef DDDD          

#endif
!#ifdef DDD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call cal_Vw(nn,PL(:)+UVL*face1(1:3),Vw)
            !call cal_Vw(PL(:),Vw0)
            call LS(Vw(:)-Vd(:),face1(1:3),Vr(:))
            !call LS(Vw0(:)-Vd(:),face1(1:3),Vr0(:))
            !Vr(:)=Vr(:)-Vr0(:)
            Fvis(:)=anu*rho0*S*Vr(:)/UVL
            
            !Vs(:)=(Vw(:)-Vd(:))-Vr(:)
            !Fdrag(:)=HALF*rho0*Cd*S*Vs(:)*abs(Vs(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#endif            
#ifdef DDD        

#endif
          
#ifdef BUGWIN2
          write(6,*) Vr,Fvis,pd*rho0*S*(-face1(:))
          continue
#endif
!          
#ifdef FLAPGATE
            PL(:)=PL(:)-fixp(kdc1,1:3)       
#else
            PL(:)=PL(:)-GDp(kdc1,1:3)   ! èdêSÇå¥ì_Ç∆çlÇ¶ÇÈ
#endif

            pdd=ZERO;pdk=ZERO;pd_c=ZERO
            pdd=pd
            pdk=pd
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          
          ! ï¿êiÇÃçáóÕÇÃåvéZ
!          FFD(kdc1,:)=FFD(kdc1,:)+pd*rho0*S*(-face1(1:3))+Fvis(:) !+phy*rho0*S*(-face1(1:3))+rho0*dsign(pdy2(:),Vs(:))*Se(:)+Fvis(:) !+Fdrag(:)
          FFD1(:)=pd*rho0*S*(-face1(1:3))+Fvis(:)
!          Fvis_all(kdc1, 1:3) = Fvis(1:3)
          
 
        
          if(PL(2).gt.ZERO)then
              pdr(1:3)=pdr(1:3)+pdd*rho0*S*(-face1(1:3))  !âE
          elseif(PL(2).lt.ZERO)then
              pdl(1:3)=pdl(1:3)+pdd*rho0*S*(-face1(1:3))   !ç∂
          endif
                    
!          write(times,'(f6.3)') t
          986 format(2(e12.5,','),e12.5)          
          987 format(e12.5,',',e12.5)          
          988 format(4(e12.5,','),e12.5)          
          989 format(2(i2,','),2(e12.5,','),e12.5)          
          990 format(3(I3,','),i3)
          991 format(8(e12.5,','),e12.5)
          998 format(3(i3,','),i3) 
          999 format(2(e12.5,','),e12.5) 
              
        
#ifdef BUGWIN2

if(ffd(kdc1,3).ne.ZERO) then
continue
endif
          FFD2(kdc1,:)=FFD2(kdc1,:)+pd*rho0*S*(-face1(1:3))
#endif
 !         rot0=rot(kdc1,2)
          ! äeé≤âÒÇËÇÃà≥óÕÇ∆îSê´óÕÇ…ÇÊÇÈâÒì]ÉÇÅ[ÉÅÉìÉgÇÃåvéZ          
          call RotationalMoment(PL,-face1(1:3)*pd*rho0*S+Fvis(:),kdc1)   !+rho0*dsign(pdy2(:),Vs(:))*Se(:)+Fvis(:),kdc) !,drot) !+rho0*dsign(pdy2(:),Vs(:))*Se(:)+Fvis+Fdrag,kdc,drot)   

#ifdef DDDD          

#endif

!
!if(t.ge.3.0d0)  write(49,"(1h ,e12.5,',',i9,',',i4,11(',',e12.5))") t,nn,kk,pd,s,face1,rho0,pl,rot(kdc,2)-rot0
!        write(*,*) 'kdc,rot=',kdc,rot(1:3)
#ifdef BUGWIN2
if(kdc1.eq.1.and.abs(rot(kdc1,1)).ne.ZERO) then
continue
endif        
        continue
#endif

    end subroutine
    
    
subroutine cal_Vw(nn,Pin,Vout)
use variables, only: ZERO
use arrays,only:u,mn,x,dx,in,jn,kn
!use variables,only:nn
!        use UNVNWN
        implicit none
        integer,intent(in)::nn
        doubleprecision,dimension(3),intent(in)::Pin
        doubleprecision,dimension(3),intent(out)::Vout
        doubleprecision::um,up,vm,vp,wm,wp
        integer::i,k,j
        um=ZERO;up=ZERO;vm=ZERO;vp=ZERO;wm=ZERO;wp=ZERO
        i=in(nn);j=jn(nn);k=kn(nn)
        !if(ax(nn).ne.ZERO.and.ax(mn(i+1,k,j)).ne.ZERO) then
          um=u(0,nn)
          up=u(0,mn(i+1,k,j))
        !elseif(ax(nn).eq.ZERO.and.ax(mn(i+1,k,j)).ne.ZERO) then
        !  um=u(mn(i+1,k,j))
        !  up=u(mn(i+1,k,j))
        !elseif(ax(nn).ne.ZERO.and.ax(mn(i+1,k,j)).eq.ZERO) then
        !  um=u(nn)
        !  up=u(nn)
        !endif
        !
        !if(ay(nn).ne.ZERO.and.ay(mn(i,k,j+1)).ne.ZERO) then
          vm=u(1,nn)
          vp=u(1,mn(i,k,j+1))
        !elseif(ay(nn).eq.ZERO.and.ay(mn(i,k,j+1)).ne.ZERO) then
        !  vm=v(mn(i,k,j+1))
        !  vp=v(mn(i,k,j+1))
        !elseif(ay(nn).ne.ZERO.and.ay(mn(i,k,j+1)).eq.ZERO) then
        !  vm=v(nn)
        !  vp=v(nn)
        !endif
        !
        !if(az(nn).ne.ZERO.and.az(mn(i,k+1,j)).ne.ZERO) then
          wm=u(2,nn)
          wp=u(2,mn(i,k+1,j))
        !elseif(az(nn).eq.ZERO.and.az(mn(i,k+1,j)).ne.ZERO) then
        !  wm=w(mn(i,k+1,j))
        !  wp=w(mn(i,k+1,j))
        !elseif(az(nn).ne.ZERO.and.az(mn(i,k+1,j)).eq.ZERO) then
        !  wm=w(nn)
        !  wp=w(nn)
        !endif
        
        Vout(1)=((x(0,i+1)-Pin(1))*um+(Pin(1)-x(0,i))*up)/dx(0,i)
        Vout(2)=((x(1,j+1)-Pin(2))*vm+(Pin(2)-x(1,j))*vp)/dx(1,j)
        Vout(3)=((x(2,k+1)-Pin(3))*wm+(Pin(3)-x(2,k))*wp)/dx(2,k)

!        Vout(1)=((x(0,i+1)-Pin(1))*u(nn)+(Pin(1)-x(0,i))*u(mn(i+1,k,j)))/dx(0,i)
!        Vout(2)=((x(1,j+1)-Pin(2))*v(nn)+(Pin(2)-x(1,j))*v(mn(i,k,j+1)))/dx(1,j)
!        Vout(3)=((x(2,k+1)-Pin(3))*w(nn)+(Pin(3)-x(2,k))*w(mn(i,k+1,j)))/dx(2,k)
        !write(*,*) 'Vout from cal_Vw'
        !write(*,*) Vout
end subroutine

    
    