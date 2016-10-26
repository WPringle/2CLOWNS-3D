module local_module
real*8,allocatable,dimension(:) :: xstart,beta
integer::flat            ! flat=1 セル内勾配水平、＝0　セル内勾配が水平でない
type adj
real*8  :: hiniadj,xs,xe,ys,ye
integer :: type
endtype adj
type(adj),allocatable,dimension(:) :: adjinfo
integer :: slopenum, adjnum
contains
subroutine read_local_data
end subroutine
end module local_module
    
program make_chikei_polygon
USE variables,only:ndri,refval,order,fmt1,inputxyz,is,js,ks,ie,je,ke,openbound,title,input,inputf,inns,inne,inne2,fg_xyz
USE arrays,only:OBJECTshape,xi,y,z,nf3,nfb3,ax,ay,az,fb3,f3,x,in,jn,kn,mn,u,p,f,fb,a,ivt,nff
USE local_module
implicit none
integer :: i,j,k,nx,inn,nn,icc
logical :: exist
integer :: n1,n2,na1,na2,dxnum,dynum,dznum,ii,iimax ,ivt0,ist

real*8 :: hini
real*8,allocatable,dimension(:) :: xx,yy,xxst,yyst,zzst,dxx,dyy,dzz
character(len=255)::outputshape
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%              Read the input text file               %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef BUGWIN
!open(5,file='cwp_depth_1350-01.txt',status='old')
!open(5,file='cwp_depth_0450-04.txt',status='old')
open(5,file='cwp_tohoku3D_GC_zl.txt',status='old')
#endif
read(5,*) ndri
allocate(OBJECTshape(max(1,ndri))) ; exist = .false.
do i=1,max(1,ndri)
    read(5,255) OBJECTshape(i)
    !Check if odata exists
    n1 = len_trim(OBJECTshape(i))
    if (OBJECTshape(i)(n1-4:n1) == 'odata') then
        exist = .true.
    endif
    if (i.gt.1.and.i.eq.ndri) read(5,*) refval
enddo
read(5,*) order,fmt1
    write(6,*) 'order to read:',order
    write(6,*) 'format to read:',fmt1
read(5,*) dxnum,dynum,dznum
if (dxnum.gt.0) then
    allocate(dxx(dxnum),dyy(dynum),dzz(dznum),&
        xxst(dxnum+1),yyst(dynum+1),zzst(dznum+1))
    read(5,*) (xxst(ii),ii=1,dxnum+1)
    read(5,*) (yyst(ii),ii=1,dynum+1)
    read(5,*) (zzst(ii),ii=1,dznum+1)
    read(5,*) (dxx(ii),ii=1,dxnum)
    read(5,*) (dyy(ii),ii=1,dynum)
    read(5,*) (dzz(ii),ii=1,dznum)
else
    read(5,255) inputxyz
    if(dxnum.ne.0) then ! 1)  dxnum == 0 use inputxyz only  2) dxnum < 0 use inputxyz maindata fdata
    read(5,255) input
    endif
endif
read(5,255) title
read(5,*) hini
read(5,*) openbound
read(5,*) adjnum
if (adjnum.gt.0) then
    allocate(adjinfo(adjnum))
    read(5,*) adjinfo
endif
read(5,*) slopenum
if (slopenum.gt.0) then
    allocate(xstart(slopenum),beta(slopenum))
    read(5,*) xstart
    read(5,*) beta
endif
read(5,*,iostat=ist) flat        
if(ist.ne.0) then
    write(6,*) 'Please input flat or not'
    write(6,*) ' flat=1 セル内勾配水平、＝0　セル内勾配が水平でない'
    stop
endif
#ifdef BUGWIN
close(5)
#endif
255 format(A255)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (exist) then
    !open the different polygons 
    call read_objectshape
    ivt0 = ivt(1)
endif    
       ! dxnum>0 create grid data
       ! dxnum<0 use full data( xyz main fdata)
       ! dxnum=0 use xyz data  
if(dxnum.gt.0) then
    !Create grid data
    call set_grid(dxx,dyy,dzz,dxnum,dynum,dznum,xxst,yyst,zzst,exist)
!
else
    !Use existing grid data

      if(dxnum.lt.0) then          
        call owari(inputxyz,n1)
        if(inputxyz(n1:n1).eq.'z') then
            call read_xyz3
            call read_main3
        else
            call read_data
            allocate(xi(1:ie+1),y(1:je+1),z(1:ke+1))
            xi(is-1:ie+1)=x(0,is-1:ie+1);y(js-1:je+1)=x(1,js-1:je+1);z(ks-1:ke+1)=x(2,ks-1:ke+1)
            deallocate(x,mn,u,f,p)
        endif
      
        call read_fdata3
                    deallocate(in,jn,kn)
      else
        call read_xyz3

      endif
endif
!Allocate the other arrays

allocate(f3(is-2:ie+2,1:ke+2,js-2:je+2));f3=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 地形ポリゴンの上部のセルを空セルに
call cal_fluid_area(z(ke),exist)  !,slopenum)  !,xstart,beta)

if(exist) then
    if(ivt(1).lt.ivt0) then
        n1 = len_trim(OBJECTshape(1))
        if (OBJECTshape(1)(n1-4:n1) == 'odata') then
            outputshape=OBJECTshape(1)(1:n1-6)//"_cut.odata"
            call output_odata_bin(outputshape,1)
        endif
    endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Call the setting of the f and boundary openness
call set_chikei(hini)  !,adjnum,adjinfo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

iimax=0
do i=is,ie
do j=js,je
do k=ks,ke
if(nf3(i,k,j).ge.0) iimax=max(i,iimax)
enddo
enddo
enddo
write(6,*) 'iimax,ie=',iimax,ie
ie=iimax+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call owari(title,n2)
if(n2.lt.3) then
    call owari(OBJECTshape(1),n2);call hajime(OBJECTshape(1),n1);
    if (ks == ke-1) then
        if (exist) title=OBJECTshape(1)(n1:n2-6)//'2D.t0.000'
        if (.not.exist) title=OBJECTshape(1)(n1:n2-4)//'2D.t0.000'
    else
        if (exist) title=OBJECTshape(1)(n1:n2-6)//'.t0.000'
        if (.not.exist) title=OBJECTshape(1)(n1:n2-4)//'.t0.000'
    endif
endif
call owari(title,nx)
write(6,*) 'bi-data write-start [',title(1:nx),'] nx=',nx
allocate(x(0:2,is-2:max(ie,je,ke)+1));
x(0,is-2:ie+1)=xi(is-2:ie+1);x(1,js-2:je+1)=y(js-2:je+1);x(2,ks-2:ke+1)=z(ks-2:ke+1)

!deallocate(xi,y,z)
!inn=0
!do j=js-1,je-1;do k=ks-1,ke-1;do i=is-1,ie-1
!if(nf(i,k,j).eq.-1.and.nfb(i,k,j).eq.0) then
!if(nf(i-1,k,j).ge.0.or.nf(i,k-1,j).ge.0.or.nf(i,k,j-1).ge.0 ) then
!    nfb(i,k,j)=200;inn=inn+1
!    cycle
!endif
!if(nfb(i,k,j).eq.200) then
!    continue
!endif
!endif
!enddo;enddo;enddo
!write(6,*) 'count of nfb=200 =',inn

call output_xyzn(1)

inne2=ie*je*ke;inn=0;icc=0
allocate(in(inne2),jn(inne2),kn(inne2),mn(ie+1,ke+1,je+1))
do j=js-1,je;do k=ks-1,ke;do i=is-1,ie
      if(nf3(i,k,j).ge.0.or.&
         (nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1&
          .and.nfb3(i,k,j).lt.200)) then
        inn=inn+1;mn(i,k,j)=inn;in(inn)=i;jn(inn)=j;kn(inn)=k
        if(nf3(i,k,j).eq.-1) icc=icc+1
          else
         mn(i,k,j)=0
      end if
enddo;enddo;enddo
inns=1;inne=inn
write(6,*) 'inne=',inne
write(6,*) 'count of nf=-1 ',icc
do j=js-1,je;do k=ks-1,ke;do i=is-1,ie
        if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.200) then
          inn=inn+1;mn(i,k,j)=inn;in(inn)=i;jn(inn)=j;kn(inn)=k
        end if
enddo;enddo;enddo

      inne2=inn
write(6,*) 'inne2=',inne2,'inne2-inne=',inne2-inne
fg_xyz=1
allocate(nff(0:inne2+1))
      do j=js-1,je
          do k=ks-1,ke;
              do i=is-1,ie
                  nn=mn(i,k,j)
                  if(nn.eq.0) then 
                  if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).eq.0)  cycle
                  if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).eq.-1)  then
                      mn(i,k,j)=inne2+1;cycle
                  else
                      continue
                  endif    
                  endif
                  nff(nn)%f=nf3(i,k,j)
                  nff(nn)%b=nfb3(i,k,j)
                  nff(nn)%fp=nf3(i,k,j)
                  nff(nn)%bp=nfb3(i,k,j)
          enddo
enddo
enddo
      nff(0)%f=-1;      nff(0)%b=0;
      nff(0)%fp=-1;      nff(0)%bp=0;   
      nff(inne2+1)%f=-1;      nff(inne2+1)%b=-1;
      nff(inne2+1)%fp=-1;      nff(inne2+1)%bp=-1;   
      
allocate(f(1:inne2),fb(1:inne2))
do nn=1,inne;
    f(nn)=f3(in(nn),kn(nn),jn(nn))
    fb(nn)=fb3(in(nn),kn(nn),jn(nn))
;enddo
deallocate(f3,fb3)
allocate(u(0:2,1:inne2),p(1:inne2));u=0.0d0;p=0.0d0
call output_data
deallocate(u,p)
allocate(a(0:2,1:inne2));a=0.0d0
do nn=1,inne2
    a(0,nn)=ax(in(nn),kn(nn),jn(nn))
    a(1,nn)=ay(in(nn),kn(nn),jn(nn))
    a(2,nn)=az(in(nn),kn(nn),jn(nn))
;enddo
call output_fdata
call write_godata
contains
    


#ifdef DDDDD
subroutine hyoko(zh)
doubleprecision::zh,zv(0:2)
doubleprecision,dimension(0:1,0:1)::abi,ab0
integer::ii,kdc
zv(2)=zmax+1.0d0
do ii=1,ipl(kdc)
a1(:)=vt(1,npo(kdc,ii,2),1:3)-vt(1,npo(kdc,ii,1),1:3)
b1(:)=cc(:)-vt(1,npo(kdc,ii,1),1:2)
c1=a1(0)*b1(1)-a1(1)*b1(0)
if(c1.lt.0.0d0) cycle

a2(:)=vt(1,npo(kdc,ii,3),1:3)-vt(1,npo(kdc,ii,2),1:3)
b2(:)=cc(:)-vt(1,npo(kdc,ii,2),1:2)
c2=a2(0)*b2(1)-a2(1)*b2(0)
if(c2.lt.0.0d0) cycle

if(ipo(kdc,ii).eq.3) then

a3(:)=vt(1,npo(kdc,ii,1),1:3)-vt(1,npo(kdc,ii,3),1:3)
b3(:)=cc(:)-vt(1,npo(kdc,ii,3),1:2)
c3=a3(0)*b3(1)-a3(1)*b3(0)
if(c3.lt.0.0d0) cycle

ab0(0,0)=a1(0);ab0(0,1)=(-1.0d0)*a3(0)
ab0(1,0)=a1(1);ab0(1,1)=(-1.0d0)*a3(1)
det=a1(0)*(-1.0d0)*a3(1)-a1(1)*(-1.0d0)*a3(0)
abi(0,0)=(-1.0d0)*a3(1);abi(0,1)=a3(0)
abi(1,0)=(-1.0d0)*a1(1);abi(1,1)=a1(0)
abi=abi/det
alp(:)=matmul(abi,b1)
zv(:)=alp(0)*a1-alp(1)*a3+vt(1,npo(kdc,ii,1),1:3)

else

a3(:)=vt(1,npo(kdc,ii,4),1:3)-vt(1,npo(kdc,ii,3),1:3)
b3(:)=cc(:)-vt(1,npo(kdc,ii,3),1:2)
c3=a3(0)*b3(1)-a3(1)*b3(0)
if(c3.lt.0.0d0) cycle

a4(:)=vt(1,npo(kdc,ii,1),1:3)-vt(1,npo(kdc,ii,4),1:3)
b4(:)=cc(:)-vt(1,npo(kdc,ii,4),1:2)
c4=a4(0)*b4(1)-a4(1)*b4(0)
if(c4.lt.0.0d0) cycle


ab0(0,0)=a1(0);ab0(0,1)=(-1.0d0)*a4(0)
ab0(1,0)=a1(1);ab0(1,1)=(-1.0d0)*a4(1)
det=a1(0)*(-1.0d0)*a4(1)-a1(1)*(-1.0d0)*a4(0)
abi(0,0)=(-1.0d0)*a4(1);abi(0,1)=a4(0)
abi(1,0)=(-1.0d0)*a1(1);abi(1,1)=a1(0)
abi=abi/det
alp(:)=matmul(abi,b1)
zv(:)=alp(0)*a1-alp(1)*a4+vt(1,npo(kdc,ii,1),1:3)
endif




!if(c1.ge.0.0d0.and.c2.ge.0.0d0.and.c3.ge.0.0d0.and.c4.ge.0.0d0) then

continue
!endif
enddo
zh=zv(2)
end subroutine


subroutine ploygon2grid
do i=is,ie-1
do j=js,je-1
    
if(i.eq.82.and.j.eq.73) then
    continue
endif 
if(i.eq.4.and.j.eq.56) then
    continue
endif    
    
    

cc(:)=(/xi(i),y(j)/)
call hyoko(zh(0))
cc(:)=(/xi(i+1),y(j)/)
call hyoko(zh(1))
cc(:)=(/xi(i),y(j+1)/)
call hyoko(zh(2))
cc(:)=(/xi(i+1),y(j+1)/)
call hyoko(zh(3))
zhmin=minval(zh)

do ii=1,ivt(1)
cc(:)=vt(1,ii,1:2)

a1(:)=(/xi(i+1),y(j),0.0d0/)-(/xi(i),y(j),0.0d0/)
b1(:)=cc(:)-(/xi(i),y(j)/)
c1=a1(0)*b1(1)-a1(1)*b1(0)
if(c1.le.0.0d0) cycle

a2(:)=(/xi(i+1),y(j+1),0.0d0/)-(/xi(i+1),y(j),0.0d0/)
b2(:)=cc(:)-(/xi(i+1),y(j)/)
c2=a2(0)*b2(1)-a2(1)*b2(0)
if(c2.le.0.0d0) cycle

a3(:)=(/xi(i),y(j+1),0.0d0/)-(/xi(i+1),y(j+1),0.0d0/)
b3(:)=cc(:)-(/xi(i+1),y(j+1)/)
c3=a3(0)*b3(1)-a3(1)*b3(0)
if(c3.le.0.0d0) cycle

a4(:)=(/xi(i),y(j),0.0d0/)-(/xi(i),y(j+1),0.0d0/)
b4(:)=cc(:)-(/xi(i),y(j+1)/)
c4=a4(0)*b4(1)-a4(1)*b4(0)
if(c4.le.0.0d0) cycle

!if(c1.ge.0.0d0.and.c2.ge.0.0d0.and.c3.ge.0.0d0.and.c4.ge.0.0d0) then
zhmin=min(vt(1,ii,3),zhmin)
!endif
enddo


do k=ks,ke-1
if(z(k+1).gt.zhmin) then
if(nf(i,k,j).eq.-1) then
nf(i,k,j)=0
nfb(i,k,j)=0
fb(i,k,j)=1.0d0
endif
endif
enddo
enddo
enddo
end subroutine ploygon2grid
#endif  
end program
#ifdef DDDDD
subroutine cross(A,B,C)
        implicit none
        doubleprecision,dimension(3),intent(in)::A,B
        doubleprecision,dimension(3),intent(out)::C
        doubleprecision::CL
        
        C(:)=0.0d0
        C(1)=A(2)*B(3)-A(3)*B(2)
        C(2)=A(3)*B(1)-A(1)*B(3)
        C(3)=A(1)*B(2)-A(2)*B(1)      
        
        CL=sqrt(sum(C(:)**2))
        if(CL.ne.0.0d0) then
          if(abs(C(1)/CL).lt.1.0d-10) C(1)=0.0d0
          if(abs(C(2)/CL).lt.1.0d-10) C(2)=0.0d0
          if(abs(C(3)/CL).lt.1.0d-10) C(3)=0.0d0
        endif
        
end subroutine
#endif   


subroutine write_godata
use variables,only:title,ks,ke,js,je,is,ie
use arrays,only:z,y,xi,nf3,fb,mn
implicit none
integer::mmd,ndmx,nx,k,j,i
doubleprecision::z1,z2

ndmx=4
mmd=0
do i=is,ie-1
    do j=js,je-1
        do k=ks,ke-1
            if(nf3(i,k,j).eq.-1) cycle
            if(nf3(i-1,k,j).eq.-1.and.i.gt.is) mmd=mmd+1
            if(nf3(i+1,k,j).eq.-1.and.i.lt.ie-1) mmd=mmd+1
            if(nf3(i,k-1,j).eq.-1) mmd=mmd+1
            if(nf3(i,k+1,j).eq.-1.and.k.lt.ke-1) mmd=mmd+1
            if(nf3(i,k,j-1).eq.-1.and.j.gt.js) mmd=mmd+1
            if(nf3(i,k,j+1).eq.-1.and.j.lt.je-1) mmd=mmd+1
        enddo
    enddo
enddo
call owari(title,nx)
mmd=0
call flatplain(mmd)

OPEN(9,FILE=title(1:nx)//'.3d.godata',form='unformatted',STATUS= 'UNKNOWN')
write(9) mmd,ndmx 
call flatplain(mmd)
CLOSE(9)        
       write(6,*) 'mmd,ndmx=',mmd,ndmx
contains
subroutine flatplain(iflag)
integer,intent(in)::iflag
mmd=iflag
do i=is,ie-1
    do j=js,je-1
        do k=ks,ke-1
            if(nf3(i,k,j).eq.-1) cycle
            
            if(nf3(i-1,k,j).eq.-1.and.i.gt.is) then
            if(iflag.eq.0)    mmd=mmd+1
            if(iflag.ne.0)    write(9) ndmx,xi(i),y(j),z(k),xi(i),y(j+1),z(k),xi(i),y(j+1),z(k+1),xi(i),y(j),z(k+1)
            endif
            if(nf3(i,k,j-1).eq.-1.and.j.gt.js) then
            if(iflag.eq.0)    mmd=mmd+1
            if(iflag.ne.0)    write(9) ndmx,xi(i),y(j),z(k),xi(i),y(j),z(k+1),xi(i+1),y(j),z(k+1),xi(i+1),y(j),z(k)
            endif
            z1=z(k)+(z(k+1)-z(k))*(1.0-fb(mn(i,k,j)))
            continue
            if(nf3(i,k-1,j).eq.-1) then 
            if(iflag.eq.0)    mmd=mmd+1
            if(iflag.ne.0)    write(9) ndmx,xi(i),y(j),z1,xi(i+1),y(j),z1,xi(i+1),y(j+1),z1,xi(i),y(j+1),z1
            endif
            if(nf3(i+1,k,j).eq.-1.and.i.lt.ie-1) then
            if(iflag.eq.0)    mmd=mmd+1
            if(iflag.ne.0)    write(9) ndmx,xi(i+1),y(j),z(k),xi(i+1),y(j),z(k+1),xi(i+1),y(j+1),z(k+1),xi(i+1),y(j+1),z(k)
            endif
            if(nf3(i-1,k,j).ge.0) then
                z2=z(k)+(z(k+1)-z(k))*(1.0-fb(mn(i-1,k,j)))                
                if(z1-z2.gt.1.0d-2) then 
                    if(iflag.ne.0)    write(9) ndmx,xi(i),y(j),z2,xi(i),y(j),z1,xi(i),y(j+1),z1,xi(i),y(j+1),z2
                    if(iflag.eq.0)    mmd=mmd+1
                endif
            endif
            if(nf3(i+1,k,j).ge.0) then
                z2=z(k)+(z(k+1)-z(k))*(1.0-fb(mn(i+1,k,j)))                
                if(z1-z2.gt.1.0d-2) then 
                    if(iflag.ne.0)    write(9) ndmx,xi(i+1),y(j),z2,xi(i+1),y(j+1),z2,xi(i+1),y(j+1),z1,xi(i+1),y(j),z1
                    if(iflag.eq.0)    mmd=mmd+1
                endif
            endif
                
            if(nf3(i,k,j+1).eq.-1.and.j.lt.je-1) then
            if(iflag.eq.0)    mmd=mmd+1
            if(iflag.ne.0)    write(9) ndmx,xi(i),y(j+1),z(k),xi(i+1),y(j+1),z(k),xi(i+1),y(j+1),z(k+1),xi(i),y(j+1),z(k+1)
            endif
            if(nf3(i,k,j-1).ge.0) then
                z2=z(k)+(z(k+1)-z(k))*(1.0d0-fb(mn(i,k,j-1)))                
                if(z1-z2.gt.1.0d-2) then 
                    if(iflag.ne.0)   write(9) ndmx,xi(i),y(j),z2,xi(i+1),y(j),z2,xi(i+1),y(j),z1,xi(i),y(j),z1
                    if(iflag.eq.0)    mmd=mmd+1
                endif
            endif
            if(nf3(i,k,j+1).ge.0) then
                z2=z(k)+(z(k+1)-z(k))*(1.0d0-fb(mn(i,k,j+1)))                
                if(z1-z2.gt.1.0d-2) then 
                    if(iflag.ne.0)   write(9) ndmx,xi(i),y(j+1),z2,xi(i),y(j+1),z1,xi(i+1),y(j+1),z1,xi(i+1),y(j+1),z2
                    if(iflag.eq.0)    mmd=mmd+1
                endif
            endif
            if(nf3(i,k+1,j).eq.-1.and.k.lt.ke-1) then 
            if(iflag.eq.0)    mmd=mmd+1
            if(iflag.ne.0)    write(9) ndmx,xi(i),y(j),z(k+1),xi(i),y(j+1),z(k+1),xi(i+1),y(j+1),z(k+1),xi(i+1),y(j),z(k+1)
            endif
        enddo
    enddo
enddo  
end subroutine flatplain
end subroutine write_godata     
