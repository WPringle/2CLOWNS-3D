subroutine pritime2(date_time_start,date_time_old)
implicit none
integer::date_time_old(1:8),date_time_start(1:8)
integer::date_time(8)
real::s_int,s_start,sec(8)
CHARACTER*8 idate
CHARACTER*10 ymd,ctime,zone
CHARACTER*255::fm="(1h ,'*** ',a4,'-',a2,'-',a2,' ',a2,':',a2,':',a2,&
', elapse=',f8.2,'[m] interval=',f8.2,'[s]')"

CALL DATE_AND_TIME(ymd,ctime,zone,date_time)
if(date_time_start(1).ne.0) then
    sec(1:8)=float(date_time(1:8)-date_time_old(1:8))
    s_int=sec(5)*3600+sec(6)*60+sec(7)+sec(8)/1000.
    sec(1:8)=float(date_time(1:8)-date_time_start(1:8))
    s_start=sec(5)*60+sec(6)+sec(7)/60.
else 
    date_time_start(1:8)=date_time(1:8)
    s_int=0;s_start=0
endif


write(6,fm) ymd(1:4),ymd(5:6),ymd(7:8),ctime(1:2),ctime(3:4),ctime(5:6),s_start,s_int
   
date_time_old(1:8)=date_time(1:8)

end subroutine pritime2
!
subroutine kai2(a,b,c,x1,x2,icon)
doubleprecision,intent(in)::a,b,c
doubleprecision::x1,x2
integer::icon
icon=0

if(b**2-4.0d0*a*c.gt.0.0d0.and.a.ne.0.0d0) then
x1=(sqrt(b**2-4.0d0*a*c)-b)/2.0d0/a
x2=(sqrt(b**2-4.0d0*a*c)+b)/2.0d0/a*(-1.0d0)
else if(abs(b**2-4.0d0*a*c).lt.1.0d-8.and.a.ne.0.0d0) then
x1=(-b)/2.0d0/a
x2=(+b)/2.0d0/a*(-1.0d0)
else if(a.eq.0.0d0.and.b.ne.0.0d0) then
x1=-c/b
x2=-c/b
else
x1=-1.d4
x2=-1.d4
icon=1000
endif
end subroutine
!
doubleprecision function renzoku_fv(axp,axm,ayp,aym,azp,azm,up,um,vp,vm,wp,wm,dx,dy,dz,fb)
doubleprecision,intent(in)::axp,axm,ayp,aym,azp,azm,up,um,vp,vm,wp,wm,dx,dy,dz,fb
renzoku_fv=( (up*axp-um*axm)/dx+(vp*ayp-vm*aym)/dy+(wp*azp-wm*azm)/dz)/fb
end function renzoku_fv
!
doubleprecision function renzoku_fv3(axp,axm,ayp,aym,azp,azm,up,um,vp,vm,wp,wm,dx,dy,dz,fb,fbd,dt)
implicit none
doubleprecision,intent(in)::axp,axm,ayp,aym,azp,azm,up,um,vp,vm,wp,wm,dx,dy,dz,fb,fbd,dt
doubleprecision::fb2,renzoku_fv
fb2=(fb+fbd)*0.5d0
renzoku_fv3=renzoku_fv(axp,axm,ayp,aym,azp,azm,up,um,vp,vm,wp,wm,dx,dy,dz,fb2)+(fbd-fb)/dt/fb2
end function renzoku_fv3
!
doubleprecision function DISTA(vec)
doubleprecision,dimension(0:2),intent(in)::vec
DISTA=sqrt(sum(vec**2))
end function DISTA

integer function sweap(a,n)
implicit none
integer,intent(in)::n
doubleprecision::a(n,n+n)
doubleprecision::t,akk,aik
integer::k,i,j,m
m=n+n
do i=1,n
a(i,n+i)=1.0d0
enddo
sweap=0
do k=1,n
    do I=k,n
        if(a(i,k).ne.0) goto 25
    enddo
    sweap=-1
    return

25 if(i.ne.k) then 
    do j=k,m
         t=a(k,j)
         a(k,j)=a(i,j)
         a(i,j)=t
    enddo
    endif
    akk=a(k,k)
    a(k,k:m)=a(k,k:m)/akk
    do i=1,n
        if(i.eq.k) cycle
        aik=a(i,k)
        a(i,1:m)=a(i,1:m)-a(k,1:m)*aik
    enddo
enddo
a(1:n,1:n)=a(1:n,n+1:n+n)
end function
    
function crossV(A,B)
        implicit none
        doubleprecision,dimension(3)::crossV        
        doubleprecision,dimension(3),intent(in)::A,B
        doubleprecision::CL
        
        crossV(:)=0.0d0
        crossV(1)=A(2)*B(3)-A(3)*B(2)
        crossV(2)=A(3)*B(1)-A(1)*B(3)
        crossV(3)=A(1)*B(2)-A(2)*B(1)
        
 !       crossV(1)=A(3)*B(2)-A(2)*B(3)
 !       crossV(2)=A(1)*B(3)-A(3)*B(1)
 !       crossV(3)=A(2)*B(1)-A(1)*B(2)        
        
        CL=sqrt(sum(crossV(:)**2))
        if(CL.ne.0.0d0) then
          if(abs(crossV(1)/CL).lt.1.0d-15) crossV(1)=0.0d0
          if(abs(crossV(2)/CL).lt.1.0d-15) crossV(2)=0.0d0
          if(abs(crossV(3)/CL).lt.1.0d-15) crossV(3)=0.0d0
        endif
        
end function    
    
doubleprecision function crossL(A,B)
    implicit none
    doubleprecision,dimension(3),intent(in)::A,B
    doubleprecision,dimension(3)::C
    call cross(A,B,C)
    crossL=sqrt(sum(C(:)**2))
end function crossL
    
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
    
function centroid2(PS,n)
    implicit none
  INTERFACE 
    FUNCTION crossV(V1,V2)
      doubleprecision,DIMENSION(0:2):: crossV,V1,V2
    END FUNCTION
  END INTERFACE
    integer,intent(in)::n
    doubleprecision,dimension(0:2)::centroid2
    doubleprecision,dimension(1:n,0:2),intent(in)::PS
    doubleprecision,dimension(1:n-2,0:2)::G
    doubleprecision,dimension(1:n-2)::S
    doubleprecision::SSUM
    integer::i
    !
    centroid2=0.0d0
    SSUM=0.0d0
    do i=2,n-1
      G(i-1,:)=(PS(1,:)+PS(i,:)+PS(i+1,:))/3.0d0
      S(i-1)=sqrt(sum(crossV(PS(i,:)-PS(1,:),PS(i+1,:)-PS(1,:))**2))*0.5d0
 
    enddo
    do i=1,n-2
      centroid2(:)=centroid2(:)+S(i)*G(i,:)
      SSUM=SSUM+S(i)
    enddo
     if (SSUM.eq.0.0d0) then !Updated by Will Aug 19
        write(6,*) 'SSUM is zero in centroid' ,n
        centroid2=-10000.0d0
        return
     endif    
     centroid2(:)=centroid2(:)/SSUM

end function centroid2
    
function bibun(up,uc,um,dxp,dxc,dxm)
    implicit none
doubleprecision::bibun
doubleprecision,intent(in)::up,uc,um,dxp,dxc,dxm
    if((up-uc)*(uc-um).lt.0.0d0) then
      bibun=( (up-um) ) / ( dxc+0.5d0*(dxm+dxp) )
    else
      bibun=( (dxc+dxm)/(dxc+dxp)*(up-uc)+(dxc+dxp)/(dxc+dxm)*(uc-um) ) / ( dxc+0.5d0*(dxm+dxp) )
    endif
    end function

function center_diff2(Vp,vc,Vm,anup,anum,dt)
    implicit none
    doubleprecision::center_diff2
    TYPE velocity
    doubleprecision::v,x,am,ap
    integer::np,npb,nm,nmb,nnp,nnm
    END TYPE velocity
    type(velocity)::Vp,vc,Vm
    doubleprecision,intent(in)::anup,anum,dt
    doubleprecision::kzp,kzm
    doubleprecision::para=100000d0  !1.6d0

    kzp=anup*(Vp%v/Vc%v-1.0d0)/(Vp%x-Vc%x)/(Vp%x-Vm%x)*2.d0
    kzm=anum*(1.0d0-Vm%v/Vc%v)/(Vc%x-Vm%x)/(Vp%x-Vm%x)*2.d0

    center_diff2=Vc%v*( &
                       dsign(1.0d0,kzp)*min((Vp%x-Vc%x)*para/dt,abs(kzp)) &
                       -dsign(1.0d0,kzm)*min((Vc%x-Vm%x)*para/dt,abs(kzm)) &
                       )
end function
function seikika(vec)
    implicit none
    doubleprecision,dimension(0:2)::seikika   
    doubleprecision,dimension(0:2),intent(in)::vec
    doubleprecision::LL
    LL=sqrt(sum(vec(:)**2))
    seikika=vec/LL
    where(abs(seikika)<1.0d-12) seikika=0.0d0
end function seikika
    
function momentum(uc,ox,pm,pc,dx,f,dt)
    implicit none
    doubleprecision::momentum
    doubleprecision::uc,ox,pm,pc,dx,f,dt
    momentum=uc+dt*(OX+(Pm-Pc)/dx-f)
end function
doubleprecision function momentum2(uc,ox,ox2,pm,pc,dx,f,dt)
    implicit none
    interface
        function momentum(uc,ox,pm,pc,dx,f,dt)
            doubleprecision::momentum
            doubleprecision::uc,ox,pm,pc,dx,f,dt
        end function    
    end interface 
    doubleprecision,intent(in)::uc,ox,ox2,pm,pc,dx,f,dt
    doubleprecision::oxx
    oxx=(ox+ox2)*0.5d0
    momentum2=momentum(uc,oxx,pm,pc,dx,f,dt)
end function momentum2    
function momentum_d(uc,pm,pc,dx,dt)
    doubleprecision::momentum_d
    doubleprecision::uc,pm,pc,dx,dt
    momentum_d=uc+dt*(Pm-Pc)/dx
    end function
    
    
  doubleprecision function momentum_rho(uc,rc,rcn,ox,pm,pc,dx,f,dt)
    implicit none
    doubleprecision,intent(in)::uc,rc,rcn,ox,pm,pc,dx,f,dt
    momentum_rho=(uc*rc+dt*((Pm-Pc)/dx-f+OX))/rcn
end function


doubleprecision function momentum_rho2(uc,rc,rcn,ox,ox2,pm,pc,dx,f,dt)
   implicit none
   interface
     function momentum_rho(uc,rc,rcn,ox,pm,pc,dx,f,dt)
    doubleprecision::momentum_rho
        doubleprecision,intent(in)::uc,rc,rcn,ox,pm,pc,dx,f,dt
    end function  
    end interface 
 
    doubleprecision,intent(in)::uc,rc,rcn,ox,ox2,pm,pc,dx,f,dt
    doubleprecision::oxx
!    OXx=1.5d0*OX-0.5d0*ox2
    OXx=0.5d0*OX+0.5d0*ox2
    momentum_rho2=momentum_rho(uc,rc,rcn,oxx,pm,pc,dx,f,dt)
end function


doubleprecision function momentum_rho_d(uc,rcn,pm,pc,dx,dt)  !,rc,
    implicit none
    doubleprecision,intent(in)::uc,rcn,pm,pc,dx,dt !,rc
    momentum_rho_d=uc+(dt*(Pm-Pc)/dx)/rcn
end function
    
function nodeofLineandPlain(A,B,Face)
doubleprecision,dimension(0:2)::nodeofLineandPlain
doubleprecision,dimension(0:2),intent(in)::A,B
doubleprecision,dimension(0:3),intent(in)::face
doubleprecision,dimension(0:2)::BA
BA=B-A
nodeofLineandPlain=A - BA(0:2)*(  face(3) + dot_product( face(0:2),  A )) /dot_product( face(0:2), BA ) 
end function
    
    
    
    
    
    
    
    
function det(A)
    implicit none
    doubleprecision::det
    doubleprecision,dimension(2,2)::A
    det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
end function det  !determinant        

subroutine InverseMatrix3(A,C)
    implicit none
interface
function det(A1)
    implicit none
    doubleprecision::det
    doubleprecision,dimension(2,2)::A1
end function
end interface
    doubleprecision,dimension(3,3)::A,C
    doubleprecision,dimension(2,2)::B
    doubleprecision::det_a
    integer::i,j,ia,ja,ib,jb
    
    do i=1,3
      do j=1,3
        ib=0
        do ia=1,3
          if(ia.eq.i) cycle
          ib=ib+1;jb=0
          do ja=1,3
            if(ja.eq.j) cycle
            jb=jb+1
            b(ib,jb)=a(ia,ja)
          enddo
        enddo
        c(j,i)=(-1.0d0)**(j+i)*det(b)
#ifdef BUGWIN
        continue
#endif
      enddo
    enddo
    
    det_a=sum(c(1,1:3)*a(1:3,1))
    if(det_a.ne.0.0d0) then
      C=C/det_a
    else
      C=0.0d0
    endif
!    d=matmul(C,A)/det_a
#ifdef BUGWIN
    continue
#endif
 end subroutine InverseMatrix3
    
function cal_area(pp,ipp)
implicit none
  INTERFACE 
    FUNCTION crossV(V1,V2)
      doubleprecision,DIMENSION(0:2):: crossV,V1,V2
    END FUNCTION
  END INTERFACE
doubleprecision::cal_area
doubleprecision,dimension(0:2)::tmpV  !,tmpV1
doubleprecision,intent(in)::pp(1:ipp,0:2)
integer,intent(in)::ipp
integer::ib
tmpV=0.0d0
do ib=2,ipp-1
tmpV=crossV(pp(ib,0:2)-pp(1,0:2),pp(ib+1,0:2)-pp(1,0:2))+tmpV
enddo
cal_area=sqrt(sum(tmpv**2))
end function cal_area
