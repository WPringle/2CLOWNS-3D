subroutine set_grid(dxx,dyy,dzz,dxnum,dynum,dznum,xxst,yyst,zzst,exist)
USE arrays,only:vt,ivt,xi,y,z,OBJECTshape
USE variables,only:ndri,is,js,ks,ie,je,ke,order,fmt1
use interface_list,only:round
implicit none       
integer :: i,j,k,unit,ierr,itemp
logical :: exist
integer :: dxnum,dynum,dznum,n,gn,ix,iix,dec,izh
integer,dimension(4) :: gnum
real*8,dimension(12) :: xxst,yyst,zzst,dxx,dyy,dzz,galpha
real*8 :: mindxx, zhmax=-100000.0d0, zhmin= 100000.0d0
real*8,allocatable::zh(:)
!計算領域入力
if (exist) then
if(xxst(dxnum+1).eq.-10000d0) xxst(dxnum+1)=max(maxval(vt(1,1:ivt(1),1)),maxval(vt(ndri,1:ivt(ndri),1)))    !ポリゴンデータを使う場合
if(yyst(dynum+1).eq.-10000d0) yyst(dynum+1)=max(maxval(vt(1,1:ivt(1),2)),maxval(vt(ndri,1:ivt(ndri),2))) 
if(zzst(dznum+1).eq.-10000d0) zzst(dznum+1)=max(maxval(vt(1,1:ivt(1),3)),maxval(vt(ndri,1:ivt(ndri),3))) 
if(xxst(1).eq.-10000d0) xxst(1)=min(minval(vt(1,1:ivt(1),1)),minval(vt(ndri,1:ivt(ndri),1)))    !ポリゴンデータを使う場合
if(yyst(1).eq.-10000d0) yyst(1)=min(minval(vt(1,1:ivt(1),2)),minval(vt(ndri,1:ivt(ndri),2))) 
if(zzst(1).eq.-10000d0) zzst(1)=min(minval(vt(1,1:ivt(1),3)),minval(vt(ndri,1:ivt(ndri),3))) 
else
    if (xxst(dxnum+1).eq.-10000d0.or.xxst(1).eq.-10000d0.or. &
        yyst(dynum+1).eq.-10000d0.or.xxst(1).eq.-10000d0) then
        write(6,*) 'We cannot use -10000d0 to set the domain size here'
        stop 'since the bounds are unknown'
    endif
    if (zzst(dznum+1).eq.-10000d0.or.zzst(1).eq.-10000d0) then
        
         if (fmt1(1:1).eq.'*') then
             izh=int((xxst(dxnum+1)-xxst(dxnum))/dxx(1))
         else
             izh=10
         endif
             allocate(zh(1:izh)  )       
        
        
        do i=1,max(1,ndri)
        write(6,*) OBJECTshape(i) 
        open(file=OBJECTshape(i) ,status='old',newunit=unit)
        do
        if (fmt1(1:1).eq.'*') then
          read(unit,*,end=200,iostat=ierr) zh(1:izh)           
        else
        read(unit,'(10f8.2)',end=200,iostat=ierr) zh(1:izh) 
        endif
        if (order(3:3).eq.'N') zh=-1.0d0*zh
        zhmax=max(zhmax,maxval(zh(1:izh)))
        zhmin=min(zhmin,minval(zh(1:izh)))
        enddo
200    continue 
        close(unit)
201     continue
        if(ierr.gt.0.or.ierr.eq.-2) then
            write(6,*) 'read error ierr=',ierr
            stop
        endif
        enddo
        if (zzst(dznum+1).eq.-10000d0) zzst(dznum+1)=zhmax
        if (zzst(1).eq.-10000d0) zzst(1)=zhmin
        continue
    endif        
endif
write(6,*) 'xmax,xmin=',xxst(dxnum+1),xxst(1)
write(6,*) 'ymax,ymin=',yyst(dynum+1),yyst(1)
write(6,*) 'zmax,zmin=',zzst(dznum+1),zzst(1)

! Determine decimal place for rounding. Determined by 
! rounding to one decimal place more than the min dx
mindxx = minval(dxx(1:dxnum),dxx>0)
do dec = 0,10
    if (mindxx/10d0**(3-dec).gt.1d0) exit
enddo
!格子間隔設定
is=3;js=3;ks=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ie=is            !X方向メッシュ間隔
gn=0;galpha=0.0d0
do i=1,dxnum
    if (xxst(i).eq.-10000d0) cycle
    if (xxst(i+1).ne.-10000d0) then
        !Static cell size between locations
        ie=ie+nint((xxst(i+1)-xxst(i))/dxx(i))          
    else
        !Gradual increase or decrease in 
        !cell size between the locations 
        gn=gn+1; ix = 0
        !In somecases we need to shift backwards one 
        if (dxx(i).eq.-10000d0) ix = -1
        !Call the subroutine gradinc to evaluate n and alpha 
        call gradinc(gnum(gn),galpha(gn),dxx(i+ix),dxx(i+2+ix),xxst(i),xxst(i+2))
        ie=ie+gnum(gn)  
    endif
enddo
allocate(xi(is-2:ie+2))
xi(is)=xxst(1)
xi(is-2)=xi(is)-dxx(1)*2.0d0;xi(is-1)=xi(is)-dxx(1)
k=0;gn=1
do i=is+1,ie
    do j=dxnum,1,-1
        if (xxst(j).eq.-10000d0) cycle
        if (xi(i-1).ge.xxst(j)) then
            exit
        endif
    enddo
    if (xxst(j+1).ne.-10000d0) then
        !Static cell size between locations
        xi(i)=round(xi(i-1)+dxx(j),dec)
    else
        !Gradual increase or decrease in 
        !cell size between the locations
        k=k+1
        if (k.gt.gnum(gn)) then
            k=1;gn=gn+1
        endif
        if (k.eq.gnum(gn)) then
            xi(i)=xxst(j+2)
        else
            ix = 0
            !In somecases we need to shift backwards one 
            if (dxx(j).eq.-10000d0) ix = -1
            !Find the new x value
            xi(i)=round(xi(i-1)+dxx(j+ix)*galpha(gn)**(k-1),dec)
        endif
    endif
enddo
xi(ie+1)=round(xi(ie)+dxx(dxnum),dec)
write(6,*) 'coefficient value(s) for graduation x:',galpha(1:gn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine decimal place for rounding. Determined by 
! rounding to one decimal place more than the min dx
mindxx = minval(dyy(1:dxnum),dyy>0)
do dec = 0,10
    if (mindxx/10d0**(3-dec).gt.1d0) exit
enddo
je=js            !Y方向メッシュ間隔
gn=0;galpha=0.0d0
do i=1,dynum
    if (yyst(i).eq.-10000d0) cycle
    if (yyst(i+1).ne.-10000d0) then
        !Static cell size between locations
        je=je+nint((yyst(i+1)-yyst(i))/dyy(i))      
    else
        !Gradual increase or decrease in 
        !cell size between the locations 
        gn=gn+1; ix = 0
        !In somecases we need to shift backwards one 
        if (dyy(i).eq.-10000d0) ix = -1
        !Call the subroutine gradinc to evaluate n and alpha 
        call gradinc(gnum(gn),galpha(gn),dyy(i+ix),dyy(i+2+ix),yyst(i),yyst(i+2))
        je=je+gnum(gn)  
    endif
enddo
allocate(y(js-2:je+2))
y(js)=yyst(1)
y(js-2)=y(js)-dyy(1)*2.0d0;y(js-1)=y(js)-dyy(1)
k=0;gn=1
do i=js+1,je
    do j=dynum,1,-1
        if (yyst(j).eq.-10000d0) cycle
        if (y(i-1).ge.yyst(j)) then
            exit
        endif
    enddo
    if (yyst(j+1).ne.-10000d0) then
        !Static cell size between locations
        y(i)=round(y(i-1)+dyy(j),dec)
    else
        !Gradual increase or decrease in 
        !cell size between the locations
        k=k+1
        if (k.gt.gnum(gn)) then
            k=1;gn=gn+1
        endif
        if (k.eq.gnum(gn)) then
            y(i)=yyst(j+2)
        else
            ix = 0
            !In somecases we need to shift backwards one 
            if (dyy(j).eq.-10000d0) ix = -1
            !Find the new x value
            y(i)=round(y(i-1)+dyy(j+ix)*galpha(gn)**(k-1),dec)
        endif
    endif
enddo
y(je+1)=round(y(je)+dyy(dynum),dec)
write(6,*) 'coefficient value(s) for graduation y:',galpha(1:gn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine decimal place for rounding. Determined by 
! rounding to one decimal place more than the min dx
if (dzz(1).ne.-1d4) then
    mindxx = minval(dzz(1:dxnum),dzz>0)
    do dec = 0,10
        if (mindxx/10d0**(3-dec).gt.1d0) exit
    enddo
else
    dzz(1) = zzst(2)-zzst(1)
endif
ke=ks            !Z方向メッシュ間隔
gn=0;galpha=0.0d0
do i=1,dznum
    if (zzst(i).eq.-10000d0) cycle
    if (zzst(i+1).ne.-10000d0) then
        !Static cell size between locations
        ke=ke+nint((zzst(i+1)-zzst(i))/dzz(i))    
    else
        !Gradual increase or decrease in 
        !cell size between the locations 
        gn=gn+1; ix = 0
        !In somecases we need to shift backwards one 
        if (dzz(i).eq.-10000d0) ix = -1
        !Call the subroutine gradinc to evaluate n and alpha 
        call gradinc(gnum(gn),galpha(gn),dzz(i+ix),dzz(i+2+ix),zzst(i),zzst(i+2))
        ke=ke+gnum(gn)
        endif
enddo
allocate(z(1:ke+2))
z(3)=zzst(1)
z(1)=z(3)-dzz(1)*2.0d0;z(2)=z(3)-dzz(1)
k=0;gn=1
do i=ks+1,ke
    do j=dznum,1,-1
        if (zzst(j).eq.-10000d0) cycle
        if (z(i-1).ge.zzst(j)) then
            exit
        endif
    enddo
    if (zzst(j+1).ne.-10000d0) then
        !Static cell size between locations
        z(i)=round(z(i-1)+dzz(j),dec)
    else
        !Gradual increase or decrease in 
        !cell size between the locations
        k=k+1
        if (k.gt.gnum(gn)) then
            k=1;gn=gn+1
        endif
        if (k.eq.gnum(gn)) then
            z(i)=zzst(j+2)
        else
            ix = 0
            !In somecases we need to shift backwards one 
            if (dzz(j).eq.-10000d0) ix = -1
            !Find the new x value
            z(i)=round(z(i-1)+dzz(j+ix)*galpha(gn)**(k-1),dec)
        endif
    endif
enddo
z(ke+1)=round(z(ke)+dzz(dznum),dec)
write(6,*) 'coefficient value(s) for graduation z:',galpha(1:gn)
write(6,*)'ie,je,ke,z(ke+1)=', ie,je,ke,z(ke+1)
!
    contains
    
!real*8 function round2(num,dec)
!implicit none!
!real*8,intent(in) :: num
!integer,intent(in) :: dec
!real*8 :: A
!Round a number to the desired number of decimal places
! num :: the number
! dec :: number of decimal places for rounding
!A = 10d0**dec
!round2 = num  !float(nint(num*A))/A
!
!endfunction round2    
end subroutine
!    
subroutine gradinc(n,alpha,dzm,dzp,zm,zp)    
implicit none
integer,intent(out) :: n
real*8,intent(out) :: alpha
real*8,intent(in) :: dzm,dzp,zm,zp
real*8 :: sum,diff,eps,alphad,nd,modr
integer :: j
parameter(eps = 1d-6)
!We know relationship alpha^n = (dzp/dzm)
!and that the sum of the cell sizes should be 
!equal to (zp-zm) so we can solve the equations
!simultaneously for n and alpha..
alpha = 1d0 + ( dzp - dzm ) / ( zp - zm )
nd    = ( log(dzp) - log(dzm) ) / log(alpha)
!changing n to an integer
n     = int(nd)
modr  = mod(nd,real(n))
!when n is already an integer..
if (modr.lt.eps) return
!When n is not an integer recalculate alpha so that
!the sum becomes more exactly equal to (zzst(i+2)-zzst(i))
!Use Netwon raphson iterations...
3 alpha = alpha - ( dzm * ( 1 - alpha**n ) + ( 1 - alpha ) * ( zm - zp ) ) /  &
                  ( zp - zm - dzm * real(n) * alpha**(n-1) )
sum  = dzm * ( 1 - alpha**n ) / ( 1 - alpha ) 
diff = sum - (zp - zm )
if (abs(diff).gt.eps) goto 3
    return !return the found 'n' and 'alpha' 
endsubroutine
    
