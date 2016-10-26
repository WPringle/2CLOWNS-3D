program add_teibou
use variables,only: TITLE,order,fmt,fmt1,inputxyz,input,inputf0
use arrays,only: one_dim_array
implicit none
type wall
   real*8 :: height,xstart,xend,ystart,yend
end type wall
type change2
   integer :: ii,jj
   real*8  :: newval
   integer :: ix,iy
end type change2
type(wall),allocatable :: bwval(:)
type(change2),allocatable :: valchange(:)
CHARACTER(LEN=255)::TEIBOU,MANNINGS
integer,allocatable,dimension(:,:) :: teibou_i
real*8,allocatable,dimension(:,:) :: teibou_h
integer :: nbw,nchange,i,nx,nx1,imod,nvchange
real*8  :: nval, oldval

#ifdef BUGWIN
open(5,file='ad_wall.txt',status='old')
#endif
!Grid data file
    read(5,255) inputxyz
!Teibou (wall) data file
    read(5,255) TEIBOU
!Mannings data file
    read(5,255) MANNINGS
    read(5,*) order,fmt,imod,fmt1 ! Direction to read data (xy) and read data format
    read(5,*) nbw                 ! Allocate number of teibous to add by hand
    if (nbw.ne.0) then
        allocate(bwval(1:nbw))
        do i=1,nbw
            read(5,*) bwval(i)
        enddo
    endif
    read(5,*) nchange
    if (nchange.ne.0) then
        allocate(valchange(1:nchange))
        do i=1,nchange
            read(5,*) valchange(i)
        enddo
    endif
    read(5,*) nvchange
    if (nvchange.ne.0) then
        read(5,*) oldval,nval
    endif
#ifdef BUGWIN
close(5)
#endif
255 format(A255)

! Read grid data
call owari(inputxyz,nx);call hajime(inputxyz,nx1)
if(inputxyz(nx:nx)=="z") then
    nx=nx-4
else if(inputxyz(nx:nx)=="n") then
    nx=nx-5    
endif
TITLE=inputxyz(nx1:nx)
input=inputxyz(nx1:nx)
inputf0=inputxyz(nx1:nx)//'.f0data'
!
! Read main data and fdata 
one_dim_array =.false.  ! éOéüå≥îzóÒÇ≈åüì¢ÅB
call read_data
call read_fdata3

! Read teibou (wall) data
if (teibou(1:3).ne.'NON') call read_teibou

! Alter teibou data based on teibou input by hand
if (nbw.ne.0) call teibou_hand

! Change main data and fdata based on teibou data
if (teibou(1:3).ne.'NON'.or.nbw.ne.0) call alter_w_teibou

call owari(title,nx)

! Add Mannings data
if (mannings(1:3).ne.'NON') call read_mannings
  !  if (teibou(1:3).ne.'NON'.or.nbw.ne.0) then
  !      title = title(1:nx)//'.fmtb' !;nx=nx+2
  !  else
  !      title = title(1:nx)//'.fm' !;nx=nx+2
  !  endif
  !else
  !  title = title(1:nx)//'.tb' !;nx=nx+2
  !endif

call owari(title,nx)

! Re-write fdata
write(6,*) 'fdata write-start',title(1:nx),'nx=',nx
call write_fdata3
!
contains
!
subroutine read_teibou
use variables,only:is,ie,js,je
implicit none
integer :: na1,iee,jee,iss,jss,n,m,i,j
real*8 :: dummy,xc,xb,yc,yb,dxx,dyy
real*8,allocatable,dimension(:,:) :: htemp
call owari(teibou,na1)
OPEN(9,FILE=teibou(1:na1),status='old')
    !read(9,*) iee, jee
    !read(9,*) dummy, dummy, dxx, dyy
    !iss=3;jss=3;jee=jee+jss;iee=iee+iss
    iss=is;jss=js;jee=je;iee=ie
    do i=1,10
        if (fmt(i:i).eq.'i') then
            !Read as integer and change after 
            allocate(teibou_i(iss:iee-1,jss:jee-1))
            exit
        endif
    enddo
    if (.not.allocated(teibou_i)) allocate(teibou_h(iss:iee-1,jss:jee-1))
    if (order(2:2).eq.'P') then
        do j=jss,jee-1
            if (fmt(1:1).eq.'*') then
                if (.not.allocated(teibou_i)) then
                    if (order(1:1).eq.'P') read(9,*) (teibou_h(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,*) (teibou_h(i,j),i=iee-1,iss,-1)
                else
                    if (order(1:1).eq.'P') read(9,*) (teibou_i(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,*) (teibou_i(i,j),i=iee-1,iss,-1)
                endif
            else
                if (.not.allocated(teibou_i)) then
                    if (order(1:1).eq.'P') read(9,fmt) (teibou_h(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,fmt) (teibou_h(i,j),i=iee-1,iss,-1)
                else
                    if (order(1:1).eq.'P') read(9,fmt) (teibou_i(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,fmt) (teibou_i(i,j),i=iee-1,iss,-1)                    
                endif
            endif  
        enddo
    elseif (order(2:2).eq.'N') then
        do j=jee-1,jss,-1
            if (fmt(1:1).eq.'*') then
                if (.not.allocated(teibou_i)) then
                    if (order(1:1).eq.'P') read(9,*) (teibou_h(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,*) (teibou_h(i,j),i=iee-1,iss,-1)
                else
                    if (order(1:1).eq.'P') read(9,*) (teibou_i(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,*) (teibou_i(i,j),i=iee-1,iss,-1) 
                endif
            else
                if (.not.allocated(teibou_i)) then
                    if (order(1:1).eq.'P') read(9,fmt) (teibou_h(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,fmt) (teibou_h(i,j),i=iee-1,iss,-1)
                else
                    if (order(1:1).eq.'P') read(9,fmt) (teibou_i(i,j),i=iss,iee-1)
                    if (order(1:1).eq.'N') read(9,fmt) (teibou_i(i,j),i=iee-1,iss,-1)                    
                endif
            endif           
        enddo
    endif
CLOSE(9)
!
endsubroutine read_teibou
!
subroutine read_mannings
use variables,only:is,ie,js,je,ks,ke
use arrays,only:az,xi,y
implicit none
integer :: na1,iee,jee,iss,jss,n,m,i,j
real*8 :: dummy,xc,xb,yc,yb,dxx,dyy
real*8,allocatable,dimension(:,:) :: aztemp
call owari(mannings,na1)
OPEN(9,FILE=mannings(1:na1),status='old',err=200)
    !read(9,*) iee, jee
    !read(9,*) dummy, dummy, dxx, dyy
    !iss=3;jss=3;jee=jee+jss;iee=iee+iss
    iss=is;jss=js;jee=je;iee=ie
    if (order(2:2).eq.'P') then
        do j=jss,jee-1
            if (fmt1(1:1).eq.'*') then
                if (order(1:1).eq.'P') read(9,*) (az(i,ks,j),i=iss,iee-1)
                if (order(1:1).eq.'N') read(9,*) (az(i,ks,j),i=iee-1,iss,-1)
            else
                if (order(1:1).eq.'P') read(9,fmt1) (az(i,ks,j),i=iss,iee-1)
                if (order(1:1).eq.'N') read(9,fmt1) (az(i,ks,j),i=iee-1,iss,-1)
            endif  
        enddo
    elseif (order(2:2).eq.'N') then
        do j=jee-1,jss,-1
            if (fmt1(1:1).eq.'*') then
                if (order(1:1).eq.'P') read(9,*) (az(i,ks,j),i=iss,iee-1)
                if (order(1:1).eq.'N') read(9,*) (az(i,ks,j),i=iee-1,iss,-1)
            else
                if (order(1:1).eq.'P') read(9,fmt1) (az(i,ks,j),i=iss,iee-1)
                if (order(1:1).eq.'N') read(9,fmt1) (az(i,ks,j),i=iee-1,iss,-1)
            endif         
        enddo
    endif
CLOSE(9) 
!
if (iee.ne.ie.or.jee.ne.je) then
!Expand mannings data to fit the mesh size
    allocate(aztemp(iss:iee-1,jss:jee-1))
    aztemp(:,:)=az(iss:iee-1,ks,jss:jee-1)
    do i=is,ie-1
        do j=js,je-1
            xc=0.5d0*(xi(i)+xi(i+1));yc=0.5d0*(y(j)+y(j+1))  
            do n=iss,iee-1
                xb=dxx*(n-iss+1)
                if (xc.lt.xb.and.xc.ge.xb-dxx) then
                    do m=jss,jee-1
                        yb=dyy*(m-jss+1)
                        if (yc.lt.yb.and.yc.ge.yb-dyy) goto 4
                    enddo
                endif
            enddo
4           az(i,ks,j)=aztemp(n,m)               
        enddo
    enddo
endif
return
200 continue
az(is:ie-1,ks,js:je-1)=0.025d0
write(6,*) 'no fm_data, set to n = 0.025 by default'
endsubroutine read_mannings
!
subroutine teibou_hand
use variables,only: is,ie,js,je, ZERO
use arrays,only:    x,ax,ay,az,xi,y
implicit none
integer :: mmn,i,j
real*8 :: slope,slopeu,slopel
!real*8,allocatable::y(:),xi(:)
!Allocate teibou_h matrix if not already
if (.not.allocated(teibou_h)) then
allocate(teibou_h(is:ie,js:je));teibou_h = ZERO
endif
!Loop over all cells
do i=is,ie
    do j=js,je
        !Loop over number of breakwaters
        do mmn=1,nbw
            !
        if (xi(i).ge.bwval(mmn)%xstart.and.xi(i).le.bwval(mmn)%xend) then 
            if (bwval(mmn)%xstart.eq.bwval(mmn)%xend) then
                !When we have infinite slope
                if (y(j).ge.bwval(mmn)%ystart.and.y(j).le.bwval(mmn)%yend) then
                    teibou_h(i,j)=bwval(mmn)%height
                endif
            else
                slope = (bwval(mmn)%yend-bwval(mmn)%ystart)/&
                        (bwval(mmn)%xend-bwval(mmn)%xstart)
                if (slope.eq.ZERO) then
                    if (y(j).eq.bwval(mmn)%ystart) then
                        teibou_h(i,j) = bwval(mmn)%height 
                    endif
                elseif (slope.gt.0.0d0) then
                    !When the slope is positive
                    if (y(j).ge.bwval(mmn)%ystart.and.y(j).le.bwval(mmn)%yend) then
                    !Within breakwater region
                        if (xi(i).eq.bwval(mmn)%xstart.and.y(j).eq.bwval(mmn)%ystart) goto 113
                        if (xi(i).eq.bwval(mmn)%xstart) cycle
                        slopeu=(y(j+1)-bwval(mmn)%ystart)/(xi(i+1)-bwval(mmn)%xstart)
                        slopel=(y(j)-bwval(mmn)%ystart)/(xi(i)-bwval(mmn)%xstart)
                        if (slopeu.gt.slope.and.slopel.lt.slope) then
                        !On breakwater line
                113     teibou_h(i,j) = bwval(mmn)%height
                        endif
                    endif
                else
                !When the slope is negative
                    if (y(j).ge.bwval(mmn)%yend.and.y(j).le.bwval(mmn)%ystart) then
                    !Within breakwater region
                        if (xi(i).eq.bwval(mmn)%xstart.and.y(j+1).eq.bwval(mmn)%ystart) goto 114
                        if (xi(i).eq.bwval(mmn)%xstart) cycle
                        slopeu=(y(j+1)-bwval(mmn)%ystart)/(xi(i+1)-bwval(mmn)%xstart)
                        slopel=(y(j)-bwval(mmn)%ystart)/(xi(i)-bwval(mmn)%xstart)
                        if (slopeu.gt.slope.and.slopel.lt.slope) then
                        !On breakwater line
                114     teibou_h(i,j) = bwval(mmn)%height 
                        endif
                    endif
                endif
            endif
        endif
        enddo
    enddo
enddo
endsubroutine teibou_hand
!
subroutine alter_w_teibou
use variables, only: is,ie,js,je,ks,ke,ZERO,VSMALL,TENTH
use arrays, only:    ax,ay,az,z,nff,mn
implicit none
integer :: i,j,n,k,ix,iy
real*8  :: htemp,zktemp
!if integer lets allocate teibou_h and put in there..
if (.not.allocated(teibou_h)) then
    allocate(teibou_h(is:ie-1,js:je-1))
    do i=is,ie-1
        do j=js,je-1
            teibou_h(i,j) = ZERO
            imod = teibou_i(i,j) / 1000
            !convert integer to real teibou height if exists
            if (imod.ne.0) teibou_h(i,j) = real ( teibou_i(i,j) - (imod * 1000) ) * TENTH
        enddo
    enddo
endif
if (nvchange.ne.0) then
    where(teibou_h.gt.oldval-VSMALL.and.teibou_h.lt.oldval+VSMALL) teibou_h = nval
endif
!
ix=0;iy=0
do i=is,ie
    if (allocated(teibou_i).and.i.eq.ie) cycle
    do j=js,je
        if (allocated(teibou_i).and.j.eq.je) cycle
        if (nchange.ne.0) then            
            do n=1,nchange
                if (i.eq.valchange(n)%ii.and.j.eq.valchange(n)%jj) then
                    teibou_h(i,j) = valchange(n)%newval
                    ix = valchange(n)%ix ; iy = valchange(n)%iy; exit
                endif
            enddo
        endif        
        if (teibou_h(i,j).ne.ZERO) then
            do k=ks,ke-1
                if (nff(mn(i,k,j))%f.eq.-1.and.nff(mn(i,k,j))%b.eq.0) cycle
                if (z(k).gt.teibou_h(i,j)) then
                    ! Nothing
                elseif (z(k).lt.teibou_h(i,j)) then
                    if (.not.allocated(teibou_i)) then
                        !We do not know whether the teibou is on south side or west side
                        !and we guess from the surrounding teibou values
45                      if (ix+iy.gt.0) then
                            if (ix.eq.1)  ax(i,k,j) = (z(k+1)-teibou_h(i,j))/(z(k+1)-z(k))
                            if (iy.eq.1)  ay(i,k,j) = (z(k+1)-teibou_h(i,j))/(z(k+1)-z(k))   
                            ix = 0 ; iy = 0
                        else
                            if (j.ne.js.and.j.ne.je) then
                                if (i.eq.is) then
                                    ax(i,k,j) = (z(k+1)-teibou_h(i,j))/(z(k+1)-z(k))
                                elseif (teibou_h(i-1,j).eq.ZERO) then !West side
                                    ax(i,k,j) = (z(k+1)-teibou_h(i,j))/(z(k+1)-z(k))
                                endif
                            endif
                            if (i.ne.is.and.i.ne.ie) then
                                if (j.eq.js) then
                                    ay(i,k,j) = (z(k+1)-teibou_h(i,j))/(z(k+1)-z(k))
                                elseif (teibou_h(i,j-1).eq.ZERO) then !South side
                                    ay(i,k,j) = (z(k+1)-teibou_h(i,j))/(z(k+1)-z(k))
                                endif
                            endif
                        endif
                        if (ay(i,k,j).lt.ZERO) ay(i,k,j) = ZERO
                        if (ax(i,k,j).lt.ZERO) ax(i,k,j) = ZERO
                    else
                        !It is specified whether teibou is on north side or east side
                        !and we can enter in the value as so...
                        imod = teibou_i(i,j) / 1000 ; if (imod.eq.0) goto 45 !Maybe we added by hand                      
                        if (imod.eq.1.or.imod.eq.3) then !East side
                            ax(i+1,k,j) = max((z(k+1)-teibou_h(i,j))/(z(k+1)-z(k)),ZERO)
                        endif                           
                        if (imod.eq.2.or.imod.eq.3) then  !West side    
                            ay(i,k,j+1) = max((z(k+1)-teibou_h(i,j))/(z(k+1)-z(k)),ZERO)
                        endif
                    endif  
                endif
            enddo
        endif
    enddo
enddo
endsubroutine alter_w_teibou
end program add_teibou
