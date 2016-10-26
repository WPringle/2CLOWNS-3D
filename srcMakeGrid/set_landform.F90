subroutine set_chikei(hini) !,adjnum,adjinfo)
USE variables,only:is,ie,js,je,ks,ke,openbound
USE arrays,only:nf3,nfb3,ax,ay,az,fb3,xi,y,z,f3
USE local_module
implicit none
integer :: i,j,k
real*8 :: hini,hinitemp
integer :: ii  !adjnum
!type adj
!real*8  :: hiniadj,xs,xe,ys,ye
!integer :: type
!endtype adj
!type(adj) :: adjinfo(adjnum)
#ifdef MESHKAT
integer::outletbound(4)
#endif
!
if(.not.allocated(ax)) then
    allocate(ax(is-2:ie+2,1:ke+2,js-2:je+2));ax=0.0d0
    allocate(ay(is-2:ie+2,1:ke+2,js-2:je+2));ay=0.0d0
    allocate(az(is-2:ie+2,1:ke+2,js-2:je+2));az=0.0d0
    continue
endif
!Open the aperture ratios where fluid exists
do i=is,ie-1
do j=js,je-1
do k=ks,ke-1
if(nf3(i,k,j).eq.-1) cycle
if (ks == ke-1) then
    if(nf3(i-1,k,j).ge.0)  ax(i,k,j)=1.0d0
    if(nf3(i,k,j-1).ge.0)  ay(i,k,j)=1.0d0
    if(nf3(i,k-1,j).ge.0)  az(i,k,j)=1.0d0
else
    if(nf3(i-1,k,j).ge.0)  ax(i,k,j)=min(fb3(i,k,j),fb3(i-1,k,j))
    if(nf3(i,k,j-1).ge.0)  ay(i,k,j)=min(fb3(i,k,j),fb3(i,k,j-1))
    if(nf3(i,k-1,j).ge.0)  az(i,k,j)=1.0d0
endif
enddo
enddo
enddo
!
!==============================================!
!Calculate the F values where fluid exists     !
!==============================================!
do i=is,ie-1 !Loop through all
do j=js,je-1 !cells in X-Y plane
    hinitemp = hini
    if (adjnum.gt.0) then
        do ii=1,adjnum
            if (xi(i).ge.adjinfo(ii)%xs.and.xi(i).le.adjinfo(ii)%xe.and.&
                y(j).ge.adjinfo(ii)%ys.and.y(j).le.adjinfo(ii)%ye) then
                hinitemp = adjinfo(ii)%hiniadj
                exit
            endif
        enddo
    endif
    do k=ks,ke-1 !Loop over all the cells in the vertical direction
        if (nf3(i,k,j).ne.0) cycle !Cycle at object cell
        if (z(k).ge.hinitemp) then
            !Set to Air cell
            f3(i,k,j) = 0.0d0
        elseif (z(k+1).le.hinitemp) then
            !Set to Fluid Cell
            f3(i,k,j) = 1.0d0
            nf3(i,k,j) = 1
        else
            if (fb3(i,k,j).gt.( z(k+1) - hinitemp )/ ( z(k+1) - z(k) )) then
                !Set to Surface Cell
                if (adjnum.gt.0) then
                    if (adjinfo(ii)%type.eq.1) then
                        !When we want to only remove water where ground level is very high
                        f3(i,k,j) = ( hini - z(k+1) + fb3(i,k,j) * ( z(k+1) - z(k) ) )&
                                 / ( ( z(k+1) - z(k) ) * fb3(i,k,j) )
                    elseif (adjinfo(ii)%type.eq.0) then
                        !When we just want to adjust the water level in a certain region uniformly
                        f3(i,k,j) = ( hinitemp - z(k+1) + fb3(i,k,j) * ( z(k+1) - z(k) ) )&
                                 / ( ( z(k+1) - z(k) ) * fb3(i,k,j) )
                    endif
                else
                    !Calculate in normal way
                    f3(i,k,j) = ( hinitemp - z(k+1) + fb3(i,k,j) * ( z(k+1) - z(k) ) )&
                                 / ( ( z(k+1) - z(k) ) * fb3(i,k,j) )
                endif
                nf3(i,k,j) = 2; nfb3(i,k,j) = -3   
            else
                !Set to Air cell
                f3(i,k,j) = 0.0d0
            endif
        endif
    enddo
enddo
enddo
!
!Set the openness of the various boundaries
!!!!!!!!!!ó¨ì¸ã´äEx!!!!!!
if (openbound%West.ne.0) then
    do j=js,je
    do k=ks,ke-1
    if(nf3(is,k,j).ge.0) then
        if(ks.eq.ke-1) ax(is,k,j)=1.0d0
        if(ks.ne.ke-1) ax(is,k,j)=fb3(is,k,j)
    nfb3(is-1,k,j)=openbound%West
    endif
    enddo
    enddo
endif
!!!!!!!!!!ó¨èoã´äEx!!!!!!
if (openbound%East.ne.0) then
    do j=js,je
    do k=ks,ke-1
    if(nf3(ie-1,k,j).ge.0) then
                if(ks.eq.ke-1)  ax(ie,k,j)=1.0d0
                if(ks.ne.ke-1)  ax(ie,k,j)=fb3(ie-1,k,j)
    nfb3(ie,k,j)=openbound%East
    endif
    enddo
    enddo
endif
!!!!!!!!!!ó¨ì¸ã´äEy!!!!!!
if (openbound%South.ne.0) then
    do i=is,ie
    do k=ks,ke-1
    if(nf3(i,k,js).ge.0) then
    if(openbound%South.gt.0) then
            if(ks.eq.ke-1) ay(i,k,js)=1.0d0
            if(ks.ne.ke-1) ay(i,k,js)=fb3(i,k,js)
    endif
    nfb3(i,k,js-1)=openbound%South
    endif
    enddo
    enddo
endif
!!!!!!!!!!ó¨èoã´äEx!!!!!!
if (openbound%North.ne.0) then
    do i=is,ie
    do k=ks,ke-1
    if(nf3(i,k,je-1).ge.0) then
    if(openbound%North.gt.0) then
            if(ks.eq.ke-1) ay(i,k,je)=1.0d0
            if(ks.ne.ke-1) ay(i,k,je)=fb3(i,k,je-1)
    endif
    nfb3(i,k,je)=openbound%North
    endif
    enddo
    enddo
endif
!!!!!!!!!!Corners!!!!!!!!!!!!!!
if (ks == ke-1) then
    if (openbound%North.ne.0.and.&
        openbound%East.ne.0) then
        nfb3(ie,ks,je)=openbound%North
    endif
    if (openbound%South.ne.0.and.&
        openbound%East.ne.0) then
        nfb3(ie,ks,js-1)=openbound%South
    endif
    if (openbound%North.ne.0.and.&
        openbound%West.ne.0) then
        nfb3(is-1,ks,je)=openbound%North
    endif
    if (openbound%South.ne.0.and.&
        openbound%West.ne.0) then
        nfb3(is-1,ks,js-1)=openbound%South
    endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Bottom Outlet!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MESHKAT
open(5,file='../dat-win/bottomoutlet.txt',status='old')
    read(5,*) (outletbound(i),i=1,4)
    do k=outletbound(4),outletbound(3),-1
    read(5,*) (ax(ie,k,j),j=outletbound(1),outletbound(2))
    enddo
close(5)
nfb3(ie,outletbound(3):outletbound(4),outletbound(1):outletbound(2))=3
#endif
!
end subroutine