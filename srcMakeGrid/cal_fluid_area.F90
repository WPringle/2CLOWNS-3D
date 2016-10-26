subroutine cal_fluid_area(zmax,poly)  
use variables,only:ndri,is,js,ks,ie,je,ke,fmt1,refval,order,title,input
use arrays,only:ivt,upl,vt,xi,y,z,ipl,ipo,npo,OBJECTshape,nf3,nfb3,fb3
USE local_module
!use ojt
implicit none
integer::i,j,k,ioerr
integer,allocatable::info(:,:)
doubleprecision,allocatable,dimension(:,:)::zk,zh0,zh1
integer::i2,j2,ii,kdc,n1,n2,nx1,nx2
logical::poly,exist !,zhdata = .false.
doubleprecision::zh(ie,je),zmax,cc(0:1),c1,c2,c3,c4
doubleprecision::zhmin,dxx
doubleprecision,dimension(0:2)::a1,a2,a3,a4,tmpV,ccc,cc1,cc2,cc3,bb1,bb2,bb3
doubleprecision,dimension(0:1)::b1,b2,b3,alp,b4,bc1,bc2,bc3,bc4
integer :: is1,js1,ie1,je1,nf_poly(ie,je)
real*8 :: zg,wdiv
!
!When we define the topography bottom from data file or slope in 2D
n1=len_trim(OBJECTshape(1))
if (OBJECTshape(1)(n1-2:n1) == 'dat'.or.&
    OBJECTshape(1)(n1-2:n1) == 'txt'.or.&
    OBJECTshape(1)(n1-2:n1) == 'csv') then
    !Read from data file if exists
    call read_DAT
elseif (OBJECTshape(1)(n1-2:n1) == 'xyz') then
    !Read from xyz file
    call read_XYZ_file
elseif (OBJECTshape(1)(n1-4:n1) == 'odata') then
    if (poly) zh = -10000d0 !For polygon
else
    write(6,*) 'Warning: No input file -> will proceed with zh = z(ks) everywhere'
    zh = z(ks) !Set zh equal to bottom in first place in case of no data
endif
k=3
do i=is,ie-1
    do j=js,je-1
        if (slopenum.gt.0) then
            !Define the ground level for a plane or composite beach
            do i2=1,slopenum
                if (xi(i).lt.xstart(i2)) cycle
                zg=z(k)
                if (i2.gt.1) then
                    do j2=2,i2
                        zg=zg+wdiv(xstart(j2)-xstart(j2-1),beta(j2-1))
                    enddo
                endif
                zh(i,j)=zg+wdiv(0.5d0*(xi(i)+xi(i+1))-xstart(i2),beta(i2))
            enddo
        endif
    enddo
enddo

POLY_IF: if (poly) then
    allocate( info(ndri,maxval(ivt(1:ndri))) );info=0
    allocate( upl(ndri,maxval(ipl(1:ndri)),1:4) )
    nf_poly = 0
POLY_DO: do kdc=1,ndri
    !delete unnecessary parts of polygon 
    dxx = 0.0d0
    do i=1,ivt(kdc)-1
        if (vt(kdc,i,1).ge.xi(is).and.vt(kdc,i,1).le.xi(ie)&
            .and.vt(kdc,i,2).ge.y(js).and.vt(kdc,i,2).le.y(je)) then
            dxx=max(dxx,sqrt((vt(kdc,i,1)-vt(kdc,i+1,1))**2+&
                             (vt(kdc,i,2)-vt(kdc,i+1,2))**2))
        endif
    enddo
    if (is.ne.ie-1) then
        call deletepoly(kdc,xi(is)-dxx,1,'<')
        call deletepoly(kdc,xi(ie)+dxx,1,'>')
    endif
    if (js.ne.je-1) then
        call deletepoly(kdc,y(js)-dxx,2,'<')
        call deletepoly(kdc,y(je)+dxx,2,'>')  
    endif
!   
    write(6,*) 'Go in unv'
    call unv(kdc)
    write(6,*) 'Out unv'
    do i2=1,ipl(kdc)
        
        if(upl(kdc,i2,3).le.0.0d0) then
            continue
            ipo(kdc,i2)=0
        endif    
        if(ipo(kdc,i2).eq.0) cycle
        do j2=1,ipo(kdc,i2)
            if (info(kdc,npo(kdc,i2,j2)).eq.0) then
                info(kdc,npo(kdc,i2,j2))=1
            endif
        enddo
    enddo
enddo POLY_DO

!$omp parallel do private(j,cc,kdc,zg) !schedule(dynamic)
do i=is,ie
    do j=js,je
        if (ks == ke-1.or.flat.eq.1) then
        if (i.eq.ie.or.j.eq.je) cycle
            cc(:)= [0.5d0*(xi(i)+xi(i+1)),0.5d0*(y(j)+y(j+1))]
        else
            cc(:)= [xi(i),y(j)]
        endif
        zg = -1d4
        do kdc=1,ndri
            if (ivt(kdc).eq.0) cycle
            !Call hyoko2 when we have no value or it is equal
            !to the reference value
            if (zg.eq.-1d4.or.zg.eq.zmax) then
                call hyoko2(kdc,cc,zg,zmax)   
            endif
            if (abs(zg).ne.1d4.and.zg.gt.zh(i,j).and.zg.ne.zmax) then
                nf_poly(i,j) = 1 !zh(i,j) = zg ; 
            endif
        enddo
    enddo
enddo
!$omp end parallel do
endif POLY_IF
!
if(.not.allocated(nf3)) then
    allocate(nf3(is-2:ie+2,1:ke+2,js-2:je+2),nfb3(is-2:ie+2,1:ke+2,js-2:je+2));nf3=-1;nfb3=0
    allocate(fb3(is-2:ie+2,1:ke+2,js-2:je+2));fb3=0.d0
else
    allocate(zk(is-1:ie,js-1:je))
    do i=is,ie-1
    do j=js,je-1
        zk(i,j)=z(ke-1)
        do k=ke-1,ks,-1
            if(nf3(i,k,j).ge.0.and.nf3(i,k-1,j).eq.-1) then
                zk(i,j)=z(k)+(1.0d0-fb3(i,k,j))*(z(k+1)-z(k))
                exit
            endif
        enddo
        do k=ks,ke-1
            if(nf3(i,k,j).ge.0) then
                nf3(i,k,j)=-1;nfb3(i,k,j)=0
            endif
        enddo
    enddo
    enddo
    continue
endif
write(6,*) 'hyodo2 end'
#ifdef DDDDD
if(.not.zhdata) then
call owari(title,nx2);call hajime(title,nx1)

open(9,file=title(nx1:nx2)//'.zhdata',form='unformatted')
    write(9) is,js,ie,je
    write(9) zh(is:ie,js:je)
close(9)

write(6,*) 'output zhdata to',title(nx1:nx2)//'.zhdata'
continue
endif
#endif

!$omp parallel do private(j,kdc,a1,a2,a3,a4,bc1,bc2,bc3,bc4,cc,b1,b2,b3,b4,c1,c2,c3,c4,zhmin,ii,k) schedule(dynamic)
do i=is,ie-1
do j=js,je-1
    if (ks == ke-1.or.flat.eq.1) then
        zhmin = zh(i,j)
    else
        zhmin = min(zh(i,j),zh(i+1,j),zh(i+1,j+1),zh(i,j+1))
        if (poly) then
            a1(:)=(/xi(i+1),y(j),0.0d0/)-(/xi(i),y(j),0.0d0/)
            a2(:)=(/xi(i+1),y(j+1),0.0d0/)-(/xi(i+1),y(j),0.0d0/)
            a3(:)=(/xi(i),y(j+1),0.0d0/)-(/xi(i+1),y(j+1),0.0d0/)
            a4(:)=(/xi(i),y(j),0.0d0/)-(/xi(i),y(j+1),0.0d0/)
            bc1(:)=(/xi(i),y(j)/)
            bc2(:)=(/xi(i+1),y(j)/)
            bc3(:)=(/xi(i+1),y(j+1)/)
            bc4(:)=(/xi(i),y(j+1)/)
            !
            do kdc=1,ndri
                if (ndri.gt.1.and.kdc.eq.ndri) then
                    if (zhmin.gt.refval) exit
                endif
                do ii=1,ivt(kdc)
                if(info(kdc,ii).eq.0)  cycle
                cc(:)=vt(kdc,ii,1:2)
                if(cc(0).lt.xi(i)) cycle
                if(cc(1).lt.y(j)) cycle
                if(cc(0).gt.xi(i+1)) cycle
                if(cc(1).gt.y(j+1)) exit
                b1(:)=cc(:)-bc1(:)
                b2(:)=cc(:)-bc2(:)
                b3(:)=cc(:)-bc3(:)
                b4(:)=cc(:)-bc4(:)
                c1=a1(0)*b1(1)-a1(1)*b1(0)
                c2=a2(0)*b2(1)-a2(1)*b2(0)
                c3=a3(0)*b3(1)-a3(1)*b3(0)
                c4=a4(0)*b4(1)-a4(1)*b4(0)
                if(c1.ge.0.0d0.and.c2.ge.0.0d0.and.c3.ge.0.0d0.and.c4.ge.0.0d0) then
                zhmin=min(vt(kdc,ii,3),zhmin)
                endif
                enddo
            enddo
        endif
    endif
#ifdef DDDDD
    if(.not.zhdata) then    
        if(allocated(zk)) then
            if(zhmin.ne.10000.d0) then
            zhmin=max(min(zmax,zhmin),zk(i,j))    
            else
            zhmin=zk(i,j)  
            endif
        endif
    endif
#endif    
    do k=ks,ke-1
        if(z(k+1)-zhmin.gt.1.0d-10&
#ifdef CYLINDER
           .or.(zhmin.ne.-1d4.and.z(k)+zhmin.lt.-1.0d-10)&
#endif      
            ) then
    
            if (nf3(i,k,j).eq.-1) then
                  fb3(i,k,j)=min(1.0d0,max(0.0d0,(z(k+1)-zhmin)/(z(k+1)-z(k))))
#ifdef CYLINDER
                 if (zhmin.ne.-1d4.and.z(k)+zhmin.lt.-1.0d-10) then
                    fb3(i,k,j) = 1d0 - min(1.0d0,max(0.0d0,1d0+(z(k)+zhmin)/(z(k+1)-z(k))))
                 endif
#endif
                if (ks == ke-1.or.flat.eq.1) then
                      if(ks == ke-1.and.(z(ke)-z(ks))*fb3(i,k,j).lt.0.01d0) then
                          fb3(i,k,j)=0.0d0                   
                      else if(ks.ne.ke-1.and.fb3(i,k,j).lt.0.01d0) then
                          fb3(i,k,j)=0.0d0    
                          if (poly) then
                            if (nf_poly(i,j) == 1) then
                                nfb3(i,k,j) = 10
                            endif
                          endif
                      else
                          nf3(i,k,j)=0;nfb3(i,k,j)=0
                      endif
                else
                    nf3(i,k,j)=0;nfb3(i,k,j)=0
                endif
            endif
         else 
            if (poly) then
            if (nf_poly(i,j) == 1) then
                nfb3(i,k,j) = 10
            endif
            endif
     !       write(6,*) i,k,j
        endif
    enddo
enddo
enddo
!$omp end parallel do

contains
    
    subroutine read_DAT
        logical :: stopflag=.false.
         write(6,*) 'start read_DAT'
        open(10,file=OBJECTshape(1),status='old')
            if (order(2:2).eq.'P') then
                do j=js,je-1
                    if (order(1:1).eq.'P') then
                        if (fmt1(1:1).eq.'*') then
                            read(10,*) (zh(i,j),i=is,ie-1)
                        else
                            read(10,fmt1) (zh(i,j),i=is,ie-1)
                        endif   
                    elseif (order(1:1).eq.'N') then
                        if (fmt1(1:1).eq.'*') then
                            read(10,*) (zh(i,j),i=ie-1,is,-1)
                        else
                            read(10,fmt1) (zh(i,j),i=ie-1,is,-1)
                        endif  
                    else
                        write(6,*) 'Must define N or P at end of file name'
                        stop                    
                    endif
                enddo
            elseif (order(2:2).eq.'N') then
                do j=je-1,js,-1
                    if (order(1:1).eq.'P') then
                        if (fmt1(1:1).eq.'*') then
                            read(10,*) (zh(i,j),i=is,ie-1)
                        else
                            read(10,fmt1) (zh(i,j),i=is,ie-1)
                        endif   
                    elseif (order(1:1).eq.'N') then
                        if (fmt1(1:1).eq.'*') then
                            read(10,*) (zh(i,j),i=ie-1,is,-1)
                        else
                            read(10,fmt1) (zh(i,j),i=ie-1,is,-1)
                        endif  
                    else
                        write(6,*) 'Must define N or P at end of file name'
                        stop                    
                    endif
                enddo     
            else
                write(6,*) 'Must define N or P at end of file name'
                stop
            endif
        close(10)
        if (order(3:3).eq.'N') zh=-zh
        !Check minimum depths with defined z
        do i=is,ie-1
            do j=js,je-1
                if (zh(i,j) < z(ks)) then
                    write(6,*) 'Read data less than defined min z=>',zh(i,j)
                    stopflag = .true.
                endif
            enddo
        enddo
         write(6,*) 'end read_DAT'  !,stopflag
        if (stopflag) stop
    endsubroutine read_DAT

    subroutine read_XYZ_file
        real*8 :: xx, yy
        integer :: length, count, m, tot
        integer,allocatable :: loop(:,:)
        write(6,*) 'start read_XYZ_file'
        !Read loop data
        !open(11,file=OBJECTshape(1)(n1:n2-6)//'.loop.dat',status='old') 
        !    read(11,*) length
        !    allocate(loop(js:je-1,1:length))
        !    do j=js,je-1
        !        read(11,*) (loop(j,i),i=1,length)
        !    enddo
        !close(11)
        !if (OBJECTshape(1)(n2-6:n2-6).eq.'N') then
        !    zh = -z(ke)
        !else
        !    zh = z(ke)   
        !endif
        !loop(320,5) = loop(320,5) - 10
        open(10,file=OBJECTshape(1),status='old')
            if (order(2:2).eq.'P') then
                do j=js,je-1
                    !count = 0
                    !write(6,*) 'j=',j
                    if (order(1:1).eq.'P') then
                        !if (j.ge.37) then
                        !    do i=66,ie-1
                        !        read(10,*) xx,yy,zh(i,j)
                        !    enddo
                        !    do i=is,65
                        !        read(10,*) xx,yy,zh(i,j)
                        !    enddo
                        !    cycle
                        !endif
                        do i=is,ie-1
                            !if (j.ge.14.and.j.lt.37.and.i.lt.ie-18) then
                            !    cycle
                            !endif
                           ! tot = 0
                            !do m=1,length,2
                            !    tot = tot + loop(j,m)
                            !enddo
                            !if (j.lt.151.and.i.lt.ie-tot) then
                            !    !We need to skip the first few values
                            !    !Keep ground level as equal to max value
                            !    cycle
                            !elseif (j.ge.326) then
                            !    if (i.lt.ie-loop(j,1)) cycle
                            !else
                            !    do m=2,length-1,2
                            !        if (loop(j,m).eq.0) exit
                            !        if (j.ge.151) then
                            !            if (m.eq.2) then
                            !                !We need to count the ground level first
                            !                if (i.ge.is+loop(j,2).and.i.lt.ie-tot+sum(loop(j,2:3))) goto 45
                            !            else
                            !                !Skip at the specified locations
                            !                if (i.ge.ie-tot+sum(loop(j,2:m))-sum(loop(j,2:m-2:2)).and.i.lt.ie-tot+sum(loop(j,2:m+1))-sum(loop(j,2:m-2:2))) goto 45                                            
                            !            endif
                            !        else
                            !            !Skip at the specified locations
                            !            if (i.ge.ie-tot+sum(loop(j,2:m))-sum(loop(j,2:m-2:2)).and.i.lt.ie-tot+sum(loop(j,2:m+1))-sum(loop(j,2:m-2:2))) goto 45
                            !        endif
                            !    enddo
                            !endif  
                            !Read the ground level
!                            count = count + 1
                            if (fmt1(1:1).eq.'*') then
                                read(10,*) xx,yy,zh(i,j)
                            else
                                read(10,fmt1) xx,yy,zh(i,j)
                            endif
                            !read(10,*) count,xx,yy,zh(i,j)
!45                          continue                            
                        enddo
                        !if (count.ne.loop(j,1)) then
                        !    continue
                        !endif
                    elseif (order(1:1).eq.'N') then
                        do i=ie-1,is,-1
                            if (fmt1(1:1).eq.'*') then
                                read(10,*) xx,yy,zh(i,j)
                            else
                                read(10,fmt1) xx,yy,zh(i,j)
                            endif
                        enddo
                    else
                        write(6,*) 'Must define N or P at end of file name'
                        stop                    
                    endif
                enddo
            elseif (order(2:2).eq.'N') then
                do j=je-1,js,-1
                    if (order(1:1).eq.'P') then
                        do i=is,ie-1
                            if (fmt1(1:1).eq.'*') then
                                read(10,*) xx,yy,zh(i,j)
                            else
                                read(10,fmt1) xx,yy,zh(i,j)
                            endif
                        enddo
                    elseif (order(1:1).eq.'N') then
                        do i=ie-1,is,-1
                            if (fmt1(1:1).eq.'*') then
                                read(10,*) xx,yy,zh(i,j)
                            else
                                read(10,fmt1) xx,yy,zh(i,j)
                            endif
                        enddo
                    else
                        write(6,*) 'Must define N or P at end of file name'
                        stop                    
                    endif
                enddo     
            else
                write(6,*) 'Must define N or P at end of file name'
                stop
            endif
        close(10)
        if (order(3:3).eq.'N') zh=-zh  
        write(6,*) 'end read_XYZ_file'
    endsubroutine read_XYZ_file 
    
    end subroutine
    

subroutine deletepoly(kdc,refval,refint,mathsym)
use variables,only:ndri
use arrays,only:ipo,ipl,vt,ivt,npo
implicit none
real*8,intent(in)::refval
integer,intent(in)::refint,kdc
character*1,intent(in)::mathsym
real*8::cc(0:1),zh
integer::ivtp,iplp,i,j,ipmax
integer,allocatable::ipop(:,:),npop(:,:,:),ipoc(:)
doubleprecision,allocatable::vtp(:,:,:)
!Delete the region of the polygon where the value is '=' or '<' or '>' to 
!the refval, of the refint dimension of the polygon and remake polygon data

!or

!Set refval = -10000d0 and Delete the region of polygon 'kdc' where it 
!overlaps the polygon 'refint' on the x-y plane
!===========================================================================
!Check inputs
if (refval.ne.-10000d0.and.(refint.gt.3.or.refint.lt.1.or.&
   (mathsym.ne.'='.and.mathsym.ne.'>'.and.mathsym.ne.'<'))) then 
        write(6,*) 'inputs to deletepoly are not supported:',refint,mathsym
        return
endif    
ipmax=3    
do i=1,ndri    
ipmax=max(ipmax,maxval(ipo(i,1:ipl(i))))
enddo
allocate(vtp(ndri,maxval(ivt(1:ndri)),1:3),ipop(ndri,maxval(ipl(1:ndri))),&
         npop(ndri,maxval(ipl(1:ndri)),ipmax),ipoc(ivt(kdc)))
vtp=vt;ipop=ipo;npop=npo
!
ivtp=0;ipoc=0
do i=1,ivt(kdc)
    if (refval.eq.-10000d0) then
        !Remove the vertices of the overlapping region
        cc(:)=(/vt(kdc,i,1),vt(kdc,i,2)/);zh=-10000d0
        call hyoko2(refint,cc,zh,zh)  
        if (zh.ne.-10000d0) cycle
    else
        !Remove the vertices of the specified region
        if (mathsym.eq.'='.and.vt(kdc,i,refint).eq.refval) cycle
        if (mathsym.eq.'<'.and.vt(kdc,i,refint).lt.refval) cycle
        if (mathsym.eq.'>'.and.vt(kdc,i,refint).gt.refval) cycle
    endif
    ivtp = ivtp+1
    vtp(kdc,ivtp,:) = vt(kdc,i,:)
    ipoc(i) = ivtp !To convert the new vertex number with the old
enddo
!Remove planes that contains vertices that no longer exist (were deleted)
iplp=0
do i=1,ipl(kdc)
    do j=1,ipo(kdc,i)
        if (ipoc(npo(kdc,i,j)).eq.0) goto 3
    enddo
        iplp=iplp+1
    do j=1,ipo(kdc,i)
        npop(kdc,iplp,j)=ipoc(npo(kdc,i,j))
    enddo
        ipop(kdc,iplp)=ipo(kdc,i)
3 continue
enddo
deallocate(vt,ipo,npo,ipoc)
!Reallocate and write polygon arrays
ivt(kdc)=ivtp;ipl(kdc)=iplp
write(6,*) 'kdc=',kdc,'New ivt(kdc)=',ivt(kdc),'New ipl(kdc)=',ipl(kdc)
allocate(vt(ndri,maxval(ivt(1:ndri)),1:3))
allocate(ipo(ndri,maxval(ipl(1:ndri))))
allocate(npo(ndri,maxval(ipl(1:ndri)),ipmax))
vt=vtp;ipo=ipop;npo=npop
deallocate(vtp,ipop,npop)
end subroutine deletepoly
