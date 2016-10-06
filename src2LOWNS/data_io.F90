!%%%%%%%%%%%%%%%%%%%% FILE: data_io.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads input data and control files and allocates arrays
subroutine read_control_file
    use variables
    implicit none
#ifdef BUGWIN
    open(5,file='debug.txt',status='old')
#endif
      read(5,*) IM; print *,'IM=',IM
      if (IM.le.0.or.IM.gt.3) stop 'Error: IM model type is not valid'
      read(5,*) nkai;print *,'nkai=',nkai
      read(5,*) mwr,dwr;print *,'mwr=',mwr,'dwr',dwr
      read(5,*) tw0;print *,'tw0=',tw0
      read(5,*) dtw;print *,'dtw=',dtw
      read(5,*) ko_tmp
      read(5,*) ko_out; print *,'output sfdata (yes =1)', ko_out
      read(5,*) t_end;print *,'t_end=',t_end
      read(5,*) g_read;print *,'g_read=',g_read
      if (IM.eq.1.or.IM.eq.3) then
          read(5,*) nu_read;print *,'nu_read=',nu_read
          read(5,*) dt,dtmax,cl_set;print *,'dt=',dt,'dtmax=',dtmax
          read(5,*) d0;print *,'d0=',d0
          read(5,*) eps0;print *,'eps0=',eps0
          read(5,'(A255)') input;print *,'input=',trim(input)
          read(5,'(A255)') inputxyz;print *,'input_xyz=',trim(inputxyz)
          read(5,'(A255)') inputf0;print *,'inputf0=',trim(inputf0)    
#ifdef DAKU
          read(5,'(A255)') inputtc;print *,'input_tc=',trim(inputtc)
#endif
      endif
      read(5,'(A255)') testa; print *,'testa=',trim(testa)  
end subroutine read_control_file
    
subroutine read_rdata
    use variables, only: akrate, erate ,ike, inne2, t0, akmin, emin, inputke, &
                         fg_xyz, inns, inne, input, ks_r, r_lim, pro_lim,     &
                         vt_op, ZERO, TEN
    use arrays, only: r,prop,in,kn,jn,alloc_aked,d,nff
    implicit none
    integer :: nn, na1, na2
    read(5,*) akrate; print *,'akrate=',akrate
    read(5,*) erate;  print *,'erate=',erate
    read(5,*) ike;    print *,'ike=',ike
    read(5,*) ks_r(0),ks_r(10); print *,'roughness height=',ks_r(0),ks_r(10)
#ifndef NY
    read(5,*) pro_lim;print *,'production limiter type=',pro_lim
    read(5,*) vt_op;print *,'turbulent viscosity model=',vt_op
#endif
    call alloc_aked(inne2)
    if (t0.ne.ZERO.and.ike.eq.1) then
        r(:)%k  = akmin
        r(:)%e  = emin
        prop(:) = ZERO
        call owari(input,na2); call hajime(input,na1)
        INPUTke=input(na1:na2)//'.rdata'
        call owari(inputke,na2); call hajime(inputke,na1)
        write(6,*) '[HOST]read KEdata',INPUTke(na1:na2)
        if (fg_xyz.eq.1) then
            OPEN(8,FILE=INPUTke(na1:na2),status='old',form='unformatted',err=201)
            do nn = 1,inne
            if(nff(nn)%f.ne.0) read(8) r(nn)%k,r(nn)%e,d(nn),prop(nn)
            enddo
        else
            OPEN(8,FILE=INPUTke(na1:na2),status='old',form='unformatted',err=201)
            do nn=1,inne
                read(8) r(nn)%k,r(nn)%e
            enddo
            do nn=1,inne
                read(8) d(nn)
            enddo
#ifndef NOPRO
            read(8,end=200) prop(inns:inne)
#endif
        endif
        CLOSE(8)
        goto 202
    201   write(6,*) 'no rdata',INPUTke(na1:na2)
        stop
    200   continue
        prop = ZERO
    202   write(6,*) '[HOST]end KEdata',INPUTke(na1:na2)
    endif
    r_lim%k = TEN
    r_lim%e = TEN
end subroutine read_rdata
!    
subroutine read_data
use variables,only:inputxyz,input,t0,nu,nu_read,is,js,ks,ie,je,ke,fg_xyz,inns,inne&
    ,g,g_read,inne2,dt0,dt,inn2d
use arrays,only:ls,le,x,u,mn,in,jn,kn,nff,p,f,ui,v,w,p3,f3,mn2d&
    ,alloc_dxdydz,alloc_uvwfp,xi,y,z,one_dim_array,xyzn,alloc_uivwp3,in2d,nf3,nfb3
implicit none
!integer,allocatable,dimension(:,:,:)::nf,nfb
integer::na1,na2,i,k,j,ii,icc,nn,icc2=0,ist
  !    integer::i,j,k,icc=0,ii,iflag,nn
      integer::inn=0,inn2=0,inn3=0

      call owari(inputxyz,na2);call hajime(inputxyz,na1)
      write(6,*) 'read xyz-data',na1,inputxyz(na1:na2)
      if(inputxyz(na2:na2).eq."z") then
        fg_xyz=0; xyzn=.false.
      else
        fg_xyz=1
      endif

      OPEN(8,FILE=inputxyz(na1:na2),status='old',form='unformatted')
      read(8) is,js,ks,IE,je,KE
      write(6,*) is,js,ks,IE,JE,KE
      ls(0:2)=(/is,js,ks/)
      le(0:2)=(/ie,je,ke/)
      call alloc_dxdydz(is,js,ks,ie+2,je+2,ke+2)
      allocate(nf3(is-2:ie+2,ks-2:ke+2,js-2:je+2),nfb3(is-2:ie+2,ks-2:ke+2,js-2:je+2)); nf3=0;nfb3=0 ! nf,nfbの初期化
      if(one_dim_array) then
          do k=0,2
              read(8) (X(k,I),I=ls(k)-2,le(k)+1)
          enddo
      else
          read(8) (Xi(I),I=ls(0)-2,le(0)+1)
          read(8) (Y(I),I=ls(1)-2,le(1)+1)
          read(8) (Z(I),I=ls(2)-2,le(2)+1)
      endif
      if(fg_xyz.eq.0) goto 128
      read(8) icc

      do ii=1,icc
        read(8) i,k,j,nf3(i,k,j),nfb3(i,k,j)  
      enddo

128  CLOSE(8)
      
      call owari(input,na2);call hajime(input,na1)
      OPEN(8,FILE=INPUT(na1:na2),status='old',form='unformatted')
      read(8) T0,g,nu;print *,T0,g,nu
      if(t0.eq.0.0d0) then; g=g_read;nu=nu_read; endif
      if(fg_xyz.eq.0) then
        read(8) (((NF3(I,k,j),K=ks-1,KE),J=js-1,je),I=is-1,IE);print *, 'pass nf'
        read(8) (((NFB3(I,k,j),K=ks-1,KE),J=js-1,je),I=is-1,IE);print *, 'pass nfb'
        if(t0.eq.0.0d0) then
          OPEN(9,FILE=INPUT(na1:na2)//'.xyzn',status='unknown',form='unformatted')
            write(9) is,js,ks,IE,je,KE
            write(6,*) is,js,ks,IE,JE,KE
            write(9) (X(0,I),I=is-2,IE+1)
            write(9) (X(1,J),J=js-2,je+1)
            write(9) (X(2,K),K=ks-2,KE+1)
            icc=0
            do i=is-1,ie
              do j=js-1,je
                do k=ks-1,ke
                  if(nf3(i,k,j).eq.-1) icc=icc+1
                enddo
              enddo
            enddo 
            write(9) icc    
            write(6,*) 'num of nf=',icc
            do i=is-1,ie
              do j=js-1,je
                do k=ks-1,ke
                  if(nf3(i,k,j).eq.-1) then
                    write(9) i,k,j,nf3(i,k,j),nfb3(i,k,j)
                  endif
                enddo
              enddo
            enddo
          close(9)
        endif
      endif

139   icc=0;inn=0
      do j=js-1,je;do k=ks-1,ke;do i=is-1,ie
          
        if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.200) then
        if(nf3(i-1,k,j).ge.0.or.nf3(i,k-1,j).ge.0.or.nf3(i,k,j-1).ge.0 ) then
            !nfb(i,k,j)=200
            icc=icc+1
        endif
        endif           

        if(nf3(i,k,j).ge.0.or.&
                 (nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1&
                  .and.nfb3(i,k,j).lt.200)) then
                    
                if(nf3(i,k,j).eq.2) inn2=inn2+1
                if(nf3(i,k,j).eq.-1) inn3=inn3+1
                inn=inn+1;mn(i,k,j)=inn;in(inn)=i;jn(inn)=j;kn(inn)=k
        else
                mn(i,k,j)=0
        end if
      enddo;enddo;enddo
      inns=1;inne=inn

      do j=js-1,je;do k=ks-1,ke;do i=is-1,ie
        if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.200) then
          inn=inn+1 ;mn(i,k,j)=inn
          if(nfb3(i,k,j).eq.300) then
            nfb3(i,k,j)=-1
          else
            nfb3(i,k,j)=0
          endif
          in(inn)=i;jn(inn)=j;kn(inn)=k
        end if
      enddo;enddo;enddo
      inne2=inn
      
      ! for 2D counting
      allocate(in2d(0:1,(ie-is)*(je-js)),mn2d(0:ie,0:je))
      in2d = 0; inn2d = 0; mn2d = 0
      do j=js,je-1;do i=is,ie-1
        inn2d = inn2d + 1; in2d(0:1,inn2d) = [i,j]
        mn2d(i,j) = inn2d
      enddo;enddo
      
      write(6,*) 'inne=',inne
      write(6,*) 'inne2=',inne2,inne2-inne
      write(6,*) 'inn2=',inn2,inn3
      write(6,*) 'count of nfb=200',icc
        if(one_dim_array) then
              call alloc_uvwfp(inne2)
        else
              call alloc_uivwp3(ie,je,ke)
        endif
      if(t0.eq.0.0d0.and.fg_xyz.eq.-1) then  ! nagashima added after ".and."  2014.08.08
        u=0.0d0;p=0.0d0;f=0.0d0
      else
        if(fg_xyz.eq.0) then
          do inn=1,inne2
            if(inn.le.inne) then
              read(8) U(0:2,inn),F(inn),P(inn)
            else
#ifndef D2DH
              read(8) U(0:2,inn)
#endif
            endif
          enddo
        else
          read(8) icc,icc2
        if(one_dim_array) then
          do inn=1,icc
            read(8) nn,U(0:2,nn),F(nn),P(nn),nf3(in(nn),kn(nn),jn(nn)),nfb3(in(nn),kn(nn),jn(nn))
          enddo
          do inn=1,icc2
            read(8) nn,U(0:2,nn)
          enddo
        else
          do inn=1,icc
            read(8) nn,Ui(in(nn),kn(nn),jn(nn)),V(in(nn),kn(nn),jn(nn)),W(in(nn),kn(nn),jn(nn)),F3(in(nn),kn(nn),jn(nn))&
            ,P3(in(nn),kn(nn),jn(nn)),nf3(in(nn),kn(nn),jn(nn)),nfb3(in(nn),kn(nn),jn(nn))
          enddo
          do inn=1,icc2
            read(8) nn,Ui(in(nn),kn(nn),jn(nn)),V(in(nn),kn(nn),jn(nn)),W(in(nn),kn(nn),jn(nn))
          enddo
        endif
        endif
      endif

      do j=js-1,je
          do k=ks-1,ke;
              do i=is-1,ie
                  nn=mn(i,k,j)
                  if(nn.eq.0) then 
                  if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).eq.0)    cycle
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

      if(t0.gt.0.0d0) then
        read(8,iostat=ist) dt0
        if(ist.ne.0) goto 321
        write(6,*) 'read dt0=',dt0
        if(dt.eq.0.0d0) dt=dt0
      else if(dt.eq.0.0d0) then ! t0=0.0d0 && dt=0.0d0 のとき、nf3,nfb3を保存しておく。
        goto 321
      endif
      deallocate(nf3,nfb3)
321   CLOSE(8)
      
      
   write(6,*) '[HOST]read end VPFdata',INPUT(na1:na2)
end subroutine read_data 
!!!!!!!!
subroutine read_fdata
use variables,only:inne,inne2,INPUTf0,inputf,input
use arrays,only:a0,fb0,alloc_a0fb0,a,fb,alloc_afb,xyzn,drift
implicit none
      integer::nn,na1,na2,icc,icc1,ii
      if(INPUTf0.ne.'NONE') then
          call alloc_a0fb0(inne2,inne2)
          call alloc_afb(inne2,inne2)
          call owari(inputf0,na2);call hajime(inputf0,na1)
          write(6,*) '[HOST]read f0data',INPUTf0(na1:na2)
          OPEN(8,FILE=INPUTf0(na1:na2),form='unformatted', STATUS= 'old')
          if(xyzn) then
              a0(:,1:inne)=1.0d0;fb0(1:inne)=1.0d0
              a0(:,inne+1:inne2)=0.0d0  
              read(8) icc,icc1
              write(6,*) icc,icc1
              do ii=1,icc
                  read(8) nn,a0(0:2,nn),fb0(nn)
              enddo
              do ii=1,icc1
                  read(8) nn,a0(0:2,nn)  ! 必要に応じてコメントアウト
              enddo
          else
              do nn=1,inne2
                if(nn.le.inne) then
                  read(8) a0(0:2,nn),fb0(nn)
                else
                  read(8) a0(0:2,nn)  ! 必要に応じてコメントアウト
                endif
              enddo
          endif
          close(8)  
      endif
       call owari(input,na2);call hajime(input,na1)
      INPUTf=INPUT(na1:na2)//'.fdata'
      if (drift) then
          call owari(inputf,na2);call hajime(inputf,na1)
          write(6,*) '[HOST]read fdata',INPUTf(na1:na2)
          OPEN(8,FILE=INPUTf(na1:na2),form='unformatted', STATUS= 'old')
          
          if(xyzn) then
              a(:,1:inne)=1.0d0;fb(1:inne)=1.0d0
              a(:,inne+1:inne2)=0.0d0              
              read(8) icc,icc1
              write(6,*) icc,icc1
              do ii=1,icc
                  read(8) nn,a(0:2,nn),fb(nn)
              enddo
              do ii=1,icc1
                  read(8) nn,a(0:2,nn)  ! 必要に応じてコメントアウト
              enddo
          else
              do nn=1,inne2
                if(nn.le.inne) then
                  read(8) a(0:2,nn),fb(nn)
                else
                  read(8) a(0:2,nn)  ! 必要に応じてコメントアウト
                endif
              enddo
          endif
          close(8)  
      else
          a(:,1:inne2)=a0(:,1:inne2)
          fb(1:inne)=fb0(1:inne)
!          deallocate(a0,fb0)
      endif
      if(INPUTf0.eq.'NONE'.and.INPUTf.eq.'NONE') then
            write(6,*) '[HOST]Dont read fdata',INPUTf0(na1:na2)
            a0(0:2,1:inne)=1.0d0;fb0(1:inne)=1.0d0;return
      endif
end subroutine read_fdata
!!!!!!!!!
subroutine read_nsdata
      use arrays,only:sp,nor,alloc_norsp,nsdata
      use variables,only:input,sdum,mms,inne,inne2,t0,t,fg_xyz
      implicit none
      integer::nn,i,na1,na2 !,icc
        mms=-1
       call alloc_norsp(inne2+1)  
!      if(t0.ne.0.0d0) then
        call owari(input,na2);call hajime(input,na1)
!        if(t.ge.999.9d0.and.input(na2:na2).ne.".") then
!          input(na2+1:na2+1)='.'; na2=na2+1
!        endif
        OPEN(8,FILE=input(na1:na2)//'.nsdata',form='unformatted', STATUS= 'old',err=3001)
       
        if(fg_xyz.eq.1) then
          sp(1:inne,0)=sdum
          read(8) mms
          do i=1,mms  ! mms:水面セルの総数mms
            read(8) nn,nor(nn,0:3),sp(nn,0:2)
          enddo
        else
          do nn=1,inne
            read(8) nor(nn,0:3),sp(nn,0:2)
          enddo
        endif
        CLOSE(8)
        write(6,*) 'read end =',input(na1:na2)//'.nsdata'
        nsdata=.true.
!      endif
3001  sp(0,0)=sdum
       sp(inne+1:inne2+1,0)=sdum
end subroutine

      subroutine read_odata_bin(kdc)
      use variables,only:ndri,inputshape
      use arrays
      implicit none
      integer,intent(in)::kdc
      integer::ipmax,i,j,ios1,ios2
      doubleprecision::xmove(1:3)

      open(8,file=inputshape,status='old',form='unformatted',err=201)  !地形ポリゴンの読み込み
!      kdc=1
      if (.not.allocated(rhod)) then
        allocate(rhod(ndri),ivt(ndri),ipl(ndri))
        rhod=0.0d0; ivt=0; ipl=0
      endif
      read(8) ivt(kdc),ipl(kdc),ipmax   !!頂点の数 面の数
      write(6,*) 'kdc=',kdc,'ivt(kdc)=',ivt(kdc),'ipl(kdc)=',ipl(kdc)
      if (.not.allocated(vt)) then
        allocate(vt(ndri,1:ivt(kdc),1:3))
        vt=0.0d0
      endif
      if (.not.allocated(ipo)) then
        allocate(ipo(ndri,1:ipl(kdc)),npo(ndri,1:ipl(kdc),1:ipmax))
        ipo=0; npo=0
      endif
      do i=1,ivt(kdc)
!        read(8) vt(kdc,i,1:3)
        read(8) vt(kdc,i,1),vt(kdc,i,2),vt(kdc,i,3)
!        write(6,*) i, vt(1,i,1:3)
      enddo
      do i=1,ipl(kdc)
        read(8) ipo(kdc,i),(npo(kdc,i,j),j=1,ipo(kdc,i))  !面iのj番目の頂点番号（時計回り）
!        write(6,*) i,ipo(1,i),(npo(1,i,j),j=1,ipo(1,i))  !面iのj番目の頂点番号（時計回り）
      enddo
      read(8,iostat=ios1) rhod(kdc) !密度
      read(8,err=90,iostat=ios2) xmove(1:3)   ! 頂点の平行移動
      close(8)
      if(ios1.ne.0) rhod(kdc)=-400.0d0
      if(ios2.ne.0) xmove=0.0d0
      ! 平行移動がある場合
      if(xmove(1).ne.0.0d0) vt(kdc,1:ivt(kdc),1)=vt(kdc,1:ivt(kdc),1)+xmove(1)
      if(xmove(2).ne.0.0d0) vt(kdc,1:ivt(kdc),2)=vt(kdc,1:ivt(kdc),2)+xmove(2)
      if(xmove(3).ne.0.0d0) vt(kdc,1:ivt(kdc),3)=vt(kdc,1:ivt(kdc),3)+xmove(3)
      return
  201 continue
      write(6,*) 'read err=',inputshape
      stop
   90 continue
      return

      end subroutine read_odata_bin

subroutine read_objectshape()
use arrays,only:rhod,ivt,ipl,vtplus,OBJECTshape,vt,ipo,npo,GD
use variables,only:ndri,iplmax
implicit none
integer::i,j,n1
integer::ipomax(1:ndri),ii,ivtmax,ipomaxmax !,iplmax
allocate(ivt(1:ndri),ipl(1:ndri))
if(.not.allocated(vtplus)) then
    allocate(vtplus(1:ndri,0:2));vtplus=0.0d0
endif
do ii=1,ndri
    n1 = len_trim(OBJECTshape(ii))
    if (OBJECTshape(ii)(n1-4:n1) .ne. 'odata') then
        ivt(ii) = 0; ipl(ii) = 0; ipomax(ii) = 0
    else
        open(8,file=OBJECTshape(ii),status='old',form='unformatted')  !地形ポリゴンの読み込み
        read(8) ivt(ii),ipl(ii),ipomax(ii)   !!頂点の数 面の数     
        close(8)
    endif
enddo        
ivtmax=maxval(ivt(1:ndri)) !+10000
iplmax=maxval(ipl(1:ndri))
ipomaxmax=maxval(ipomax(1:ndri))

        allocate(vt(1:ndri,1:ivtmax,1:3))      
        allocate(ipo(1:ndri,1:iplmax),npo(1:ndri,1:iplmax,1:ipomaxmax)) 
        
do ii=1,ndri
    if (ivt(ii).eq.0) cycle
    if(allocated(rhod)) then
    if(rhod(ii).lt.0.0d0) then
        vtplus(ii,:)=GD(ii,:);GD(ii,:)=0.0d0
        continue
    endif
    endif
    
open(8,file=OBJECTshape(ii),status='old',form='unformatted')  !地形ポリゴンの読み込み

     
        read(8) !  読み飛ばし                     ivt(1),ipl(1),ipmax   !!頂点の数 面の数      

        do i=1,ivt(ii)
            read(8)  vt(ii,i,1:3)
            vt(ii,i,1:3)=vt(ii,i,1:3)+vtplus(ii,0:2)
        enddo
     
        do i=1,ipl(ii)
            read(8) ipo(ii,i),(npo(ii,i,j),j=1,ipo(ii,i))  !面iのj番目の頂点番号（時計回り）           
        enddo
  !      read(8) rhod(ii) !密度   
         
close(8)
enddo
end subroutine read_objectshape

subroutine read_o2data
use variables,only:ndri,inputxyz,inputshape,mm0,mm0mx&
    ,ncpl02max,ncpl01max,ncpl00max,ipl0,title,inns,inne
use arrays,only:objnoi,pp_0,ncpl0,info_0,upl0,objno
!      use object_module
      implicit none
      integer::na1,na2,na3,mm,ic1,ic2,ic3 !,ic1_0,ic2_0,ic3_0

      if(ndri.eq.0) then
        call owari(inputxyz,na2);call hajime(inputxyz,na1)
        do mm=na2,na1,-1
          if(inputxyz(mm:mm).eq.".") then
            na3=mm-1
            exit
          endif
        enddo
        open(9,file=inputxyz(na1:na3)//".o2data",status='old',form='unformatted',err=234)
      else
        call owari(inputshape,na2);call hajime(inputshape,na1)
        open(9,file=inputshape(na1:na2-4)//"2data",status='old',form='unformatted')
      endif

      read(9) mm0
      mm0mx=mm0
      allocate(OBJNOI(1:mm0))

      read(9) OBJNOI(1:mm0)
      read(9) ncpl00max,ncpl02max,ncpl01max
      write(6,*) ncpl00max,ncpl02max,ncpl01max
      allocate(pp_0(1:mm0mx,0:ncpl01max,0:2),ncpl0(mm0mx,0:ncpl00max,-2:ncpl02max),info_0(mm0mx,1:2,1:ncpl01max))

      do mm=1,mm0
          
#ifdef WANKYOKU
        read(9) ic1,ic2,ic3,ncpl0(mm,0:ic1,-2:ic2)   !,info_0(mm,1:2,1:ic3)
#else
        read(9) ic1,ic2,ic3,ncpl0(mm,0:ic1,-2:ic2),info_0(mm,1:2,1:ic3)
#endif
        read(9) pp_0(mm,1:ic3,0:2)
      enddo

      read(9) ipl0
      allocate(upl0(1:ipl0,1:4))

      do mm=1,ipl0
        read(9) upl0(mm,1:4)
      enddo
      close(9)

      allocate(OBJNO(inns:inne));OBJNO=0
      do mm=1,mm0
        if(OBJNOI(mm).le.inne.and.OBJNOI(mm).ge.1) then
            OBJNO(OBJNOI(mm))=mm
        else
            OBJNOI(mm)=0
        endif
      enddo

      if(ndri.eq.0) then
        write(6,*) 'read end ',inputxyz(na1:na3)//".o2data"
      else
        write(6,*) title(na1:na2)//".o2data"
      endif
      return
  234 continue
      mm0=0 ! メッシュ内地形なし
      allocate(OBJNO(1:inne));objno=0
end subroutine read_o2data

subroutine read_dridata
use variables,only:DRIINIT,ndri
use arrays,only:OBJECTshape,rhod,drift,GD
implicit none
integer::ii

      read(5,'(A255)') DRIINIT;print *,'DRIFILE=',DRIINIT
    if(DRIINIT=='NODRI') return
    drift=.true.
    open(55,file=DRIINIT,status='old',action='read',err=2001)
    read(55,*) ndri;print *,'ndri=',ndri
    allocate(OBJECTshape(1:ndri),GD(1:ndri,0:2),rhod(1:ndri))
    do ii=1,ndri
          read(55,'(A255)') OBJECTshape(ii)
          read(55,*) GD(ii,0:2)              !初期重心位置
          read(55,*) rhod(ii)
    enddo 
        call read_objectshape()
      return
2001  continue
      write(6,*) 'Cannot open file',DRIINIT
end subroutine
    
subroutine output_odata_bin(outputshape,kdc)
!use variables,only:ndri
use arrays
implicit none
CHARACTER(LEN=255),intent(in)::outputshape
integer,intent(in)::kdc
integer::ipmax,i,j
doubleprecision::rhod0=-400.0d0
if(.not.allocated(vtplus_S)) then
    allocate(vtplus_s(0:2))
    vtplus_s=0.0d0
endif
ipmax=maxval(ipo(kdc,1:ipl(kdc))) ! 面の頂点数最大値を設定
open(9,file=outputshape,status='unknown',form='unformatted',err=201)  !地形ポリゴンの読み込み
!        kdc=1

        write(9) ivt(kdc),ipl(kdc),ipmax   !!頂点の数 面の数  
        write(6,*) 'kdc=',kdc,'ivt(kdc)=',ivt(kdc),'ipl(kdc)=',ipl(kdc),'ipmax=',ipmax       
        do i=1,ivt(kdc)
            write(9) vt(kdc,i,1:3)
       !     write(6,*) i, vt(1,i,1:3)            
        enddo
     
        do i=1,ipl(kdc)
            write(9) ipo(kdc,i),(npo(kdc,i,j),j=1,ipo(kdc,i))  !面iのj番目の頂点番号（時計回り）
        !    write(6,*) i,ipo(1,i),(npo(1,i,j),j=1,ipo(1,i))  !面iのj番目の頂点番号（時計回り）           
        enddo
!        write(6,*) kdc,npo(kdc,1101,1:ipo(kdc,1101))
        if(allocated(rhod)) then 
            write(9) rhod(kdc) !密度   
        else
            write(9) rhod0
        endif
        write(9) vtplus_s(0:2)   ! 頂点の平行移動         
close(9)
!平行移動がある場合
return
201 continue
    write(6,*) 'read err=',outputshape
    stop
end subroutine output_odata_bin

subroutine output_xyz
USE variables,only:nx1,nx2,title
USE arrays,only:ls,le,x
implicit none
integer::i,k
OPEN(9,FILE=title(nx1:nx2)//'.xyz',form='unformatted', STATUS= 'UNKNOWN') 
       write(9) ls(0:2),le(0:2)
       write(6,*) ls(0:2),le(0:2)
       do k=0,2
       write(9) (X(k,I),I=ls(k)-2,le(k)+1)
       enddo
      CLOSE(9)    
end subroutine  output_xyz    

subroutine output_o2data(flag)
use variables,only:inputxyz,mm0,title
use arrays,only:objnoi,ncpl0,info_0,ipl,pp_0,upl
implicit none
integer::na1,na2,mm,ic1,ic2,ic3,flag

if(flag.eq.0) then  !地形変化なし
 call owari(inputxyz,na2);call hajime(inputxyz,na1)
open(9,file=inputxyz(na1:na2-4)//".o2data",status='unknown',form='unformatted') !時間的に移動しない
else
call owari(title,na2);call hajime(title,na1)
open(9,file=title(na1:na2)//".o2data",status='unknown',form='unformatted')
endif
write(9) mm0
write(9) OBJNOI(1:mm0)

ic1=maxval(ncpl0(1:mm0,0,0))  !セルあたりの面の最大数
ic3=maxval(ncpl0(1:mm0,0,1))  !セルあたりの頂点の最大数
!ic2=0
!do mm=1,mm0
!ic2=max(ic2, maxval(ncpl0(mm,1:ncpl0(mm,0,0),0)) )*2   !セルあたりの面の頂点の最大数
!enddo
ic2=maxval(ncpl0(1:mm0,0,2))
write(9) ic1,ic2,ic3
write(6,*) ic1,ic2,ic3
do mm=1,mm0
ic1=ncpl0(mm,0,0)
ic2=ncpl0(mm,0,2)
ic3=ncpl0(mm,0,1)
#ifdef WANKYOKU
if(mm.eq.10330) then
continue
endif
write(9) ic1,ic2,ic3,ncpl0(mm,0:ic1,-2:ic2)
#else
write(9) ic1,ic2,ic3,ncpl0(mm,0:ic1,-2:ic2),info_0(mm,1:2,1:ic3)
#endif
write(9) pp_0(mm,1:ic3,0:2)
enddo
!deallocate(pp_0)

write(9) ipl(1)
do mm=1,ipl(1)
write(9) upl(1,mm,1:4)
enddo
deallocate(upl)
close(9)
if(flag.eq.0) then
write(6,*) inputxyz(na1:na2-4)//".o2data"  !地形変化なし
else
write(6,*) title(na1:na2)//".o2data"  !地形変化あり
endif
end subroutine output_o2data           
!
subroutine output_data
use variables,only:title,nx1,nx2,is,js,ks,ie,ke,je,fg_xyz,g,nu,t,inne,inne2,dt
use arrays,only:in,jn,kn,u,f,p,xyzn,nff,mn
implicit none
integer::nn,i,k,j,icc,icc2,fnum
integer,allocatable::nfb1(:,:,:)
fnum=9
call owari(title,nx2);call hajime(title,nx1)
  !    OPEN(newunit=fnum,FILE=title(nx1:nx2)//'.txt', STATUS= 'UNKNOWN')
      OPEN(newunit=fnum,FILE=title(nx1:nx2),form='unformatted', STATUS= 'UNKNOWN')
      write(fnum) T,g,nu
  !    write(fnum,*) T,g,nu
      if(.not.xyzn) then
        write(fnum) (((NFF(mn(I,k,j))%f,K=ks-1,KE),J=js-1,JE),I=is-1,IE)
        write(6,*) 'pass nf'
        allocate(nfb1(is-1:ie,ks-1:ke,js-1:JE))
        do i=is-1,ie
            do k=ks-1,ke
                do j=js-1,je
                    nfb1(i,k,j)=nff(mn(i,k,j))%b
                enddo
            enddo
        enddo
       ! nfb1(is-1:ie,ks-1:ke,js-1:JE)=nfb(is-1:ie,ks-1:ke,js-1:JE)
!$omp parallel do
        do nn=inne+1,inne2
          if(nff(nn)%b.eq.-1) then
            nfb1(in(nn),kn(nn),jn(nn))=300
          else
            nfb1(in(nn),kn(nn),jn(nn))=200
          endif
        enddo
!$omp end parallel do
        write(fnum) (((NFB1(I,k,j),K=ks-1,KE),J=js-1,JE),I=is-1,IE)
        write(6,*) 'pass nfb'
      else
        icc=0;icc2=0
        do nn=1,inne2
          if(nff(nn)%f.eq.0.and.sum(u(0:2,nn)**2).eq.0.0d0.and.f(nn).lt.1.0d-13) cycle
          if(nn.le.inne) then
            icc=icc+1
          else
            icc2=icc2+1
          endif
        enddo
        write(fnum) icc,icc2
  !      write(fnum,*) icc,icc2
        write(6,*) 'icc,icc2(in output_data)=',icc,icc2        
      endif

      if (xyzn) then
        do nn=1,inne2
          if(nfF(nn)%f.eq.0.and.sum(u(0:2,nn)**2).eq.0.0d0.and.f(nn).lt.1.0d-13) cycle
          
          if(nn.le.inne) then
            write(fnum) nn,U(0:2,nn),F(nn),P(nn),nff(nn)%f,nff(nn)%b
   !         write(fnum,888) nn,U(0:2,nn),F(nn),P(nn),nff(nn)%f,nff(nn)%b
          else
            write(fnum) nn,U(0:2,nn)
          endif
        enddo
      else
        do nn=1,inne2
          if(nn.le.inne) then
            write(fnum) U(0:2,nn),F(nn),P(nn)
          else
            write(fnum) U(0:2,nn)
          endif
        enddo
      endif
      write(fnum) dt
      close(fnum)

      write(6,*) 'write end =',title(nx1:nx2)
888 format(1h ,i8,5(',',e12.5),',',i3,',',i3)

end subroutine output_data  
!
subroutine output_fdata
use variables,only:inne,inne2,title
use arrays,only:a,fb
implicit none
integer::nn,na1,na2,icc,icc2  !,icc=0
call owari(title,na2);call hajime(title,na1)
icc=0;icc2=0
    do nn=1,inne
        if(a(0,nn).eq.1.0d0.and.a(1,nn).eq.1.0d0.and.a(2,nn).eq.1.0d0.and. fb(nn).eq.1.0d0)  cycle
        icc=icc+1
    enddo
    do nn=inne+1,inne2    
        if(a(0,nn)+a(1,nn)+a(2,nn).eq.0.0d0) cycle
        icc2=icc2+1
 !           write(9) a(0:2,nn)  ! 必要に応じてコメントアウト
    enddo
open(9,file=title(na1:na2)//".fdata",status='unknown',form='unformatted')
     write(9) icc,icc2
     write(6,*) 'icc,icc2=',icc,icc2
    do nn=1,inne2
#ifdef BUGWIN
          if(nn.eq.400932) then
              continue
          endif
#endif
        if(nn.le.inne) then
        if(a(0,nn).eq.1.0d0.and.a(1,nn).eq.1.0d0.and.a(2,nn).eq.1.0d0.and.fb(nn).eq.1.0d0) cycle
            write(9) nn, a(0:2,nn),fb(nn)
        else
            if(a(0,nn)+a(1,nn)+a(2,nn).eq.0.0d0) cycle
            write(9) nn, a(0:2,nn)  ! 必要に応じてコメントアウト 
        endif
    enddo
close(9)
write(6,*) 'write fdata',title(na1:na2)//".fdata", 'icc,icc2=',icc,icc2
end subroutine output_fdata
!
subroutine output_godata
use variables,only:nx1,nx2,title
use arrays,only:mmd,nd,dp,mmd2,nd2,dp2,idp2
implicit none
integer::i,j
OPEN(9,FILE=title(nx1:nx2)//'.3d.godata',form='unformatted',STATUS= 'UNKNOWN')
        write(9) mmd(1),maxval(nd(1,1:mmd(1)))       
        do i=1,mmd(1)  ! mmd:面の総数
          write(9) nd(1,i),(dp(1,i,j,:),j=1,nd(1,i))
        enddo
CLOSE(9)        
       write(6,*) 'mmd,nd_max=',mmd(1),maxval(nd(1,1:mmd(1))) 
OPEN(9,FILE=title(nx1:nx2)//'.2d.godata',form='unformatted',STATUS= 'UNKNOWN')           
        write(9) mmd2(1),maxval(nd2(1,1:mmd2(1)))
        do i=1,mmd2(1)
          write(9) (IDP2(1,i,j),j=1,2),nd2(1,i),(dp2(1,i,j,:),j=1,nd2(1,i))
        enddo
CLOSE(9)        
      write(6,*) 'mmd2,nd2_max=',mmd2(1),maxval(nd2(1,1:mmd2(1))) 
      write(6,*) 'end godata'
end subroutine output_godata
!
subroutine output_drdata
use variables,only:nx1,nx2,title,t,ndri
use arrays,only:mmd,nd,dp,nd,mmd2,idp2,nd2,dp2
implicit none
integer::ii,i,j
call owari(title,nx2);call hajime(title,nx1)
#ifndef HENDO
if(ubound(nd,1).eq.0) return
OPEN(9,FILE=title(nx1:nx2)//'.3d.drdata',form='unformatted',STATUS= 'UNKNOWN')
write(9) t
write(9) ndri
write(9) mmd(1:ndri) 
do ii=1,ndri
        do i=1,mmd(ii)  ! mmd:面の総数
          write(9) nd(ii,i),(dp(ii,i,j,:),j=1,nd(ii,i))
        enddo
enddo        
CLOSE(9)
#endif
!if(ubound(nd2,1).eq.1) return
OPEN(9,FILE=title(nx1:nx2)//'.2d.drdata',form='unformatted',STATUS= 'UNKNOWN')
write(9) ndri
do ii=1,ndri    
        write(9) mmd2(ii)
        do i=1,mmd2(ii)
          write(9) (IDP2(ii,i,j),j=1,2),nd2(ii,i),(dp2(ii,i,j,:),j=1,nd2(ii,i))          
        enddo
enddo     
CLOSE(9)
end subroutine output_drdata
! 
subroutine output_nsdata
use variables,only:isfn,fg_xyz,mms,inne,nx1,nx2,title
use arrays,only:sfn,sfnoi,sp,nor,sfno2,sfno2i,surfp
      implicit none
      INTEGER::nn,i,fnum

      OPEN(newunit=fnum,FILE=title(nx1:nx2)//'.nsdata',form='unformatted', STATUS= 'UNKNOWN')
      if(fg_xyz.eq.1) then
        write(fnum) mms
        do i=1,mms  ! mms:水面セルの総数
          nn=SFNOI(i)
          write(fnum) nn,nor(nn,0:3),sp(nn,0:2)
        enddo
      else
        do nn=1,inne
          write(fnum) nor(nn,0:3),sp(nn,0:2)
        enddo
      endif
      CLOSE(fnum)
!#ifdef BUGWIN
      write(6,*) 'write end =',title(nx1:nx2)//'.nsdata mms=',mms
!#endif>	naga02.exe!OUTPUT_XYZ()  行 106	Fortran
end subroutine output_nsdata
subroutine output_sfdata
      use variables,only:isfn,inne,nx1,nx2,title,t,mms2  !,fg_xyz
      use arrays,only:sfn,f,in,jn,kn,nds,surfp,SFNO2I,dim3d,nff,mn
      implicit none
      integer::mms3
      integer::i,j,i0,nn,fnum !,k,ii,mm
      ! 水面の頂点座標（3次元可視化用）
#ifdef FLAT
      doubleprecision::pp1(1:4,1:3)
      integer::i4=4
#endif
    if(.not.dim3d) return
      i0=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!追加
      write(6,*) 'mms2=',mms2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!追加
      mms3=0
#ifdef FLAT
      do mm=1,isfn
        nn=sfn(mm)
        if(fb(nn).ne.1.0d0) cycle
        i=in(nn);j=jn(nn);k=kn(nn)
        if(nf(i,k,j).ne.2) cycle
        mms3=mms3+1
      enddo
#endif
      OPEN(newunit=fnum,FILE=title(nx1:nx2)//'.sfdata',form='unformatted',STATUS= 'UNKNOWN')
      write(fnum) t
       write(6,*) 'ndsmax=',maxval(nds(1:mms2)),ubound(surfp,2)     
      write(fnum) mms2+mms3,maxval(nds(1:mms2))
      do i=1,mms2  ! mms2:水面の総数mms2
        nn=SFNO2I(i)
        if(nn.eq.0) then
          write(fnum) i0              
        else if(f(nn).lt.1.0d-2.and.nff(mn(in(nn),kn(nn)-1,jn(nn)))%f.eq.-1) then
          write(fnum) i0    !nds(i),(surfp(i,j,:),j=1,nds(i))
        else
          write(fnum) nds(i),(surfp(i,j,:),j=1,nds(i))
        endif
      enddo
#ifdef FLAT
      do mm=1,isfn
        nn=sfn(mm)
        if(fb(nn).ne.1.0d0) cycle
        i=in(nn);j=jn(nn);k=kn(nn)
        if(nf(i,k,j).ne.2) cycle
        if(nfb(i,k,j).eq.-3) then
          pp1(1:4,3)=z(k)+(z(k+1)-z(k))*f(nn)
          pp1(1:4,1)=(/x(i),x(i+1),x(i+1),x(i)/)
          pp1(1:4,2)=(/y(j),y(j),y(j+1),y(j+1)/)
        else if(nfb(i,k,j).eq.3) then
          pp1(1:4,3)=z(k+1)-(z(k+1)-z(k))*f(nn)
          pp1(1:4,1)=(/x(i),x(i),x(i+1),x(i+1)/)
          pp1(1:4,2)=(/y(j),y(j+1),y(j+1),y(j)/)
        else if(nfb(i,k,j).eq.-2) then
          pp1(1:4,2)=y(j)+(y(j+1)-y(j))*f(nn)
          pp1(1:4,1)=(/x(i),x(i),x(i+1),x(i+1)/)
          pp1(1:4,3)=(/z(k),z(k+1),z(k+1),z(k)/)
        else if(nfb(i,k,j).eq.2) then
          pp1(1:4,2)=y(j+1)-(y(j+1)-y(j))*f(nn)
          pp1(1:4,1)=(/x(i),x(i+1),x(i+1),x(i)/)
          pp1(1:4,3)=(/z(k),z(k),z(k+1),z(k+1)/)
        else if(nfb(i,k,j).eq.-1) then
          pp1(1:4,1)=x(i)+(x(i+1)-x(i))*f(nn)
          pp1(1:4,2)=(/y(j),y(j+1),y(j+1),y(j)/)
          pp1(1:4,3)=(/z(k),z(k),z(k+1),z(k+1)/)
        else if(nfb(i,k,j).eq.1) then
          pp1(1:4,1)=x(i+1)-(x(i+1)-x(i))*f(nn)
          pp1(1:4,2)=(/y(j),y(j),y(j+1),y(j+1)/)
          pp1(1:4,3)=(/z(k),z(k+1),z(k+1),z(k)/)
        endif
        write(fnum) i4,(pp1(j,:),j=1,i4)
      enddo
#endif
!#ifdef BUGWIN
      write(6,*) 'write end =',title(nx1:nx2)//'.sfdata'
!#endif
      CLOSE(fnum)
    end subroutine output_sfdata 
 subroutine read_xyz_cnv
 use variables,only:is,js,ks,ie,je,ke,inputxyz,nx1,nx2
 use arrays,only:xi,y,z
 implicit none
 integer::i,j,k
 call owari(inputxyz,nx2);call hajime(inputxyz,nx1)
OPEN(9,FILE=inputxyz(nx1:nx2),form='unformatted')
read(9) is,js,ks,IE,JE,KE
write(6,*) is,js,ks,IE,JE,KE
if(ubound(xi,1).eq.0) allocate(xi(1:ie+1))
read(9) (Xi(I),I=is-2,IE+1)
if(ubound(y,1).eq.0) allocate(y(1:je+1))
read(9) (Y(J),J=js-2,JE+1)
if(ubound(z,1).eq.0) allocate(z(1:ke+1))
read(9) (Z(K),K=ks-2,KE+1)
CLOSE(9) 
end subroutine

subroutine read_main_cnv
 use variables,only:is,js,ks,ie,je,ke,input,nx1,nx2,t,g,nu,dt0
 use arrays,only:nf3,nfb3,ui,v,w,f3,p3
 implicit none
 integer::i,j,k,iccc
 doubleprecision::UU
 call owari(input,nx2);call hajime(input,nx1)
OPEN(8,FILE=input(nx1:nx2),form='unformatted', STATUS= 'UNKNOWN')       
read(8) T,g,nu
write(6,*) T,g,nu

do j=js,je
do k=ks,ke
do i=is,ie
if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).eq.0) then
!  if((nf(i-1,k,j).ge.0.and.fb(i-1,k,j).lt.1.0d0)&
!.or.(nf(i,k-1,j).ge.0.and.fb(i,k-1,j).lt.1.0d0)&
!.or.(nf(i,k,j-1).ge.0.and.fb(i,k,j-1).lt.1.0d0)&
!) then

if(nf3(i-1,k,j).ge.0.or.nf3(i,k-1,j).ge.0.or.nf3(i,k,j-1).ge.0 ) then
!write(6,*) 'nfb=200',i,k,j,nf(i,k,j),nf(i,k,j-1)
nfb3(i,k,j)=200
endif

endif

enddo
enddo
enddo

read(8) (((NF3(I,k,j),K=ks-1,KE),J=js-1,JE),I=is-1,IE)
read(8) (((NFB3(I,k,j),K=ks-1,KE),J=js-1,JE),I=is-1,IE)

iccc=0
UU=0.0d0
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).ge.0.or.(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1.and.nfb3(i,k,j).lt.200) ) then
read(8) ui(i,k,j),v(i,k,j),w(i,k,j),F3(i,k,j),p3(i,k,j)
iccc=iccc+1
endif
enddo
enddo
enddo
write(6,*) 'iccc=',iccc
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).eq.200) then
read(8) ui(i,k,j),v(i,k,j),w(i,k,j)
iccc=iccc+1
endif
enddo
enddo
enddo

      if(t.gt.0.0d0) then
       read(8) dt0
       write(6,*) 'read dt0=',dt0   
      endif

CLOSE(8)
write(6,*) 'iccc2=',iccc
end subroutine
   
subroutine read_fdata_cnv
 use variables,only:is,js,ks,ie,je,ke,inputf,nx1,nx2,t,g,nu,dt0
 use arrays,only:nf3,nfb3,ax,ay,az,fb3
 implicit none
 integer::i,j,k !,iccc
! doubleprecision::UU
 call owari(inputf,nx2);call hajime(inputf,nx1)
OPEN(8,FILE=inputf(nx1:nx2),form='unformatted', STATUS= 'UNKNOWN')       

do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).ge.0.or.(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1 .and.nfb3(i,k,j).lt.200) ) then
    read(8) ax(i,k,j),ay(i,k,j),az(i,k,j),fb3(i,k,j)
endif
enddo
enddo
enddo
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.200) then
read(8) ax(i,k,j),ay(i,k,j),az(i,k,j)
endif
enddo
enddo
enddo  
write(6,*) 'write end =',inputf(nx1:nx2)
CLOSE(8)
end subroutine
    
subroutine write_xyz
use variables,only:is,js,ks,IE,JE,KE,nx1,nx2,title
use arrays,only:xi,y,z
implicit none
integer::i,j,k
call owari(title,nx2);call hajime(title,nx1)
OPEN(9,FILE=title(nx1:nx2)//'.xyz',form='unformatted')
write(9) is,js,ks,IE,JE,KE
write(6,*) is,js,ks,IE,JE,KE
write(9) (Xi(I),I=is-2,IE+1)
write(9) (Y(J),J=js-2,JE+1)
write(9) (Z(K),K=ks-2,KE+1)
CLOSE(9) 
end subroutine

subroutine write_main
use variables,only:is,js,ks,IE,JE,KE,nx1,nx2,title,T,g,nu
use arrays,only: nf3,nfb3,f3,fb3
implicit none
doubleprecision::UU,PP
integer::i,k,j,iccc
call owari(title,nx2);call hajime(title,nx1)
OPEN(9,FILE=title(nx1:nx2),form='unformatted', STATUS= 'UNKNOWN')       
write(9) T,g,nu
write(6,*) T,g,nu

do j=js,je
do k=ks,ke
do i=is,ie
if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).eq.0) then
!  if((nf(i-1,k,j).ge.0.and.fb(i-1,k,j).lt.1.0d0)&
!.or.(nf(i,k-1,j).ge.0.and.fb(i,k-1,j).lt.1.0d0)&
!.or.(nf(i,k,j-1).ge.0.and.fb(i,k,j-1).lt.1.0d0)&
!) then

if(nf3(i-1,k,j).ge.0.or.nf3(i,k-1,j).ge.0.or.nf3(i,k,j-1).ge.0 ) then
!write(6,*) 'nfb=200',i,k,j,nf(i,k,j),nf(i,k,j-1)
nfb3(i,k,j)=200
endif

endif

enddo
enddo
enddo

write(9) (((NF3(I,k,j),K=ks-1,KE),J=js-1,JE),I=is-1,IE)
write(9) (((NFB3(I,k,j),K=ks-1,KE),J=js-1,JE),I=is-1,IE)

iccc=0
UU=0.0d0;PP=0.0d0
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).ge.0.or.(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1.and.nfb3(i,k,j).lt.200) ) then
!if (ks == ke-1) then
!    if (F3(i,k,j).eq.0.0d0.and.nf(i,k,j).ge.0) then
!        PP=(1.0d0-fb3(i,k,j))*(z(k+1)-z(k))+z(k)
!    else
!        PP=0.0d0
!    endif
!endif
write(9) UU,UU,UU,F3(i,k,j),PP
iccc=iccc+1
endif
enddo
enddo
enddo
write(6,*) 'iccc=',iccc
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).eq.200) then
write(9) UU,UU,UU
iccc=iccc+1
endif
enddo
enddo
enddo

CLOSE(9)
write(6,*) 'iccc2=',iccc
end subroutine

subroutine write_fdata3
use variables,only:is,js,ks,IE,JE,KE,title,T,g,nu,inne,inne2
use arrays,only: nf3,nfb3,fb3,ax,ay,az,xyzn,in,jn,kn
implicit none
integer::i,k,j,nx1,nx2,nn,icc,icc2
call owari(title,nx2);call hajime(title,nx1)
OPEN(9,FILE=title(nx1:nx2)//'.fdata',form='unformatted', STATUS= 'UNKNOWN')       
if(xyzn) then
icc=0;icc2=0
    do nn=1,inne
        i=in(nn);j=jn(nn);k=kn(nn)
        if(ax(i,k,j).eq.1.0d0.and.ay(i,k,j).eq.1.0d0.and.az(i,k,j).eq.1.0d0.and. fb3(i,k,j).eq.1.0d0)  cycle
        icc=icc+1
    enddo
    do nn=inne+1,inne2
        i=in(nn);j=jn(nn);k=kn(nn)
        if(ax(i,k,j)+ay(i,k,j)+az(i,k,j).eq.0.0d0) cycle
        icc2=icc2+1
 !           write(9) a(0:2,nn)  ! 必要に応じてコメントアウト
    enddo
     write(9) icc,icc2
     write(6,*) icc,icc2,inne,inne2
    do nn=1,inne2
        i=in(nn);j=jn(nn);k=kn(nn)
        if(nn.le.inne) then
        if(ax(i,k,j).eq.1.0d0.and.ay(i,k,j).eq.1.0d0.and.az(i,k,j).eq.1.0d0.and. fb3(i,k,j).eq.1.0d0)  cycle
            write(9) nn,ax(i,k,j),ay(i,k,j),az(i,k,j),fb3(i,k,j)
        else  
        if(ax(i,k,j)+ay(i,k,j)+az(i,k,j).eq.0.0d0) cycle
            write(9) nn,ax(i,k,j),ay(i,k,j),az(i,k,j)! 必要に応じてコメントアウト
        endif
    enddo
else
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).ge.0.or.(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1.and.nfb3(i,k,j).lt.200)) then
write(9) ax(i,k,j),ay(i,k,j),az(i,k,j),fb3(i,k,j)
endif
enddo
enddo
enddo
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.200) then
write(9) ax(i,k,j),ay(i,k,j),az(i,k,j)
endif
enddo
enddo
enddo  
endif
write(6,*) 'write end =',title(nx1:nx2)//'.fdata'
CLOSE(9)
end subroutine write_fdata3
subroutine read_xyz3
use variables,only:inputxyz,is,js,ks,ie,je,ke
use arrays,only:xi,y,z
implicit none
integer::na1,na2,i,j,k
call owari(inputxyz,na2);call hajime(inputxyz,na1)
         write(6,*) '[HOST]read STdata',na1,inputxyz(na1:na2)
      OPEN(8,FILE=inputxyz(na1:na2),status='old',form='unformatted')
       read(8) is,js,ks,IE,je,KE
       write(6,*) is,js,ks,IE,JE,KE
       write(6,*) 'size(x)=',size(xi),' ie=',ie
       allocate(xi(1:ie+2),y(1:je+2),z(1:ke+2))
       read(8) (xi(i),I=is-2,IE+1)
       read(8) (Y(J),J=js-2,je+1)
       read(8) (Z(K),K=ks-2,KE+1)
       continue
CLOSE(8)   
end subroutine read_xyz3

subroutine read_main3
use variables,only:is,js,ks,IE,JE,KE,nx1,nx2,input,T,g,nu,inns,inne
use arrays,only: nf3,nfb3,f3,fb3,z,in,jn,kn
implicit none
!doubleprecision::UU,PP
integer::i,k,j,iccc,inn2,inn3,inne2,inn

call owari(input,nx2);call hajime(input,nx1)
OPEN(8,FILE=input(nx1:nx2),form='unformatted', STATUS= 'old')       
read(8) T,g,nu
write(6,*) T,g,nu
allocate(nf3(is-2:ie,ks-2:ke,js-2:je),nfb3(is-2:ie,ks-2:ke,js-2:je))
read(8) (((NF3(I,k,j),K=ks-1,KE),J=js-1,JE),I=is-1,IE)
read(8) (((NFB3(I,k,j),K=ks-1,KE),J=js-1,JE),I=is-1,IE)

  139 do j=js-1,je;do k=ks-1,ke;do i=is-1,ie
      if(nf3(i,k,j).ge.0.or.&
         (nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1&
          .and.nfb3(i,k,j).lt.200)) then
        if(nf3(i,k,j).eq.2) inn2=inn2+1
        if(nf3(i,k,j).eq.-1) inn3=inn3+1
        inn=inn+1;in(inn)=i;jn(inn)=j;kn(inn)=k !mn(i,k,j)=inn;
      else
        ! mn(i,k,j)=0
      end if
      enddo;enddo;enddo
      inns=1;inne=inn

      do j=js-1,je;do k=ks-1,ke;do i=is-1,ie
        if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.200) then
          inn=inn+1 ;in(inn)=i;jn(inn)=j;kn(inn)=k ! ;mn(i,k,j)=inn
          if(nfb3(i,k,j).eq.300) then
            nfb3(i,k,j)=-1
          else
            nfb3(i,k,j)=0
          endif
      !    in(inn)=i;jn(inn)=j;kn(inn)=k
        end if
      enddo;enddo;enddo
      inne2=inn


#ifdef DDDDD
do j=js,je
do k=ks,ke
do i=is,ie
if(nf(i,k,j).eq.-1.and.nfb(i,k,j).eq.0) then
if(nf(i-1,k,j).ge.0.or.nf(i,k-1,j).ge.0.or.nf(i,k,j-1).ge.0 ) then
nfb(i,k,j)=200
endif
endif
enddo
enddo
enddo


iccc=0
UU=0.0d0;PP=0.0d0
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf(i,k,j).ge.0.or.(nf(i,k,j).eq.-1.and.nfb(i,k,j).ge.1.and.nfb(i,k,j).lt.200) ) then
if (ks == ke-1) then
    if (F3(i,k,j).eq.0.0d0.and.nf(i,k,j).ge.0) then
        PP=(1.0d0-fb3(i,k,j))*(z(k+1)-z(k))+z(k)
    else
        PP=0.0d0
    endif
endif
write(9) UU,UU,UU,F3(i,k,j),PP
iccc=iccc+1
endif
enddo
enddo
enddo
write(6,*) 'iccc=',iccc
do j=js-1,je
do k=ks-1,ke
do i=is-1,ie
if(nf(i,k,j).eq.-1.and.nfb(i,k,j).eq.200) then
write(9) UU,UU,UU
iccc=iccc+1
endif
enddo
enddo
enddo
#endif
CLOSE(9)
write(6,*) 'iccc2=',iccc
end subroutine read_main3

subroutine read_fdata3
use variables,only:input,is,js,ks,ie,je,ke,inputf0
use arrays,only:ax,ay,az,fb3,nff,in,jn,kn,xyzn,mn,nf3,nfb3
implicit none
integer::nx1,nx2,i,j,k,icc,icc1,ii,nn
call owari(inputf0,nx2);call hajime(inputf0,nx1)
if(nx2.gt.1) then
OPEN(8,FILE=inputf0(nx1:nx2),form='unformatted', STATUS= 'old')  
write(6,*) 'read start ',inputf0(nx1:nx2)
else
call owari(input,nx2);call hajime(input,nx1)
OPEN(8,FILE=input(nx1:nx2)//".fdata",form='unformatted', STATUS= 'old') 
write(6,*) 'read start ',input(nx1:nx2)//".fdata"
endif
allocate(ax(is-1:ie,ks-1:ke,js-1:je),ay(is-1:ie,ks-1:ke,js-1:je),az(is-1:ie,ks-1:ke,js-1:je),fb3(is-1:ie,ks-1:ke,js-1:je))

if(xyzn) then
    ax=0.0d0;ay=0.0d0;az=0.0d0;fb3=0.0d0
    do j=js-1,je
    do k=ks-1,ke
    do i=is-1,ie
        nn=mn(i,k,j)
    if(nff(nn)%f.ge.0.or.(nff(nn)%f.eq.-1.and.nff(nn)%b.ge.1 .and.nff(nn)%b.lt.200) ) then
    ax(i,k,j)=1.0d0;ay(i,k,j)=1.0d0;az(i,k,j)=1.0d0;fb3(i,k,j)=1.0d0
    endif
    enddo
    enddo
    enddo
  !  do j=js-1,je
  !  do k=ks-1,ke
  !  do i=is-1,ie
  !      nn=mn(i,k,j)
  !  if(nff(nn)%f.eq.-1.and.nff(nn)%b.ge.200) then
  !  ax(i,k,j)=0.0d0;ay(i,k,j)=0.0d0;az(i,k,j)=0.0d0
  !  endif
  !  enddo
  !  enddo
  !  enddo 

    read(8) icc,icc1
    write(6,*) icc,icc1
    do ii=1,icc 
    read(8) nn, ax(in(nn),kn(nn),jn(nn)),ay(in(nn),kn(nn),jn(nn)),az(in(nn),kn(nn),jn(nn)),fb3(in(nn),kn(nn),jn(nn))
    enddo
    do ii=1,icc1
    read(8) nn,ax(in(nn),kn(nn),jn(nn)),ay(in(nn),kn(nn),jn(nn)),az(in(nn),kn(nn),jn(nn))
    enddo
else
    do j=js-1,je
    do k=ks-1,ke
    do i=is-1,ie
    if(nf3(i,k,j).ge.0.or.(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.1 .and.nfb3(i,k,j).lt.200) ) then
        read(8) ax(i,k,j),ay(i,k,j),az(i,k,j),fb3(i,k,j)
    endif
    enddo
    enddo
    enddo
    do j=js-1,je
    do k=ks-1,ke
    do i=is-1,ie
    if(nf3(i,k,j).eq.-1.and.nfb3(i,k,j).ge.200) then
    read(8) ax(i,k,j),ay(i,k,j),az(i,k,j)
    endif
    enddo
    enddo
    enddo  
endif
write(6,*) 'read end fdata'
CLOSE(8)
end subroutine

subroutine output_rdata
    use variables, only: title, nx1, nx2, fg_xyz, inne
    use arrays, only: nff, pro, r, in, jn, kn, d
    implicit none
    integer :: nn,fnum

    OPEN(newunit=fnum,FILE=title(nx1:nx2)//'.rdata',form='unformatted', STATUS= 'UNKNOWN')
    if (fg_xyz.eq.1) then
        !New Format
        do nn = 1,inne
            if(nff(nn)%f.ne.0) write(fnum) r(nn)%k,r(nn)%e,d(nn),pro(nn)
        enddo
    else
        !Old Format
        do nn = 1,inne
            write(fnum) r(nn)%k,r(nn)%e
        enddo
        do nn = 1,inne
            write(fnum) d(nn)
        enddo
        write(fnum) pro
    endif
    CLOSE(fnum)
    write(6,*) 'write end =',title(nx1:nx2)//'.rdata'
end subroutine output_rdata

!#ifdef DAKU
subroutine output_cdata
use arrays,only:te,c,rho
use variables,only:title,nx1,nx2,mfk1,inne
      implicit none
      INTEGER::nn,fnum

      OPEN(newunit=fnum,FILE=title(nx1:nx2)//'.cdata',form='unformatted', STATUS= 'UNKNOWN')
do nn=1,inne
#ifdef SALT
    write(fnum) te(nn),c(0,nn),c(1,nn)
#else
    write(fnum) te(nn),c(0,nn),rho(nn)
#endif    
enddo

if(mfk1.ge.2) then
do nn=1,inne
    write(fnum) c(1:mfk1,nn)
enddo
endif
      close(fnum)
!#ifdef BUGWIN
      write(6,*) 'write end =',title(nx1:nx2)//'.cdata'
!#endif
    end subroutine output_cdata
!#endif
    
#ifdef DAKU
subroutine read_cdata
use variables,only:inne2,inputtc,t0,mfk1,inne
use arrays,only:c,rho,te,alloc_c,alloc_te,alloc_rho
implicit none   
integer::nn,na1,na2
  !     USE param
  !     USE daku1
  !     USE mitsudo
  !     use mainpara
  !     use dxdydz

call alloc_c(mfk1,inne2)
call alloc_te(inne2)
call alloc_rho(inne2) 



if(t0.eq.0.0d0) return
call owari(inputtc,na2);call hajime(inputtc,na1)
write(6,*) '[HOST]read Cdata',INPUTtc(na1:na2)
call flush(6)
OPEN(8,FILE=INPUTtc(na1:na2),status='old',form='unformatted')

do nn=1,inne  !all
 read(8) te(nn),c(0,nn),rho(nn)
enddo

if(mfk1.ge.2.and.t0.ne.0.0d0) then
 do nn=1,inne  !all
  read(8) c(1:mfk1,nn)
 enddo
endif

CLOSE(8)            
write(6,*) '[HOST]end Cdata',INPUTtc(1:na1)
end subroutine
#endif

subroutine output_xyzn(iflag)
use variables,only:is,js,ks,IE,JE,KE,nx1,nx2,title
use arrays,only:x,nf3,nfb3
implicit none
integer::i,k,j,icc
integer,intent(in)::iflag

if(iflag.eq.1) then
    do j=js,je;do k=ks,ke;do i=is,ie
    if(k.eq.ke.and.ke-1.eq.ks) cycle        
    if(nf3(i,k,j).ne.-1) cycle
    if(nfb3(i,k,j).eq.0) then
    if(nf3(i-1,k,j).ge.0.or.nf3(i,k,j-1).ge.0.or.nf3(i,k-1,j).ge.0) then
 !   if(nf3(i-1,k,j).ge.0.or.nf3(i,k,j-1).ge.0) then
           nfb3(i,k,j)=200;cycle
    endif
    endif   
   ! if(nfb3(i,k,j).eq.-1) then  
      !  else if(nfb3(i,k,j).eq.-1) then
      !    nfb3(i,k,j)=0
      !      write(6,*) '300'
   ! endif

    enddo;enddo;enddo
endif    
    
call owari(title,nx2);call hajime(title,nx1)
OPEN(9,FILE=title(nx1:nx2)//'.xyzn',status='unknown',form='unformatted')
            write(9) is,js,ks,IE,je,KE
            write(6,*) is,js,ks,IE,JE,KE
            write(9) (X(0,I),I=is-2,IE+1)
            write(9) (X(1,J),J=js-2,je+1)
            write(9) (X(2,K),K=ks-2,KE+1)
            icc=0
            do i=is-1,ie
              do j=js-1,je
                do k=ks-1,ke
                  if(nf3(i,k,j).eq.-1) icc=icc+1
                enddo
              enddo
            enddo 
            write(9) icc    
            write(6,*) 'num of nf=',icc
            do i=is-1,ie
              do j=js-1,je
                do k=ks-1,ke
                  if(nf3(i,k,j).eq.-1) then
                    write(9) i,k,j,nf3(i,k,j),nfb3(i,k,j)
            !         write(99,*) i,k,j,nf3(i,k,j),nfb3(i,k,j)
                  endif
                enddo
              enddo
            enddo
          close(9)
end subroutine output_xyzn
