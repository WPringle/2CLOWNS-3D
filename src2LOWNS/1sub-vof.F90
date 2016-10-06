!%%%%%%%%%%%%%%%%%%%% FILE: 1sub-vof.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Updates the F function from the advection equation using new velocities 
#if defined(DRI) || defined(HENDO)
#define A_ a2
#define FB_ fbd
#else
#define A_ a
#define FB_ fb
#endif
subroutine cal_vof
      implicit none
      
      !call set_surface_velocity(2)
   !   call cal_fgtONE      
      call cal_quvwf
#ifdef NORMAL
      call quvwf_normal
#else
      call quvwf_kaidan
#endif
      call quvwf_modify
      
      call cal_new_f
#ifdef DAKU
!    call turbidity_transport_diffusion
      call turbidity_temperature
#endif
      ! Avoid sliding if possible:
      ! At the moment we can avoid sliding especially in the 2DV case
#ifdef SLIDE1
      ! Only slide nfb = -3, and nfb /= -3 when f next door is < Leps
      call simple_slide(1)
#elif SLIDE0
      ! Slide all nf = 2 or 0
      call simple_slide(0)
#endif
      call set_nflag
end subroutine cal_vof
!        
subroutine cal_quvwf
    use interface_list, only: qff
    use variables, only: inns, inne, ONE
    use arrays, only: in, jn, kn, mn, un, u, qf, dx, lf
    implicit none
    integer :: nnd, nna, iu, i, k, j, nn, l, lg(0:2)
!%%%%%quf‚ÌŒvŽZ%%%%%
!$omp parallel do private(i,j,k,iu,nnd,nna,l,lg)
    do nn = inns,inne
        i = in(nn) ; k = kn(nn) ; j = jn(nn)
        lg = [i, j, k]        
        do l = 0,2
            iu = ( 1 - SIGN(ONE,un(l,nn)) ) / 2
            nnd = mn(i+(iu-1)*lf(l,0),k+(iu-1)*lf(l,2),j+(iu-1)*lf(l,1))
            nna = mn(i-iu*lf(l,0),k-iu*lf(l,2),j-iu*lf(l,1))
            QF(l,nn) = QFF(un(l,nn),dx(l,lg(l)-1+iu),nnd,nna,l+1)
        enddo
    enddo
!$omp end parallel do
end subroutine cal_quvwf
! 
real*8 function QFF(QX,DXX,nnd,nna,nde)
    use arrays, only: in, jn, kn, f, nff
    use variables, only: icon, dt, ZERO, ONE
    implicit none
    real*8,intent(in) :: QX,DXX
    integer,intent(in) :: nnd,nna,nde
    real*8  :: float
    integer :: NFD1, NFD1B, NFA, NFAB

    QFF = ZERO

    if (nnd.eq.0.or.nna.eq.0.or.qx.eq.ZERO) then
        QFF = ZERO
    else
        NFD1 = nff(nnd)%f ; NFD1B = nff(nnd)%b
        NFA  = nff(nna)%f ; NFAB  = nff(nna)%b
        if (NFD1.EQ.1) THEN
            QFF = QX * f(nnd)   
        elseif (NFD1.EQ.-1) then
            if (NFD1B.eq.5) then
                QFF = QX * f(nna)
            elseif (NFD1B.eq.2) then
                QFF = QX * f(nnd)            
            else if(NFD1B.eq.3) THEN
                if (NFA.eq.1) then
                    QFF = QX
                else
                    QFF = QX * f(nnd)
                endif
            elseif(NFD1B.eq.4) THEN 
                if(NFA.eq.1) then       !  C³ŠJŽn for ”r»‚Ì—¬“ü‹«ŠE 100607
                    QFF=QX
                else
                    QFF=QX*f(nnd)
                end if                  !@C³I—¹ for ”r»‚Ì—¬“ü‹«ŠE@100607           
            else !if(NFD1B.eq.10) THEN
                QFF = QX
            endif
        elseif (NFD1.EQ.2) then
            if (ABS(NFD1B).NE.NDE) then            
                    QFF = Qx*f(nnd)
            elseif (ABS(NFD1B).EQ.NDE) then
                if (float(NFD1B)*QX.lt.ZERO) then
                    QFF = SIGN(ONE,QX) * MAX(ABS(QX)-(ONE-f(nnd))*DXX/DT,ZERO)
                else
                    QFF = SIGN(ONE,QX) * MIN(ABS(QX),f(nnd)*DXX/DT)
                endif
            endif
        elseif (NFD1.EQ.0) THEN
            QFF = Qx * f(nnd)          
        else                                   
            ICON = 100
            write(6,*) 'nnd=',nnd,'nna=',nna
        endif
    endif
endfunction QFF  
!        
#ifdef NORMAL
subroutine quvwf_normal
    use interface_list, only: cal_idou2
    use arrays, only: qf, dx, in, jn, kn, mn, a, x, un,sfnoi, nff, lf
    use variables, only: mms, ZERO
    implicit none
    integer :: i, k, j, nn, mm, l, ll, lsn, nnn 
    ! Modified William April 17 2015 - If next cell is nf = 1, the minimum 
    !                                  value of that from cal_quwvf and 
    !                                  here was used. That was inconsistent
    !                                  and has been deleted
    !                                - Rewritten to reduce size to compact 
    !                                  and consistent looping system
!$omp parallel do default(none) &
!$omp shared(mms,sfnoi,in,jn,kn,mn,un,x,dx,qf,nff,lf &
#if defined(DRI) || defined(HENDO)
!$omp          ,a2 &
#else
!$omp          ,a &
#endif
!$omp         ) &
!$omp private(i,k,j,nn,nnn,l,ll,lsn)
      do mm = 1, mms
        nn = sfnoi(mm) ; if (nn.eq.0) cycle
        i = in(nn) ; j = jn(nn) ; k = kn(nn)
        do l = 0,2 ! Loop over the dimensions
            do ll = 0,1 ! Forward and backward
                nnn = mn(i+ll*lf(l,0),k+ll*lf(l,2),j+ll*lf(l,1))
                lsn = 1 - 2*ll
                if (un(l,nnn)*real(lsn).ge.ZERO.or.A_(l,nnn).eq.ZERO) cycle
                qf(l,nnn) = cal_idou2(nff(nn)%b,lsn*(l+1),un(l,nnn),x(0,i),   &
                            x(1,j),x(2,k),dx(0,i),dx(1,j),dx(2,k),nn)/A_(l,nnn)
            enddo
        enddo
      enddo
!$omp end parallel do
end subroutine quvwf_normal
#else
!
subroutine quvwf_kaidan
    use variables, only: isfn, ZERO, ONE
    use arrays, only: sfn, in, jn, kn, un, a, fb, qf, f, mn
    implicit none
    integer :: nnf, nn, i, k, j, nxp, nyp
!$omp parallel do private(nn,i,j,k,nxp,nyp)
    do nnf = 1,isfn        
          nn = sfn(nnf)
          i=in(nn);j=jn(nn);k=kn(nn)
          if(un(0,nn).lt.ZERO.and.a(0,nn).gt.ZERO.and.a(0,nn).lt.ONE.and.a(0,nn).ne.fb(nn)) then
            qf(0,nn)=un(0,nn)*max(ZERO,a(0,nn)-fb(nn)+f(nn)*fb(nn))/a(0,nn)
            continue
          endif
          if(un(1,nn).lt.ZERO.and.a(1,nn).gt.ZERO.and.a(1,nn).lt.ONE.and.a(1,nn).ne.fb(nn)) then
            qf(1,nn)=un(1,nn)*max(ZERO,a(1,nn)-fb(nn)+f(nn)*fb(nn))/a(1,nn)      
          endif
          nxp=mn(i+1,k,j)
          if(un(0,nxp).gt.ZERO.and.a(0,nxp).gt.ZERO.and.a(0,nxp).lt.ONE.and.a(0,nxp).ne.fb(nn)) then
            qf(0,nxp)=un(0,nxp)*max(ZERO,a(0,nxp)-fb(nn)+f(nn)*fb(nn))/a(0,nxp)        
          endif
          nyp=mn(i,k,j+1)
          if(un(1,nyp).gt.ZERO.and.a(1,nyp).gt.ZERO.and.a(1,nyp).lt.ONE.and.a(1,nyp).ne.fb(nn)) then
            qf(1,nyp)=un(1,nyp)*max(ZERO,a(1,nyp)-fb(nn)+f(nn)*fb(nn))/a(1,nyp)        
          endif        
    enddo
!$omp end parallel do
end subroutine
#endif
!
subroutine quvwf_modify
    use interface_list, only: df_vof
    use arrays, only: nff, sfn, mn, qf, lf, kn, jn, in, fb, dx, A_
    use variables, only: isfn, dt, F_lim, ONE
    implicit none
    integer :: i, k, j, nn, l, nnf, nfb0, nn1, lg(0:2)
    real*8  :: dfv2
    !
    ! Modify surface cells that have a F value larger than the F_lim   
!$omp parallel do private(nn,dfv2,i,k,j,nfb0,l,lg,nn1)
    do nnf = 1, isfn
        nn = sfn(nnf)
        dfv2 = df_vof(nn)
        if (dfv2.lt.F_lim) cycle
        i = in(nn) ; k = kn(nn) ; j = jn(nn)
        nfb0 = nff(nn)%b
        l = abs(nfb0)-1; lg(0:2) = [i,j,k]
        if (nfb0.gt.0) then
            nn1 = nn
        else
            nn1 = mn(i+lf(l,0),k+lf(l,2),j+lf(l,1))
        endif
!$omp atomic update 
        qf(l,nn1) = qf(l,nn1) - dx(l,lg(l)) * (dfv2 - ONE) / dt               &
                  * FB_(nn) / A_(l,nn1) * real(sign(1,nfb0))
    enddo
!$omp end parallel do  
end subroutine quvwf_modify
!    
subroutine cal_new_f
      use interface_list,only: df_vof
      use arrays, only: F,in,jn,kn,fp,fb,nff,qf,un
      use variables,only: inns, inne, dt, t, n, icon, ZERO, TINY, F_lim, DT_lim, Leps
      implicit none
      real*8  :: dfv2
      integer :: i,k,j,nn

222 continue

!$omp parallel do default(none) &
!$omp shared(inns,inne,in,jn,kn,nff,fb,dt,f,n,fp,icon) &
!$omp private(dfv2)
      do nn = inns,inne
         if (nff(nn)%f.eq.-1) cycle       
        dfv2 = df_vof(nn)
        if (-dfv2*dt*fb(nn).gt.0.4d0) then
          ! The change in F is too large
          write(6,*) '-dfv2*dt=',-dfv2*dt,'dt=',dt
          write(6,*) 'nn=',nn,in(nn),jn(nn),kn(nn),'n=',n,fp(nn)
          icon = 30
#ifdef BUGWIN2
          dfv2 = df_vof(nn)
#endif
        elseif (abs(dfv2*dt).lt.TINY.and.dfv2.ne.ZERO) then
          ! The new F is small enough to ignore
          dfv2 = ZERO
        endif
        if (dfv2.gt.F_lim) then
          ! F value exceeds the upper limit
          if (nff(nn)%f.eq.2) then
            write(6,*) nff(nn)%f,nff(nn)%b,f(nn),fp(nn),nff(nn)%fp,nff(nn)%bp
            write(6,*)'nn,i,k,j=',nn,in(nn),jn(nn),kn(nn),"error in cal_new_f n=",n
            dfv2 = df_vof(nn)
            write(6,*) 'dfv2=',dfv2,f(nn),fp(nn)
            stop
            continue
          endif
        endif
        f(nn) = dfv2
      enddo
!$omp end parallel do
  221 if (icon.eq.30) then
        t = t - dt
        dt = dt / 1.1d0
        t = t+dt
        write(6,*) 'GENSOKU dt=',dt * 1.1d0,'to dt=',dt,'T=',t,'n=',n
        if (dt.lt.DT_lim) then
          write(6,*) 'nn=',nn,'dt',dt,'stop'
          stop
        endif
        icon = 0
        goto 222
      endif 
endsubroutine cal_new_f
!    
function df_vof(nn)
    use interface_list, only: renzoku_fv
    use arrays, only: qf, mn, in, jn, kn, fp, fb, dx, A_
    use variables, only: dt
    implicit none
    integer,intent(in) :: nn
    real*8             :: df_vof
    integer            :: i,k,j
      
    i = in(nn) ; k = kn(nn) ; j = jn(nn)
#if defined(DRI) || defined(HENDO)
    if (fbd(nn).ne.ZERO) then
        df_vof = fp(nn) * fb(nn) / fbd(nn) - dt * &
                 renzoku_fv( A_(0,mn(I+1,k,j)),A_(0,nn),A_(1,mn(I,k,j+1)),    &
                             A_(1,nn),A_(2,mn(I,k+1,j)),A_(2,nn),             &
                             qf(0,mn(I+1,k,j)),qf(0,nn),qf(1,mn(I,k,j+1)),    &
                             qf(1,nn),qf(2,mn(I,k+1,j)),qf(2,nn),             &
                             dx(0,i),dx(1,j),dx(2,k),fbd(nn) ) 
    else
#endif
        df_vof = fp(nn) - dt *                                                &
                 renzoku_fv( A_(0,mn(I+1,k,j)),A_(0,nn),A_(1,mn(I,k,j+1)),    &
                             A_(1,nn),A_(2,mn(I,k+1,j)),A_(2,nn),             &
                             qf(0,mn(I+1,k,j)),qf(0,nn),qf(1,mn(I,k,j+1)),    &
                             qf(1,nn),qf(2,mn(I,k+1,j)),qf(2,nn),             &
                             dx(0,i),dx(1,j),dx(2,k),fb(nn) ) 
#if defined(DRI) || defined(HENDO)
    endif                   
#endif
end function df_vof
!    
subroutine simple_slide(flag)
    use variables,only: inn2d, is, ie, js, je, ks, ke, ZERO, Leps
    use arrays,only: qf, f, in2d, mn, dx, a, nff,fb
    implicit none
    integer,intent(in) :: flag
    integer :: i, k, j, nn, nnm, nn2d
      
!!!!$omp parallel do private(i,k,j,nn,nnm) 
    do nn2d = 1,inn2d
        i = in2d(0,nn2d) ; j = in2d(1,nn2d)

        do k = ke-1,ks,-1
            nn = mn(i,k,j) ; nnm = mn(i,k-1,j)
            if (nff(nn)%f.eq.-1) exit
            if (nff(nn)%f.eq.1) cycle
            if (nff(nnm)%f.eq.-1) cycle
!#ifdef NORMAL            
    !        if(qf(2,nn).ge.ZERO) cycle
!#endif
            if (flag.eq.1) then
            if ((is.eq.ie-1.or.(i.gt.is.and.i.lt.ie-1)).and.                  &
                (js.eq.je-1.or.(j.gt.js.and.j.lt.je-1))) then
                !(Always slide for cell next to the boundary)
                if (nff(nn)%f.eq.2) then
                    if (nff(nn)%b.ne.-3) then
                        if (f(mn(i-1,k,j)).gt.Leps) cycle
                        if (f(mn(i+1,k,j)).gt.Leps) cycle
                        if (f(mn(i,k,j-1)).gt.Leps) cycle
                        if (f(mn(i,k,j+1)).gt.Leps) cycle
                    endif
                elseif (nff(nn)%f.eq.0) then
#ifdef SATO                    
                        if(fb(nnm).lt.0.1d0)   goto 294
#endif
                    !if(dx(0,i)/dx(2,k).lt.2.0d0) then
                        if ( & !!a(2,mn(i-1,k,j)).eq.ZERO.or.                      &
                            f(mn(i-1,k,j)).gt.Leps) cycle
                        if ( & !!a(2,mn(i+1,k,j)).eq.ZERO.or.                      &
                            f(mn(i+1,k,j)).gt.Leps) cycle
                    !endif
                    !if(dx(1,j)/dx(2,k).lt.2.0d0) then
                        if ( & !a(2,mn(i,k,j-1)).eq.ZERO.or.                      &
                            f(mn(i,k,j-1)).gt.Leps) cycle
                        if ( & !a(2,mn(i,k,j+1)).eq.ZERO.or.                      &
                            f(mn(i,k,j+1)).gt.Leps) cycle
                    !endif
                endif
            endif
            endif
294         call cal_fslide_m0(nn,nnm,-3)
        enddo
    enddo
!!!$omp end parallel do 
end subroutine simple_slide
    
subroutine cal_fslide_m0(nn,nnm,nflag)
    use interface_list, only: df_vof
    use arrays, only: qf, f, in, jn, kn, fb, dx, a, cn, cmax, cmin, nff
    use variables, only: n, dt, t, ONE, ZERO, TINY
    implicit none
    integer,intent(in) :: nn, nnm, nflag
    integer :: l, nnc, lg(0:2), lgm(0:2)
    real*8 :: dvm, f2, vmesh, ftemp, ff1, ff2, dvc, fmp 
    !  
    fmp = f(nnm)
    lg  = [ in(nn),  jn(nn),  kn(nn)  ]
    lgm = [ in(nnm), jn(nnm), kn(nnm) ]
    vmesh = dx(0,lg(0))*dx(1,lg(1))*dx(2,lg(2))
#if defined(DRI) || defined(HENDO) 
    !dvc:‘—‚è‘¤‚Ì—¬‘Ì‘ÌÏAdvm:Žó‚¯‘¤‚Ì—¬‘Ì‘ÌÏ
    dvc = vmesh * fbd(nn)
    dvm = dx(0,lgm(0)) * dx(1,lgm(1)) * dx(2,lgm(2)) * fbd(nnm)  
#else
    dvc = vmesh * fb(nn)
    dvm = dx(0,lgm(0)) * dx(1,lgm(1))  *dx(2,lgm(2)) * fb(nnm)
#endif
    if (dvc*dvm.eq.ZERO) return
    f(nn)  = df_vof(nn)
    f(nnm) = df_vof(nnm)
    !
    !f2:ˆÚ“®‰Â”\—Ê@Žó‚¯‘¤‚Ì‹ó‚«‘ÌÏA‘—‚è‘¤‚Ì‘S‘ÌÏ
    f2 = min( dvm * ( ONE - F(nnm) ) , dvc * f(nn) )  
    if (abs(f2).eq.ZERO) return
    if (dvc.eq.ZERO.or.dvm.eq.ZERO) return
    ff1 = f(nn) ; ff2 = f(nnm)       
    l = abs(nflag) - 1
    if (nflag.lt.0) nnc = nn
    if (nflag.gt.0) nnc = nnm
    qf(l,nnc) = qf(l,nnc) + sign(ONE,real(nflag)) * f2 / vmesh                &
                                                 / dt * dx(l,lg(l)) / A_(l,nnc)
    ftemp = df_vof(nn)
    if (abs(ftemp).lt.TINY) then
        f(nn) = ZERO
    else
        f(nn) = ftemp
    endif
    f(nnm) = df_vof(nnm)
#ifdef DAKU
    if (f2.gt.ZERO) then 
        if(ff2.gt.ZERO) then
            cn(0:1,nnm)=(cn(0:1,nnm)*dvm*ff2+cn(0:1,nn)*f2)/(dvm*f(nnm))
        else
            cn(0:1,nnm)=cn(0:1,nn)                
        endif
        cn(0,nnm)=min(max(cmin(0),cn(0,nnm)),cmax(0))
        cn(1,nnm)=min(max(cmin(1),cn(1,nnm)),cmax(1))
    else
        cn(0:1,nn)=(cn(0:1,nn)*(-f2)+cn(0:1,nn)*dvc*ff1)/(dvc*f(nn))
        cn(0,nn)=min(max(cmin(0),cn(0,nn)),cmax(0))
        cn(1,nn)=min(max(cmin(1),cn(1,nn)),cmax(1))
    endif
#endif
end subroutine cal_fslide_m0
!
subroutine cal_fgtONE
    use variables, only: inns, inne
    use arrays, only:f,nff,in,jn,kn,mn,un
    implicit none
    integer::nn,icc,i,j,k,mn1(4),ii
    doubleprecision::fcc,ff,dff

do nn = inns,inne
        if(nff(nn)%f.ne.0) cycle   
        if(f(nn).lt.0.5d0) cycle   
        i = in(nn) ; k = kn(nn) ; j = jn(nn)
        if(nff(mn(i+1,k,j))%f.eq.0.and.f(mn(i+1,k,j)).eq.0.0d0) then
            un(0,mn(i+1,k,j))=( un(0,nn)+ un(0,mn(i+1,k,j)))*0.5d0
        endif
        if(nff(mn(i-1,k,j))%f.eq.0.and.f(mn(i-1,k,j)).eq.0.0d0) then
            un(0,nn)=( un(0,nn)+ un(0,mn(i+1,k,j)))*0.5d0
        endif
        if(nff(mn(i,k,j+1))%f.eq.0.and.f(mn(i,k,j+1)).eq.0.0d0) then
            un(1,mn(i,k,j+1))=( un(1,nn)+ un(1,mn(i,k,j+1)))*0.5d0
        endif
        if(nff(mn(i,k,j-1))%f.eq.0.and.f(mn(i,k,j-1)).eq.0.0d0) then
            un(1,nn)=( un(1,nn)+ un(1,mn(i,k,j+1)))*0.5d0
        endif
        if(nff(mn(i,k-1,j))%f.eq.0.and.f(mn(i,k-1,j)).eq.0.0d0) then
            un(2,nn)=-0.01d0
        endif
    enddo
end subroutine