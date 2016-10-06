!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-momentum.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Updates the velocity using the SMAC algorithm
!    : Explicitly calculate advection and viscous stresses
!    : Estimate new velocity via momentum equation
!    : Correct pressure in PPE to obtain divergence free velocity field
!    : Obtain new correct velocities from corrected pressure
#if defined (DRI) || defined (HENDO)
#define A_ a2
#define FB_ fbd
#else
#define A_ a
#define FB_ fb
#endif
#ifdef NORMAL
#define FS_ fss      
#else
#define FS_ f(nn)
#endif

 subroutine cal_momentum
    use arrays, only: u, un, ox, p, xe, mno2mn
    use variables, only: TOME, mnoe, ik, ricon
    implicit none
    integer :: nnf, nn, loop
  
    call set_pressure_on_surface
    
    IMPROVED_EULER_ITER: do loop = 1,TOME
      ! Start loop of the momentum calculation with TOME (1st or 2nd) temporal order
      if (loop.eq.1) then
          call cal_advection_viscosity(u,ox(loop,:,:))
      elseif (loop.eq.2) then
          call set_surface_velocity(2)
          call cal_advection_viscosity(un,ox(loop,:,:))
      endif      
      !Update velocity
      call estimate_velocity(loop)

    enddo IMPROVED_EULER_ITER   

    ik = 0
    do !Start loop of Pressure correction to satisfy the divergence free condition
        call cal_err_of_continuity        
        call check_dmax

        if (ricon.ne.0) return
        
        call cal_matrix_element
        
        call cal_poisson_equation
        
        !Update the pressure
!$omp parallel do private(nn) 
         do nnf = 1,mnoe
            nn = mno2mn(nnf)
            p(nn) = p(nn) - xe(nnf)
         enddo
!$omp end parallel do
        call estimate_velocity2       
    enddo
      
end subroutine cal_momentum

subroutine set_pressure_on_surface
!C*********************************
!C  êÖñ Ç≈ÇÃà≥óÕÇÃã´äEê›íË
!C*********************************
    use interface_list, only: pshu
    use arrays, only: nor, sfn, sp, u, p, x, dx, in, jn, kn, f, fb, mn,  &
                      nff, d, rho, lf
    use variables,only: isfn, nu, g, inns, inne, iob, ZERO, ONE, HALF, TWO
    implicit none
    real*8 :: PS, fss, ro = ONE
    integer,dimension(0:2) :: lg
    integer :: isx, i, k, j, nn, l, nfm, nnm, loop, mm
    
!!$omp parallel do
!do nn=inns,inne
!    if(nff(nn)%f.ge.1) p(nn)=ZERO
!enddo
!!$omp end parallel do 
ITER_LOOP: do loop = 1,2
!$omp parallel do default(none) shared(d,dx,g,in,jn,kn,mn,nff,            &
!$omp nu,p,u,x,rho,sfn,sp,nor,isfn,ro,loop,lf,f,fb)                       &
!$omp private(nn,i,k,j,ps,isx,l,lg,nfm,nnm,fss) schedule(dynamic,250)
    do mm = 1,isfn
!C=================================
      nn = sfn(mm)
      if (nff(nn)%f.ne.2) cycle
      i = in(nn); j = jn(nn); k = kn(nn)
      lg = [i,j,k]
#ifdef DAKU
      ro = rho(nn)
#endif
      isx = sign(1,nff(nn)%b)
      l = abs(nff(nn)%b) - 1
      nnm = mn(i+isx*lf(l,0),k+isx*lf(l,2),j+isx*lf(l,1))
      nfm = nff(nnm)%f
      
      if (loop.eq.1) then
          ! First loop interpolating from nf =/ 2 only
          if (nfm.eq.2) cycle
      elseif (loop.eq.2) then
          ! Second loop interpolating from nf = 2 only
          if (nfm.ne.2) cycle
      endif
#ifdef NORMAL
      ! Get the normalised distance from the surface to the cell edge
      if (isx.eq.1) then
        fss = (x(l,lg(l)+1)-sp(nn,l))/dx(l,lg(l))
      else
        fss = (sp(nn,l)-x(l,lg(l)))/dx(l,lg(l))
      endif
#else
      fss = ( ONE - fb(nn) ) + f(nn) * fb(nn)
#endif
#ifdef RAN
      !Calculating the viscous stress on the surface 
      ps = ( nor(nn,0) * (u(0,mn(i+1,k,j))-u(0,nn))/dx(0,i)   &
           + nor(nn,1) * (u(1,mn(i,k,j+1))-u(1,nn))/dx(1,j)   &
           + nor(nn,2) * (u(2,mn(i,k+1,j))-u(2,nn))/dx(2,k) ) &
           * (d(nn) + nu ) * TWO
#else 
      ps = ZERO
#endif
      if (nfm.eq.iob) then
        ! Assume hydrostatic pressure
        if (l.eq.2) then
          p(nn) = ps + (fss - HALF) * dx(l,lg(l)) * g * ro
        else
          p(nn) = ZERO
        endif
      else 
        p(nn) = PSHU(p(nnm),ps,dx(l,lg(l)),dx(l,lg(l)+isx),fss)  
      endif
!C=================================
    enddo
!$omp end parallel do
enddo ITER_LOOP 

end subroutine set_pressure_on_surface       
    
subroutine estimate_velocity(loop)
    use interface_list, only : momentum, momentum_rho
    use variables, only : iob, inns, inne, g, dt, TWOTHIRD, HALF, ZERO, TINY, &
                          TOME
    use arrays, only : in, jn, kn, nff, p,mn ,a, dx, un, u, ox, rho,     &
                       rhon, fc, r
    implicit none
    integer,intent(in) :: loop
    real*8 :: pc, dxc, ff, pm, oxc
    integer :: i, k, j, nn, l
    integer,dimension(0:2) :: ic, nnp, nnm, nfm, nfbm

!$omp parallel do private(i,j,k,ic,pc,pm,nnp,nnm,nfm,nfbm,l,dxc,ff,oxc) !,rc,rcn,nnm2,rcm,vcc1
    do nn = inns,inne
      if (nff(nn)%f.eq.0) cycle 
      if (nff(nn)%f.eq.iob.and.nff(nn)%b.ne.4) cycle
        i = in(nn) ; j = jn(nn) ; k = kn(nn)    
        ic(0:2) = [i,j,k]
#ifdef RAN !Add on the contribution from the turbulent energy
        pc = p(nn) + TWOTHIRD * r(nn)%k 
#else
        pc = p(nn) 
#endif
        nnp(0:2) = [mn(i+1,k,j),mn(i,k,j+1),mn(i,k+1,j)]
        nnm(0:2) = [mn(i-1,k,j),mn(i,k,j-1),mn(i,k-1,j)]
        nfm(0:2) = [nff(mn(i-1,k,j))%f,nff(mn(i,k,j-1))%f,nff(mn(i,k-1,j))%f]
        nfbm(0:2)= [nff(mn(i-1,k,j))%b,nff(mn(i,k,j-1))%b,nff(mn(i,k-1,j))%b]
        
        if (nff(nn)%f.eq.iob.and.nff(nn)%b.eq.4) then
            ! For special case of pressure boundary cell
            do l = 0,2
                if (nfm(l).ge.1) then
                    dxc = (DX(l,ic(l))+DX(l,ic(l)-1)) * HALF
#ifdef RAN          !Add on the contribution from the turbulent energy
                    pm = p(nnm(l)) + TWOTHIRD * r(nnm(l))%k 
#else
                    pm = p(nnm(l))
#endif            
#ifndef DAKU
                    ff = MAX(u(l,nn),ZERO)*(u(l,nn)-u(l,nnm(l)))/dx(l,ic(l)-1)+fc(l)
                    un(l,nn)  =momentum(u(l,nn),ZERO,pm,pc,dxc,ff,dt)
#else
                    RCN=(RHON(nn) + RHON(nnm(l)) ) * HALF
                    RC= (RHO(nn) + RHO(nnm(l)) ) * HALF
                    RCm= (RHO(nnm(l)) + RHO(nnm2(l)) ) * HALF
                    vcc1=u(l,nn)*0.9d0+u(l,nnm(l))*0.1d0
                    ff  =MAX(vcc1,ZERO)*(rc*u(l,nn)-rcm*u(l,nnm(l)))/dx(l,ic(l)-1)+fc(l)*RC
                    un(l,nn) = momentum_rho(u(l,nn),rc,rcn,ZERO,Pm,Pc,dxc,ff,dt)
#endif
                endif
            enddo               
            
        else ! For normal case of fluid or surface cell
            do l = 0,2
                if (A_(l,nn).ne.ZERO) then
#ifdef RAN          !Add on the contribution from the turbulent energy
                    pm = p(nnm(l)) + TWOTHIRD * r(nnm(l))%k 
#else
                    pm = p(nnm(l))
#endif
                    dxc = (DX(l,ic(l))+DX(l,ic(l)-1)) * HALF
                    if (loop.eq.1) then
                        oxc = ox(1,l,nn)
                    elseif (loop.eq.2) then
                        oxc = HALF * sum(ox(1:2,l,nn))
                    endif
                    if (nfm(l).ge.1) then
#ifndef DAKU
                        un(l,nn) = momentum(u(l,nn),oxc,pm,pc,dxc,fc(l),dt)
                        if (abs(un(L,NN)) < TINY) UN(L,NN) = ZERO  ! nagashima (2015.02.24)
#else
                        RC= ( RHO(nn) + RHO(nnm(l)) ) * HALF
                        RCN= ( RHON(nn) + RHON(nnm(l)) ) * HALF
                        un(l,nn) = momentum_rho(u(l,nn),rc,rcn,oxc,Pm,pc,dxc,fc(l)*rc,dt)
#endif
                    elseif (nfm(l).eq.-1.and.nfbm(l).eq.4) then 
                        ff = (u(l,nnp(l))-u(l,nn)) / dx(l,ic(l)) * min(u(l,nn),ZERO) + fc(l)
#ifndef DAKU
                        un(l,nn) = momentum(u(l,nn),ZERO,pm,pc,dxc,ff,dt)
#else       
                        RC = RHO(nn)
                        un(l,nn) = momentum_rho(u(l,nn),rho(nn),rhon(nn),ZERO,Pm,pc,dxc,ff*rc,dt)
#endif
                    endif
                endif
            enddo
        endif
    enddo
!$omp end parallel do
    
end subroutine estimate_velocity    

subroutine cal_err_of_continuity
    use interface_list, only: renzoku_n, renzoku_p
    use variables, only: ifln, is, ie, js, je, ks, ke, iob, dt, n, m_alloc,   &
                          ik, nl, mnoe, dmax, ZERO, TINY, ONE, mp, t
    use arrays, only: fln, in, jn, kn, fb, nff, f, p, un, mn, dx, ox, u, b, xe&
                      , de, mno, mno2mn, alloc_matrix, dealloc_matrix, lf
    implicit none
    integer :: nnn, l, ll, nnd1, dr1_large
    integer :: nnd, nnf, i, j, k, nn
    real*8 :: dr1,dmax1,fbc,ff 

#ifdef BUGWIN
    if (ifln.eq.0) then
      write(6,*) 'ifln=',ifln
      continue
    endif
#endif
    ! Allocate and resize the matrices if required
    if (m_alloc.eq.0.and.ik.eq.0) then
      m_alloc = ifln + 5000
      call alloc_matrix(m_alloc,ie,je,ke)
    elseif (ifln+60.gt.m_alloc) then
      write(6,*) 'mnoe,m_alloc1=',mnoe,m_alloc,'in_calerr3'
      call dealloc_matrix
      m_alloc = ifln + 5000
      write(6,*) 'mnoe,m_alloc2=',mnoe,m_alloc,'in_calerr3'
      call alloc_matrix(m_alloc,ie,je,ke)
      write(6,*) 'mnoe,m_alloc3=',mnoe,m_alloc,'in_calerr3'
    endif
    
111 nnd = 0       
    dmax = ZERO ; dmax1 = ZERO
!$omp parallel do
    do nnf = 1,ubound(b,1)
        b(nnf)  = ZERO ; xe(nnf) = ZERO ; de(nnf) = ZERO
    enddo
!$omp end parallel do
    
    if (ik.eq.0) then
      nn = 0
!$omp parallel workshare
      mno = 0
!$omp end parallel workshare
      do nnf = 1,ifln
        nnn = fln(nnf)
        i = in(nnn) ; j = jn(nnn) ; k = kn(nnn)
#if defined(DRI) || defined(HENDO)
        fbc = fb(nnn) + fbd(nnn)
#else
        fbc = fb(nnn)
#endif
        if (nff(nnn)%f.eq.0.or.nff(nnn)%f.eq.iob.or.fbc.eq.ZERO) cycle
        if (nff(nnn)%f.eq.2) then
            l = abs(nff(nnn)%b) - 1; ll = sign(1,nff(nnn)%b)
            if (nff(mn(i+ll*lf(l,0),k+ll*lf(l,2),j+ll*lf(l,1)))%f.eq.iob) cycle
        endif    
        nn = nn + 1
        mno(i,k,j) = nn
        mno2mn(nn) = nnn
      enddo

      mnoe = nn       
      if (nn.lt.10) then
        write(6,*) 'ik=',ik,mnoe,' in_calerr1'
        return
      endif
      ! Check if the size of the matrix is reasonable
      if (mnoe+60.gt.m_alloc) then
        write(6,*) 'mnoe,m_alloc1=',mnoe,m_alloc,'in_calerr3'
        call dealloc_matrix
        m_alloc=mnoe+60+10000
        write(6,*) 'mnoe,m_alloc2=',mnoe,m_alloc,'in_calerr3'
        call alloc_matrix(m_alloc,ie,je,ke)
        write(6,*) 'mnoe,m_alloc3=',mnoe,m_alloc,'in_calerr3'
        goto 111
      endif   
    endif

    dr1_large = 0    
!$omp parallel do private(ff,nn,dr1) reduction(max:dmax) schedule(dynamic,250)
    do nnf = 1,mnoe
      nn = mno2mn(nnf)
      if (nff(nn)%f.ne.1) cycle
   !  
#if defined(DRI) || defined(HENDO)
        DR1 = cal_drr(nn)
#else

#ifndef NORMAL
        if (f(nn).lt.0.3d0) then
          ff = max(-0.1,min(1.0d0-f(nn),0.1d0))
        elseif (f(nn).lt.0.6d0) then
          ff = max(-0.01,min(1.0d0-f(nn),0.01d0))
        else
          ff = max(-0.001,min(1.0d0-f(nn),0.001d0))
        endif
        DR1 = renzoku_n(nn) + ff / dt
#else
        ff  = max(-0.001d0,min(ONE - f(nn),0.001d0))
        DR1 = renzoku_n(nn) + ff / dt
#endif

#endif
        if (dr1_large.eq.0) then
            if (abs(dr1).gt.5000d0.or.(ik.gt.1.and.abs(dr1).gt.500d0)) then
                dr1_large = 1
            endif
        endif
        if (abs(DR1/DT).le.TINY.and.abs(DR1/DT).ne.ZERO) then
        else
          b(nnf) = mp * DR1 / DT
        endif
#ifndef BUGWIN
        if (ABS(DR1).GT.DMAX) then
          nnD = nn
          DMAX = ABS(DR1)
        endif
        if (ABS(DR1).GT.DMAX1.and.k.lt.5) then
          nnD1  = nn
          DMAX1 = ABS(DR1)
        endif
#else
        DMAX = max(Dmax,abs(dr1))
#endif
    enddo
!$omp end parallel do
    if (dr1_large.eq.1) then
        write(6,*) 'dr1 is large=', dr1, 'ik loop=', ik
    endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#ifdef BUGWIN
    call showdata
#endif

contains


subroutine showdata
    use variables, only: d0, t
    use arrays, only: u
    implicit none
    integer :: id,kd,jd,nnc
    
    IF((IK.GE.6.or.dmax.gt.5000d0).and.nnd.ne.0)  then
!cc      nnd=2
      nnc=nnd
      id=in(nnc)
      jd=jn(nnc)
      kd=kn(nnc)

      write(6,*) 'ik=',ik,abs(b(mno(id,kd,jd)))*dt,dmax,nnd,id,kd,jd
      if(abs(b(mno(id,kd,jd)))*dt.gt.d0) then
        WRITE(6,341) IK,b(mno(id,kd,jd))*dt,ID,KD,JD,nff(nnc)%f,js,je,ke
        write(6,*) 'n=',n,'t=',t
        write(6,*) un(0,nnc),un(0,mn(id+1,kd,jd)),'xn'
        write(6,*) u(0,nnc),u(0,mn(id+1,kd,jd)),'x'
        write(6,*) un(1,nnc),un(1,mn(id,kd,jd+1)),'yn'
        write(6,*) u(1,nnc),u(1,mn(id,kd,jd+1)),'y'
        write(6,*) un(2,nnc),un(2,mn(id,kd+1,jd)),'wn'
        write(6,*) u(2,nnc),u(2,mn(id,kd+1,jd)),'w'
        write(6,*) p(nnc),p(mn(id,kd+1,jd)),'p,pzp'
!        write(6,*) p(mn(id+1,kd,jd)),'pxp'
        write(6,*) nff(nnc)%f,nff(mn(id,kd+1,jd))%f,'nf'
        write(6,*) f(nnc),f(mn(id,kd+1,jd)),'fz'
        write(6,*) (UN(0,mn(Id+1,kd,jd))-UN(0,nnc))/DX(0,Id) &
                   +(un(1,mn(Id,kd,jd+1))-un(1,nnc))/Dx(1,Jd) &
                   +(un(2,mn(Id,kd+1,jd))-un(2,nnc))/Dx(2,Kd) &
                  ,(1.0d0-f(nnc))*0.1d0/dt*0.01d0,dt
#if defined(DRI) || defined(HENDO)
        write(6,*) 'fb,fbd=',fb(nnc),fbd(nnc)
#endif
#ifdef BUGWIN
        continue
#endif
      end if
    end if
341 format(1h ,i4,'d=',e12.5,'ikj=',3(1x,i4),'nf=',i3 &
           ,'js=',i3,' je=',i4,' ke=',i3)

end subroutine showdata           

end subroutine cal_err_of_continuity 

subroutine check_dmax
    use variables, only: icon, d0, t, n, ricon, mnoe, ik, dmax, iter, eps
    implicit none
    real*8 :: dmax1
    
    ricon = 0
    if (icon.ne.0) then
      write(6,*) 'mnoe=',mnoe,' in_shu5',icon
      icon = 1030 ; ricon = 1
      return
    end if
    
    if (ik.eq.0.and.dmax.gt.5000d0) then  
      dmax1 = dmax
    elseif (ik.ge.900.and.dmax.gt.d0) then
      write(6,*) 'dmax=',dmax,'t=',t,'n=',n,'ik=',ik,'iter=',iter
    endif

    if (DMAX.LT.D0) THEN
      if (icon.eq.3000) icon = 0
      if(ik.gt.100)    write(6,*) 't,dmax,iter,ik,eps=',t,dmax,iter,ik,eps
      ricon = 1
      return
    else
      if(icon.eq.12) then
        icon=0
        ricon=1
        return
      endif
      IK = IK + 1
      if (ik.ge.900) then
        write(6,*) 't,dmax,iter,ik,eps=',t,dmax,iter,ik,eps
        if (ik.gt.1000) then
          ricon = 1
          return
        endif
    !  elseif (ik.eq.100) then
    !    write(6,*) 'ik=',ik
    !    icon  = 100
    !    ricon = 1
    !    return
      endif
    endif
    
endsubroutine check_dmax    
    
subroutine cal_matrix_element
    use variables, only : ifln, ifu, iob, icon, nl, is, ie, js, je, ZERO, ONE
    use arrays, only : fln, in, jn, kn, a, sp, x, mn, dx, fb, rhon, rho, nff, &
                       alu, mno, de, llu, DE_C, alu_c, ll_c
    use interface_list, only: cal_alu
    implicit none
    real*8  :: fss
    integer :: nn, i, k, j, nnn, nnf, l, np, ii, ll
    real*8,parameter :: rate = -1.0d-20
! Update on July 30th by William Pringle
! -> Changed so that initial coefficients are saved in "init_matrix" subroutine
! -> Loop through size 1 to nl (4 or 6) and insert saved coefficients if NF/= 0
! -> For surface cells we still need to calculate alu based on the current F
!    value or sp value in the direction of the surface normal

!$omp parallel do private(i,j,k,l,nn,nnn,np,ll,ii,fss) !,fb2
Fluid_cell: DO nnf = 1,ifln
        nnn = fln(nnf)
        i = in(nnn); j = jn(nnn); k = kn(nnn)
        nn  = mno(i,k,j) 
        if (nn.eq.0) cycle
        if (nff(nnn)%f.eq.ifu) then
            de(nn) = DE_C(nnn)
        elseif (nff(nnn)%f.eq.2) then
            de(nn) = ONE
        else
            write(6,*) 'Error in Cal_mat_ele: Unexpected nflag =',nff(nnn)%f
            stop
        endif
        ! Loop over the size of the matrix elements
        do ll = 1,nl
            if (LL_C(nnn,ll).eq.0) then
                ! Surrounding cell is a boundary cell (nfb /= 4)
                llu(nn,ll) = 0
                alu(nn,ll) = ZERO
            else
                NP = mno(in(LL_C(nnn,ll)),kn(LL_C(nnn,ll)),jn(LL_C(nnn,ll)))
                if (nff(LL_C(nnn,ll))%f.eq.0.or.NP.eq.0) then
                    ! Surrounding cell is an air cell
                    llu(nn,ll) = 0
                    alu(nn,ll) = ZERO 
                    ! We need to subtract from the stored amount 
                    if (nff(nnn)%f.eq.ifu.and.nff(LL_C(nnn,ll))%f.ne.2) then
                        de(nn) = de(nn) + ALU_C(nnn,ll)
                    endif     
                else! Surrounding cell is fluid or surface cell
                    ! Set the stored mapping element
                    llu(nn,ll) = NP
                    if (nff(nnn)%f.eq.ifu) then
                        ! Fluid Cell: Can use the stored values as is
                        alu(nn,ll) = ALU_C(nnn,ll)
                    else 
                        ! Surface Cell: need to check status of orientation
                        l = (ll + mod(ll,2))/2 - 1
                        if (is.eq.ie-1) l = l + 1
                        if (js.eq.je-1.and.ll > 2) l = l + 1
                        if (abs(nff(nnn)%b)-1.eq.l) then
                            ! Fluid or surface cell in normal direction
                            if (l.eq.0) ii = i
                            if (l.eq.1) ii = j
                            if (l.eq.2) ii = k
#ifdef NORMAL
                            fss = real(1-2*mod(ll,2)) * ( x(l,ii-mod(ll,2)+1) &
                                                       - sp(nnn,l) ) / dx(l,ii)
                            alu(nn,ll) = cal_alu(fss,dx(l,ii),                &
                                                 dx(l,ii-2*mod(ll,2)+1))
#else
#ifndef DAKU
                            alu(nn,ll) = cal_alu(f(nnn),dx(l,ii),             &
                                                 dx(l,ii-2*mod(ll,2)+1))
#else
                            alu(nn,ll) = cal_alu_rho(rho(nnn),rhon(np),f(nnn),&
                                               dx(l,ii),dx(l,ii-2*mod(ll,2)+1))
#endif
#endif
                        else
                            ! Fluid or surface cell not in normal direction
                            alu(nn,ll) = rate  
                        endif
                    endif
                endif
            endif
        enddo    
    enddo Fluid_cell
!$omp end parallel do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end subroutine cal_matrix_element   
    
subroutine cal_poisson_equation
    use variables, only : eps0, icon, eps, iter, nl, mm2, mnoe, ZERO
    use arrays, only : alu, llu, de, b, xe
    implicit none
    real*8,dimension(2) :: eps1
    integer :: ier = 0, n2
    real*8 :: parm = 5.d0

    iter = 10000
    eps1(1) = eps0; eps1(2) = ZERO
    n2 = mm2 - mnoe
         call PCGPME(alu,nl,mm2,llu,de,mnoe,64,b,xe(1:mm2),eps1,parm,iter,ier)
    !call PCGPME_Para(alu,nl,mm2,llu,de,mnoe,64,b,xe(1:mm2),eps1,parm,iter,ier) 
    eps = eps1(1)
    
    if (ier.eq.3000) then
      write(6,*) 'ier,mnoe,mm2=',ier,mnoe,mm2 &
                ,'iter=',iter,'eps1=',eps1(1) &
                ,'eps1/eps0=',eps1(1)/eps0
    elseif (ier.eq.2100) then
      write(6,*) 'ier=',ier
      icon=96
    elseif (ier.eq.2200) then
      write(6,*) 'ier=',ier
      icon=96
    elseif (ier.eq.2300) then
      write(6,*) 'ier=',ier,mnoe,mm2
      icon=96
    elseif (ier.eq.12) then
      write(6,*) 'ier=',ier,mnoe,mm2
      icon=12
    elseif (ier.ne.0) then
      write(6,*) 'ier=',ier
      icon=911
    endif
end subroutine cal_poisson_equation

subroutine estimate_velocity2
    use interface_list,only : momentum_d, momentum_rho_d
    use variables, only : inns, inne, iob, dt, ZERO, HALF, TINY
    use arrays, only: in, jn, kn, mn, a, dx, un, rhon, nff, mno, xe
    implicit none
    integer :: i, j, k, nn, nnn, l
    integer :: ic(0:2), nno(0:2), nnm(0:2), nfm(0:2), nfbm(0:2)
    real*8 :: xec,dxc

!$omp parallel do private(i,j,k,ic,nnn,nno,nnm,nfm,nfbm,l,dxc,xec) !,rucn
    do nn = inns,inne
      if (nff(nn)%f.eq.0) cycle 
      if (nff(nn)%f.eq.iob.and.nff(nn)%b.ne.4) cycle
        i = in(nn) ; j = jn(nn) ; k = kn(nn)    
        ic(0:2) = [i,j,k]
        nnn = mno(i,k,j)
        xec = -xe(nnn)
        nno(0:2) = [mno(i-1,k,j),mno(i,k,j-1),mno(i,k-1,j)]
        nnm(0:2) = [mn(i-1,k,j),mn(i,k,j-1),mn(i,k-1,j)]
        nfm(0:2) = [nff(mn(i-1,k,j))%f,nff(mn(i,k,j-1))%f,nff(mn(i,k-1,j))%f]
        nfbm(0:2)= [nff(mn(i-1,k,j))%b,nff(mn(i,k,j-1))%b,nff(mn(i,k-1,j))%b]
        if (nff(nn)%f.eq.iob.and.nff(nn)%b.eq.4) then
            ! For special case of pressure boundary cell
            do l = 0,2
                if (nfm(l).ge.1) then
                    dxc = ( dx(l,ic(l)) + dx(l,ic(l)-1) ) * HALF
#ifndef DAKU
                    un(l,nn) = momentum_d(un(l,nn),-xe(nno(l)),-xe(nnn),dxc,dt)
#else
                    RUCN = (RHON(nn) + RHON(nnm(l)) ) * HALF
                    UN(l,nn) = momentum_rho_d(un(l,nn),rucn,-xe(nno(l)),-xe(nnn),dxc,dt)
#endif
                endif
            enddo               
            
        else ! For normal case of fluid or surface cell
            do l = 0,2
                if (A_(l,nn).ne.ZERO) then
                    dxc = ( dx(l,ic(l)) + dx(l,ic(l)-1) ) * HALF
                    if (nfm(l).ge.1) then
#ifndef DAKU
           !         if(nnn.eq.0.or.nno(l).eq.0) cycle
                    un(l,nn) = momentum_d(un(l,nn),-xe(nno(l)),xec,dxc,dt)
                    if (abs(un(L,NN)) < TINY) UN(L,NN) = ZERO  ! nagashima (2015.02.24)
#else
                        RUCN = ( RHON(nn) + RHON(nnm(l)) ) * HALF
                        un(l,nn) = momentum_rho_d(un(l,nn),rucn,-xe(nno(l)),xec,dxc,dt)
#endif
                    elseif (nfm(l).eq.-1.and.nfbm(l).eq.4) then
#ifndef DAKU
                        un(l,nn) = momentum_d(un(l,nn),-xe(nno(l)),xec,dxc,dt)
#else
                        un(l,nn) = momentum_rho_d(un(l,nn),rhon(nn),-xe(nno(l)),xec,dxc,dt)
#endif
                    endif
                endif
            enddo
        endif
    enddo
!$omp end parallel do

end subroutine estimate_velocity2

#ifdef DAKU2
real*8 function momentum_rho(uc,rc,rcn,ox,pm,pc,dx,f,dt)
    implicit none
    real*8,intent(in)::uc,rc,rcn,ox,pm,pc,dx,f,dt
    momentum_rho=(uc*rc+dt*((Pm-Pc)/dx-f+OX))/rcn
end function


real*8 function momentum_rho2(uc,rc,rcn,ox,ox2,pm,pc,dx,f,dt)
    implicit none
    real*8,intent(in)::uc,rc,rcn,ox,ox2,pm,pc,dx,f,dt
    real*8::oxx
!    OXx=1.5d0*OX-HALF*ox2
    OXx=HALF*OX+HALF*ox2
    momentum_rho2=momentum_rho(uc,rc,rcn,oxx,pm,pc,dx,f,dt)
end function


real*8 function momentum_rho_d(uc,rcn,pm,pc,dx,dt)  !,rc,
    implicit none
    real*8,intent(in)::uc,rcn,pm,pc,dx,dt !,rc
    momentum_rho_d=uc+(dt*(Pm-Pc)/dx)/rcn
end function
#endif   
    
 
