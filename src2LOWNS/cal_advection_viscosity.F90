!%%%%%%%%%%%% FILE: cal_advection_viscosity.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the advection and viscous stress terms in the momentum equations
subroutine cal_advection_viscosity(U1,OX1)
    use type_list, only: velocity
    use arrays, only: a, d, dx, fln, in, jn, kn, mn, rho, nff
    use variables, only: aratez, inns, inne, inne2, ifln, nu, ratex,          &
                         ZERO, ONE, HALF, TWO
    use interface_list, only: iryu_term_x, iryu_term, BIBN2
    implicit none
    real*8,dimension(0:2,0:inne2+1),intent(in)    :: U1
    real*8,dimension(0:2,0:inne2+1),intent(inout) :: OX1
    real*8 :: vuc, wuc
    real*8 :: fux, fuy, fuz, visux, visuy, visuz
    real*8 :: kadox
    type(velocity),dimension(4) :: VX, VY, VZ
    type(velocity) :: VS
#ifdef DRI
    integer::mm,ii,ll,jj
    real*8::xmax,xmin,ymax,ymin,zmax,zmin
#endif
#ifdef DAKU
    type(velocity),dimension(4) :: RVX, RVY, RVZ
    type(velocity) :: RVS
    real*8 :: RX, RY, RZ
#ifdef FURYOKU
    real*8 :: bfp, bfm, akcf, ecf, rhocf
#endif
#endif
    integer :: nnf, m0
    integer :: i, k, j, nn
    integer :: l, l1, m1
    integer :: i1, k1, j1
    integer :: i2, k2, j2
    
!$omp parallel default(none)                                                  &
!$omp shared(ifln,fln,in,jn,kn,a,nff,u1,ox1,dx,mn,d,ratex,aratez,nu,inns,inne)&
!$omp private(nnf,nn,i,k,j,l,m0,i1,k1,j1,l1,m1,i2,k2,j2,vs,vx,vy,vz           & 
!$omp         ,VISUX,visuy,visuz,fux,fuy,vuc,wuc,fuz,kadox                    &
#ifdef DAKU
!$omp         ,RX,RVS,RVX,RVY,RVZ                                             &
#ifdef FURYOKU
!$omp         ,rhocf,akcf,ecf,Bfm,Bfp                                         &
#endif
#endif
!$omp         ) 
!$omp do schedule(dynamic,500)
    do nn = inns,inne
        if (nff(nn)%f.eq.0) OX1(0:2,nn) = ZERO
    enddo
!$omp end do
!$omp do schedule(dynamic,500)
ADV_LOOP: do nnf = 1,ifln
        nn = fln(nnf)
        i = in(nn) ; j = jn(nn) ; k = kn(nn)

l_LOOP: do l = 0,2
          if (l.eq.0) then
            m0 = i ; i1 = 1 ; k1 = 0 ; j1 = 0
          elseif (l.eq.1) then
            m0 = j ; i1 = 0 ; k1 = 0 ; j1 = 1
          elseif (l.eq.2) then
            m0 = k ; i1 = 0 ; k1 = 1 ; j1 = 0
          endif
          if (a(l,nn).eq.ZERO) cycle
          if (nff(mn(i-i1,k-k1,j-j1))%f.lt.1) cycle
          VS%x  = ZERO ; VS%v = u1(l,nn) ; VS%ap = a(l,nn) ;  VS%np = nff(nn)%f
#ifdef DAKU
          RVS%x = ZERO
#endif
          VX(1)%v  = u1(l,mn(i+i1,k+k1,j+j1))
          VX(1)%x  = dx(l,m0)
          VX(1)%np = nff(mn(i+i1,k+k1,j+j1))%f
          VX(1)%ap = a(l,mn(i+i1,k+k1,j+j1))
          
          VX(2)%v  = u1(l,mn(i-i1,k-k1,j-j1))
          VX(2)%x  = -dx(l,m0-1)
          VX(2)%np = nff(mn(i-i1,k-k1,j-j1))%f
          VX(2)%ap = a(l,mn(i-i1,k-k1,j-j1))
          
          VX(3)%v  = u1(l,mn(i+i1*2,k+k1*2,j+j1*2))
          VX(3)%x  = dx(l,m0) + dx(l,m0+1)
          VX(3)%np = nff(mn(i+i1*2,k+k1*2,j+j1*2))%f
          VX(3)%ap = a(l,mn(i+i1*2,k+k1*2,j+j1*2))
          
          VX(4)%v  = u1(l,mn(i-i1*2,k-k1*2,j-j1*2))
          VX(4)%x  = -dx(l,m0-1) - dx(l,m0-2)
          VX(4)%np = nff(mn(i-i1*2,k-k1*2,j-j1*2))%f
          VX(4)%ap = a(l,mn(i-i1*2,k-k1*2,j-j1*2))
          !-------------------------------------------------------------------!  
          ! Added by William 12/31 for improved performance on slopes         !
          ! -> Method can be justified by imagining a change of coordinate    !
          !    systems that runs along the slope would eliminate the problem. !
          !    Thus the diagonally opposite value is equivalent to the        ! 
          !    adjacent value in the current problem                          ! 
          !-------------------------------------------------------------------!
          if (l.ne.2) then
            !For ADV and VISCOUS term in horizontal direction
            if (VX(1)%ap == ZERO .and. nff(mn(i+i1,k+1,j+j1))%f >= 0) then
                VX(1)%v  = u1(l,mn(i+i1,k+1,j+j1))
            endif
            if (VX(2)%ap == ZERO .and. nff(mn(i-i1,k+1,j-j1))%f >= 0) then
                VX(2)%v  = u1(l,mn(i-i1,k+1,j-j1))
            endif
          endif
          
!CCC  CONVECTION TERM  OF U CCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC  DUU/DX  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          if(l.eq.0) then
            l1 = 1 ; m1 = j ; i2 = 0 ; k2 = 0 ; j2 = 1
          else if(l.eq.1) then
            l1 = 2 ; m1 = k ; i2 = 0 ; k2 = 1 ; j2 = 0
          else if(l.eq.2) then
            l1 = 0 ; m1 = i ; i2 = 1 ; k2 = 0 ; j2 = 0
          endif
          
          call cal_u_yz(U1,VS,VY(1:4),VUC,l,l1,m0,m1,i,i2,k,k2,j,j2,i1,k1,j1)
          
          if(l.eq.0) then
            l1 = 2 ; m1 = k ; i2 = 0 ; k2 = 1 ; j2 = 0
          else if(l.eq.1) then
            l1 = 0 ; m1 = i ; i2 = 1 ; k2 = 0 ; j2 = 0
          else if(l.eq.2) then
            l1 = 1 ; m1 = j ; i2 = 0 ; k2 = 0 ; j2 = 1
          endif
          
          call cal_u_yz(U1,VS,VZ(1:4),WUC,l,l1,m0,m1,i,i2,k,k2,j,j2,i1,k1,j1)
          
!CCC  VISCOS TERM OF U CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC D2U/DX2 CCCCCCCCCCCCCCCCCCCCCCCCCC
          VISUX = BIBN2(VX(1),VS,VX(2))
!CCC D2U/DY2 CCCCCCCCCCCCCCCCCCCCCCCCCC
          VISUY = BIBN2(VY(1),VS,VY(2))
!CCC D2U/DZ2 CCCCCCCCCCCCCCCCCCCCCCCCCC
          VISUZ = BIBN2(VZ(1),VS,VZ(2))
!CCCC uDu/Dx
!CCCC vDu/Dy CCCCC
#ifndef DAKU
          FUX = iryu_term(u1(l,nn),VS,VX(1),VX(2),VX(3),VX(4))
          FUY = iryu_term_x(vuc,VS,VY(1),VY(2),VY(3),VY(4))
          FUZ = iryu_term_x(wuc,VS,VZ(1),VZ(2),VZ(3),VZ(4))
#else
          RX  = ( rho(nn)+rho(mn(i-i1,k-k1,j-j1)) ) * HALF
          RVS = RCAL(VS,nn,mn(i-i1,k-k1,j-j1))
          
          RVX(1) = RCAL(VX(1),nn,mn(i+i1,k+k1,j+j1))
          RVX(2) = RCAL(VX(2),mn(i-i1,k-k1,j-j1),mn(i-2*i1,k-2*k1,j-2*j1))
          RVX(3) = RCAL(VX(3),mn(i+i1,k+k1,j+j1),mn(i+2*i1,k+2*k1,j+2*j1))
          RVX(4) = RCAL(VX(4),mn(i-2*i1,k-2*k1,j-2*j1),mn(i-3*i1,k-3*k1,j-3*j1))
          FUX    = iryu_term(u1(l,nn),RVS,RVX(1),RVX(2),RVX(3),RVX(4))
          
          RVY(1) = RCAL2(VY(1),RX)
          RVY(2) = RCAL2(VY(2),RX)
          RVY(3) = RCAL2(VY(3),RX)
          RVY(4) = RCAL2(VY(4),RX)
          FUY    = iryu_term_x(vuc,RVS,RVY(1),RVY(2),RVY(3),RVY(4))
          
          RVZ(1) = RCAL2(VZ(1),RX)
          RVZ(2) = RCAL2(VZ(2),RX)
          RVZ(3) = RCAL2(VZ(3),RX)
          RVZ(4) = RCAL2(VZ(4),RX)
          FUZ    = iryu_term_x(wuc,RVS,RVZ(1),RVZ(2),RVZ(3),RVZ(4))
#endif
#ifdef RAN
          kadoX = ( d(nn) + d(mn(i-i1,k-k1,j-j1)) ) * HALF * ( ONE + ratex ) 
#else
          kadoX = ZERO
#endif
          ! Enter into OX array
          OX1(l,nn) = ( -ONE * ( FUX + FUY + FUZ )                            &
                       + ( VISUX + VISUY + VISUZ * aratez ) * ( nu + kadox )  &
#ifdef DAKU
                       * RX &
#endif
#ifdef REY
                       + ru(nn) &
#endif
                       )
        enddo l_LOOP
      
enddo ADV_LOOP
!$omp end do
!$omp end parallel
return
end subroutine cal_advection_viscosity

   subroutine cal_u_yz(U1,VS,VVX,VUC,l,l1,m0,m1,i,i2,k,k2,j,j2,i1,k1,j1)
        use arrays, only: A, Mn, DX, nff
        use type_list, only: velocity
        use interface_list, only: UUN
        use variables, only: HALF, inne2
        implicit none
        real*8,dimension(0:2,0:inne2+1),intent(in) :: U1
        type(velocity),intent(in) :: VS
        integer,intent(in) :: l, l1, m0, m1, i, i2, k, k2, j, j2, i1, k1, j1
        type(velocity) :: VVX(4)
        real*8 :: dlc, dlp, dlm, VUC
        integer :: mn1, mn2, nnn
    
        dlc        = dx(l1,m1)  
        dlp        = dx(l1,m1+1) 
        dlm        = dx(l1,m1-1)   
        nnn        = mn(i,k,j)
        mn1        = mn(i+i2,k+k2,j+j2)
        mn2        = mn(i+i2-i1,k+k2-k1,j+j2-j1)
        
        VVX(1)%x   = dlc*HALF+dlp*HALF
        VVX(1)%np  = nff(mn(i+i2,k+k2,j+j2))%f
        VVX(1)%nm  = nff(mn(i+i2-i1,k+k2-k1,j+j2-j1))%f
        VVX(1)%npb = nff(mn(i+i2,k+k2,j+j2))%b
        VVX(1)%nmb = nff(mn(i+i2-i1,k+k2-k1,j+j2-j1))%b
        VVX(1)%ap  = a(l1,mn1)
        VVX(1)%am  = a(l1,mn2)
        VVX(1)%v   = UUN(VVX(1),VS,u1(l,mn1),dlc,dlp)


        VVX(2)%x   = -dlc*HALF-dlm*HALF
        VVX(2)%np  = nff(mn(i-i2,k-k2,j-j2))%f
        VVX(2)%nm  = nff(mn(i-i2-i1,k-k2-k1,j-j2-j1))%f
        VVX(2)%npb = nff(mn(i-i2,k-k2,j-j2))%b
        VVX(2)%nmb = nff(mn(i-i2-i1,k-k2-k1,j-j2-j1))%b
        VVX(2)%ap  = a(l1,mn(i,k,j))
        VVX(2)%am  = a(l1,mn(i-i1,k-k1,j-j1))
        VVX(2)%v   = UUN(VVX(2),VS,u1(l,mn(i-i2,k-k2,j-j2)),dlc,dlm)

        VVX(3)%x   = (dlc+dx(l1,m1+2))*HALF+dlp
        VVX(3)%np  = nff(mn(i+i2*2,k+k2*2,j+j2*2))%f
        VVX(3)%nm  = nff(mn(i+i2*2-i1,k+k2*2-k1,j+j2*2-j1))%f
        VVX(3)%npb = nff(mn(i+i2*2,k+k2*2,j+j2*2))%b
        VVX(3)%nmb = nff(mn(i+i2*2-i1,k+k2*2-k1,j+j2*2-j1))%b
        VVX(3)%ap  = a(l1,mn(i+i2*2,k+k2*2,j+j2*2))
        VVX(3)%am  = a(l1,mn(i+i2*2-i1,k+k2*2-k1,j+j2*2-j1))
        VVX(3)%v   = UUN(VVX(3),VVX(1),u1(l,mn(i+i2*2,k+k2*2,j+j2*2)),        &
                         dlp,dx(l1,m1+2))

        VVX(4)%x   = -(dlc+dx(l1,m1-2))*HALF-dlm
        VVX(4)%np  = nff(mn(i-i2*2,k-k2*2,j-j2*2))%f
        VVX(4)%nm  = nff(mn(i-i2*2-i1,k-k2*2-k1,j-j2*2-j1))%f
        VVX(4)%npb = nff(mn(i-i2*2,k-k2*2,j-j2*2))%b
        VVX(4)%nmb = nff(mn(i-i2*2-i1,k-k2*2-k1,j-j2*2-j1))%b
        VVX(4)%ap  = a(l1,mn(i-i2,k-k2,j-j2))
        VVX(4)%am  = a(l1,mn(i-i2-i1,k-k2-k1,j-j2-j1))
        VVX(4)%v   = UUN(VVX(4),VVX(2),u1(l,mn(i-i2*2,k-k2*2,j-j2*2)),        &
                         dlm,dx(l1,m1-2))
        
        VUC = (  u1(l1,mn(i-i1,k-k1,j-j1)) * dx(l,m0)                         & 
               + u1(l1,nnn) * dx(l,m0-1)                                      &     
               + u1(l1,mn(i+i2-i1,k+k2-k1,j+j2-j1)) * dx(l,m0)                &
               + u1(l1,mn(i+i2,k+k2,j+j2)) * dx(l,m0-1) )                     &
              / ( dx(l,m0) + dx(l,m0-1) ) * HALF

#ifdef DAKU
        VVX(1)%nnp = mn1
        VVX(1)%nnm = mn2
        VVX(2)%nnp = mn(i-i2,k-k2,j-j2)
        VVX(2)%nnm = mn(i-i2-i1,k-k2-k1,j-j2-j1)
        VVX(3)%nnp = mn(i+i2*2,k+k2*2,j+j2*2)
        VVX(3)%nnm = mn(i+i2*2-i1,k+k2*2-k1,j+j2*2-j1)
        VVX(4)%nnp = mn(i-i2*2,k-k2*2,j-j2*2)
        VVX(4)%nnm = mn(i-i2*2-i1,k-k2*2-k1,j-j2*2-j1)
#endif
    end subroutine cal_u_yz    

    function UUN(VV,vs,vcc,xc,xcc)
        use interface_list, only: ubnds, ubnds2
        use variables, only: ONE, ZERO
        use type_list, only: velocity   
        implicit none
        real*8 :: UUN
        type(velocity),intent(in) :: VV, vs 
        real*8,intent(in) :: vcc, xc, xcc
        real*8 :: a11 = ONE, vc

        vc = vs%v ! ã´äEì‡ÇÃó¨ë¨Å@vcc ! ã´äEäOÇÃó¨ë¨

        if (VV%np.eq.0.and.VV%nm.eq.0) then
            UUN = vc   ! ÉtÉäÅ[ÉXÉäÉbÉv
        elseif (VV%np.eq.-1.and.VV%nm.eq.-1) then
            if (VV%npb.eq.0.or.VV%nmb.eq.0) then
#ifdef WANKYOKU
                UUN = ubndS(vs%ap,vc,xc,xcc) !   UUN=vc*( (0.78d0-ONE)*vs%ap+ONE )
#else
                UUN = ubndS(vs%ap,vc,xc,xcc)
#endif    
            elseif (VV%npb.eq.2.or.VV%nmb.eq.2) then
                UUN = vcc
            elseif (VV%npb.eq.5.or.VV%nmb.eq.5) then
                UUN = vcc    
            elseif (VV%npb.eq.4.or.VV%nmb.eq.4) then
                UUN = vcc
            elseif (VV%npb.eq.10.or.VV%nmb.eq.10) then
#ifdef WANKYOKU
                UUN = ubnds2(vs%ap,vc,xc,xcc)  !UUN=vc*( (0.78d0-ONE)*vs%ap+ONE )
#else    
#ifdef DRI    
                UUN = vcc
#else
                UUN = ubnds2(vs%ap,vc,xc,xcc)    
#endif
#endif        
            else
                UUN = vc
            endif
        elseif (VV%np.eq.-1.and.VV%nm.ge.0) then
        !    if(VV%npb.ne.10) then
        !    UUN=ubnds(a11,vc,xc,xcc)  !vcc
        !    else
            if (VV%npb.eq.10) then
#ifdef DRI        
                UUN = vcc
#else
                UUN = ubnds2(vs%ap,vc,xc,xcc)  !vcc
#endif
            elseif (VV%npb.eq.4) then
                UUN = vcc        
            else
                UUN = ubnds(vs%ap,vc,xc,xcc)  !vcc    
            endif
        elseif (VV%np.ge.0.and.VV%nm.eq.-1) then
            if (VV%nmb.eq.10) then
#ifdef DRI        
                UUN = vcc
#else
                UUN = ubnds2(vs%ap,vc,xc,xcc)  !vcc
#endif
            elseif (VV%nmb.eq.4) then
                UUN = vcc     
            else
                UUN = ubnds(vs%ap,vc,xc,xcc)  !vcc        
            endif 
        elseif (VV%ap.eq.ZERO.and.VV%am.eq.ZERO) then
#ifdef NAGASA
            UUN = vc*(-ONE)  !ubnds(a11,vc,xc,xcc)
#else
            UUN = ubnds(vs%ap,vc,xc,xcc)
#endif    
        elseif (VV%ap.eq.ZERO.or.VV%am.eq.ZERO) then
#ifdef NAGASA
            UUN = vc*(-ONE)  !vcc 
#else
            UUN = vcc 
#endif 
        else
            UUN = vcc
        endif
    end function UUN
    
    real*8 function third_upwind(uc,Vp,Vc,Vm,Vmm)
        use interface_list,only:first_upwind, u3ss1
        use type_list, only: velocity        
      implicit none
	    real*8::uc  !,u3ss1
	    type(velocity)::Vp,Vc,Vm,Vmm
	    if (abs(Vm%x-Vmm%x).le.abs(Vc%x-Vm%x)*1.2d0) then
  		    third_upwind = uc*u3ss1(Vp,Vc,Vm,Vmm)
	    else
    	    if(vp%x.gt.vm%x) then
        	    third_upwind = first_upwind(uc,vc,vp,Vm)
    	    else
       		    third_upwind = first_upwind(uc,vc,vm,Vp)
    	    endif
	    endif
    end function

    real*8 function BIBN2(Vp,vc,Vm) ! ìÒäKî˜ï™
        use type_list, only: velocity   
        use variables, only: HALF
        implicit none
        type(velocity) :: Vp,vc,Vm
        bibn2 = (&
                (Vp%v-Vc%v)/(Vp%x-Vc%x)-(Vc%v-Vm%v)/(Vc%x-Vm%x)&
                )   &
              / ( (Vp%x-Vm%x) * HALF ) 
    end function

    real*8 function BIBN2f(Vp,vc,Vm,bp,bm) ! ìÒäKî˜ï™
        use type_list, only: velocity   
            use variables, only: ONE, HALF
    implicit none
    type(velocity):: Vp, vc, Vm
    real*8 :: bp, bm, stpara
    stpara = 0.1d0
    bibn2f = (&
            (Vp%v-Vc%v)/(Vp%x-Vc%x)/(ONE+stpara*bp)&
            -(Vc%v-Vm%v)/(Vc%x-Vm%x)/(ONE+stpara*bm)&
            )   &
           / ( (Vp%x-Vm%x) * HALF ) 
    end function

    real*8 function iryu_term_x(vuc,VS,VY1,VY2,VY3,VY4)
        use type_list, only: velocity   
        use interface_list,only: first_upwind, third_upwind
        use variables, only: ZERO
        implicit none
        real*8 :: vuc
        type(velocity) :: VS, VY1, VY2, VY3, VY4
      
        if (vuc.eq.ZERO) then
            iryu_term_x = ZERO
        elseif (vuc.lt.ZERO.and.VY1%np.ge.1.and.VY1%ap.gt.ZERO.and.           & 
                VY1%nm.ge.1.and.VY1%am.gt.ZERO) then
            iryu_term_x = third_upwind(vuc,VY2,VS,VY1,VY3)
        elseif (vuc.gt.ZERO.and.VY2%np.ge.1.and.VY2%ap.gt.ZERO.and.           & 
                VY2%nm.ge.1.and.VY2%am.gt.ZERO) then
            iryu_term_x = third_upwind(vuc,VY1,VS,VY2,VY4)
        else
            iryu_term_x = first_upwind(vuc,VS,VY1,VY2)
        endif
    end function
      
    real*8 function iryu_term(uc,VS,VX1,VX2,VX3,VX4)
        use type_list, only: velocity   
        use interface_list,only: first_upwind, third_upwind
        use variables, only:ZERO
        implicit none
        real*8 :: uc
        type(velocity) :: VS, VX1, VX2, VX3, VX4
        
        if (uc.ge.ZERO.and.VX2%np.ge.1.and.VX4%np.ge.1                        &
            .and.VX2%ap.gt.ZERO.and.VX4%ap.gt.ZERO) then
            iryu_term = third_upwind(uc,VX1,VS,VX2,VX4)
        elseif (uc.lt.ZERO.and.VX1%np.ge.1.and.VX1%ap.gt.ZERO) then
            iryu_term = third_upwind(uc,VX2,VS,VX1,VX3)
        else 
            iryu_term = first_upwind(uc,VS,VX1,VX2)
        endif
    end function

    real*8 function u3ss1(VP,V0,VM,VMM)
        use type_list, only: velocity   
        use variables, only: ZERO, VSMALL
        implicit none
        real*8 :: a11, a12, a13, a21, a22, a23, a31, a32, a33, detem
        real*8 :: y0, y1, y2, y3, b11, b12, b13, as1, as2, u3temp
        type(velocity) :: VP,V0,VM,VMM
        a11 = VP%x ; y0 = VP%v ; y1 = V0%v
        a21 = VM%x ; y2 = VM%v ; a31 = VMM%x ; y3 = VMM%v
               a12 = a11**2
               a22 = a21**2
               a32 = a31**2
               a13 = a11**3
               a23 = a21**3
               a33 = a31**3
               detem = a11*a22*a33&
                     + a21*a32*a13&
                     + a31*a12*a23&
                     - a11*a32*a23&
                     - a21*a12*a33&
                     - a31*a22*a13
              if (detem.ne.ZERO) then
  
                b11 = a22*a33 - a32*a23
                b12 = a32*a13 - a12*a33
                b13 = a12*a23 - a22*a13
                u3temp = (b11*(y0-y1)+b12*(y2-y1)+b13*(y3-y1))/detem
                as1 = (y1-y2)/(ZERO-a21)
                as2 = (y1-y3)/(ZERO-a31)
                if (as1*as2.ge.ZERO) then
                    if (as1*u3temp.le.ZERO) then
                        u3ss1 = as1
                    else
                        u3ss1 = u3temp
                    endif
                else 
                    u3ss1 = u3temp
                endif
              else
                    u3ss1 = - VSMALL     
#ifdef BUGWIN
                write(6,*) 'error in func u3ss1=',u3ss1
#endif
              endif
    end function

#ifdef DAKU
    function RCAL(VX,nnc,nnm)
        use type_list, only: velocity   
        implicit none
        type(velocity) :: RCAL,VX 
        integer :: nnc, nnm
        
        RCAL = VX
        if (rho(nnm).eq.ZERO) then
            RCAL%v = rho(nnc) * VX%v      
        else
            RCAL%v = ( rho(nnm) + rho(nnc) ) * HALF * VX%v
        endif
    end function

    function RCAL2(VX,Rc)
    implicit none
    type(velocity)    :: RCAL2, VX
    real*8,intent(in) :: Rc 
          RCAL2 = VX
          if (rho(VX%nnp).eq.ZERO.and.rho(VX%nnm).eq.ZERO) then
                RCAL2%v = RC * VX%v              
          elseif (rho(VX%nnp).ne.ZERO.and.rho(VX%nnm).eq.ZERO) then
                RCAL2%v = rho(VX%nnp) * VX%v  
          elseif (rho(VX%nnp).eq.ZERO.and.rho(VX%nnm).ne.ZERO) then
                RCAL2%v = rho(VX%nnm) * VX%v  
          else      
                RCAL2%v = ( rho(VX%nnp) + rho(VX%nnm) ) * HALF * VX%v
          endif
    end function
#endif

