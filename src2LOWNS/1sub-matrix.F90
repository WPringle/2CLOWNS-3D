!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-matrix.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Subroutine to initialise the matrix and functions that it uses
subroutine init_matrix
    use variables, only:  inns, inne, is, ie, js, je, ks, ifu, iob, nl, mp,   &
                          ZERO, HALF
    use arrays, only: in, jn, kn, fb, nff, a, dx, mn, lf, alu_c, ll_c, de_c,  &
                      dim3d
    use interface_list, only: cal_alu_f
    implicit none
    integer :: i,j,k,nnn,l,ii,ll,np,npp,lg(0:2)
    ! This subroutine makes an initial matrix of coefficients that form
    ! a base for the elements calculated every step in cal_matrix_element
    ! to increase speed and efficiency 
    ! (we can just retrieve these coefficients later)
    if (dim3d) then
        nl = 6
        mp = minval(dx, dx > ZERO)**2 
    else
        nl = 4
        if (js.eq.je-1) mp = minval(dx(0:2:2,:), dx(0:2:2,:) > ZERO)**2
        if (is.eq.ie-1) mp = minval(dx(1:2,:), dx(1:2,:) > ZERO)**2
    endif
    allocate(ALU_C(inns:inne,nl),LL_C(inns:inne,nl),DE_C(inns:inne))
    !
    DO nnn = inns,inne
        if (nff(nnn)%f.eq.iob.or.fb(nnn).eq.ZERO) cycle
        i = in(nnn); j = jn(nnn); k = kn(nnn)
        lg(0:2) = [ i, j, k ]
#if defined(DRI) || defined(HENDO)
        fb2 = (fb(nnn)+fbd(nnn)) * HALF
#endif      
        ll = 0 ; DE_C(nnn) = ZERO
        do l = 0,2 ! Loop over dimensions
            if (is.eq.ie-1.and.l.eq.0) cycle ! Cycle when
            if (js.eq.je-1.and.l.eq.1) cycle ! nl = 4
            ll = ll + 1
            do ii = -1,1,2 ! Loop back and forward
                np  = mn(i+ii*lf(l,0),k+ii*lf(l,2),j+ii*lf(l,1))
                npp = mn(i+max(0,ii)*lf(l,0),                                 &
                         k+max(0,ii)*lf(l,2),j+max(0,ii)*lf(l,1))
                if ((nff(np)%f.eq.iob.and.nff(np)%b.ne.4).or.                 &
                     a(l,npp).eq.ZERO.or.fb(np).eq.ZERO) then
                    ! Surrounding cell is a (non-pressure) boundary/object cell
                    LL_C (nnn,2*ll+min(0,ii)) = 0
                    ALU_C(nnn,2*ll+min(0,ii)) = ZERO
                else
                    ! Surrounding cell is a fluid or surface cell

                    LL_C (nnn,2*ll+min(0,ii)) = np
                    ALU_C(nnn,2*ll+min(0,ii)) = -mp *                         &
#if defined(DRI) || defined(HENDO)
                    cal_alu_f(dx(l,lg(l)),dx(l,lg(l)+ii),a2(l,npp),fb2)
#else
#ifndef DAKU
                    cal_alu_f(dx(l,lg(l)),dx(l,lg(l)+ii),a(l,npp),fb(nnn))
#else
                    cal_alu_f_rho(rhon(nnn),rhon(np),dx(l,lg(l)),             &
                                  dx(l,lg(l)+ii),a(l,npp),fb(nnn))
#endif
#endif      
                    DE_C(nnn) = DE_C(nnn) - ALU_C(nnn,2*ll+min(0,ii))
                endif
            enddo
        enddo
    enddo 
endsubroutine init_matrix
    
real*8 function cal_alu_rho(rho,rhom,fnn,dx,dxm)
    use variables, only: HALF, ONE
    implicit none
    real*8,intent(in)::rho,rhom,fnn,dx,dxm
    cal_alu_rho=(fnn-HALF)*dx*rho/(fnn*dx*rho+dxm*HALF*rhom)*(-ONE)
end function cal_alu_rho

real*8 function cal_alu(fnn,dx,dxm)
    use variables, only: HALF, ONE
    implicit none
    real*8,intent(in)::fnn,dx,dxm
    cal_alu = (fnn-HALF)*dx/(fnn*dx+dxm*HALF)*(-ONE)
end function cal_alu

real*8 function cal_alu_f(dx,dxm,ax,fb)
    use variables, only: TWO
    implicit none
    real*8,intent(in)::dx,dxm,ax,fb
    cal_alu_f = TWO / DX / (DX + DXm ) * ax / fb
end function cal_alu_f

real*8 function cal_alu_f_rho(rho,rhom,dx,dxm,ax,fb)
    use variables, only: TWO
    use interface_list ,only : cal_alu_f
    implicit none
    real*8,intent(in)::rho,rhom,dx,dxm,ax,fb
    cal_alu_f_rho=cal_alu_f(dx,dxm,ax,fb)/(rho+rhom)*TWO
end function  cal_alu_f_rho