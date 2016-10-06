#ifdef NORMAL
!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-normal-set_surface_normal2.F90 %%%%%%%%%%%%%%
! Calculates the plane normal based on the gradient of the F function
subroutine set_surface_normal2
    use arrays, only: nor, jn, kn, mn, nff, in, x, sfn, lf, dx
    use variables, only: n, isfn, isfn2, sdum, ZERO, ONE, HALF
    use interface_list, only: crossV, surface_point
    implicit none
    integer :: nnf, nffb, i, k, j, nn, icc,                                   &
               l, ll, lll, lg(0:2), nnn, nn1, nn2, dir, dir1
    real*8 :: nor_sum, s(1:4,0:2), cr_tmp(0:2), nor_tmp(0:3)

! Compute the free surface normal based on the cross-products of the vectors
! from the free surface centroid in the cell
! to the free surface centroids in surrounding cells
!$omp parallel default(none)                                                  &
!$omp shared(n,nff,nor,mn,x,in,jn,kn,lf,dx,isfn2,isfn,sfn)                    &
!$omp private(nn,i,k,j,nffb,l,ll,lll,lg,nnn,nn1,nn2,s,cr_tmp,dir,dir1,        &
!$omp         icc,nor_tmp,nor_sum)
!$omp do
SURF_LOOP: do nnf = 1,isfn
    nn = sfn(nnf) ; nffb = nff(nn)%b
    if (nffb.eq.0) then
        write(6,*) 'nffb=',nffb
        stop
    endif
    i = in(nn) ; j = jn(nn) ; k = kn(nn) ; lg = [i, j, k]  
#ifdef BUGWIN
    if (nn.eq.87390) then
        continue
    endif
#endif
    icc = 0 ; s = ZERO ; lll = abs(nffb) - 1 ; dir = sign(1,-nffb)
    if (lll.eq.1) then
        dir1 = dir*(-1) ! Yï˚å¸ÇÕÅAãtï˚å¸Ç…íTÇ∑ÅB
    else
        dir1 = dir
    endif
    do ll = 1,-1,-2
        do l = abs(min(0,2*dir1)),max(0,2*dir1),dir1
            if (l.eq.lll) cycle
            icc = icc + 1
            nnn = mn(i+ll*lf(l,0),k+ll*lf(l,2),j+ll*lf(l,1))
            nn1 = mn(i+ll*lf(l,0)-dir*lf(lll,0),                              &
                        k+ll*lf(l,2)-dir*lf(lll,2),                           &
                        j+ll*lf(l,1)-dir*lf(lll,1))
            nn2 = mn(i+ll*lf(l,0)+dir*lf(lll,0),                              &
                        k+ll*lf(l,2)+dir*lf(lll,2),                           &
                        j+ll*lf(l,1)+dir*lf(lll,1))
            s(icc,:) = surface_point(nn,nn1,nnn,nn2,ll*(l+1),lg)
        enddo
    enddo
    icc = 0 ; nor_tmp = ZERO
    do l = 1,4
        ll = l + 1 ; if (ll.eq.5) ll = 1
        cr_tmp = crossV(s(l,:),s(ll,:))
        if (sum(cr_tmp**2).ne.ZERO) then
            nor_tmp(0:2) = nor_tmp(0:2) + cr_tmp
            icc = icc + 1
        endif  
    enddo
    if (icc.ne.0) then
        nor_tmp = nor_tmp / dfloat(icc)
        nor_sum = sum(nor_tmp(0:2)**2)
        if (nor_sum.gt.ZERO) then
            nor(nn,0:2) = nor_tmp(0:2) / sqrt(nor_sum)
        else
!$omp critical 
            !write(6,*) n,'nn',nn,'nf=',nff(nn)%f,nff(nn)%b,'icc=',icc
            !write(6,*) 'nor_tmp=',nor_tmp
            !write(6,*) 'spn=',spn(nn,:)
            !write(6,*) s
            !write(6,*) 'Nor sum is <= 0: In surf loop'
            !nor(nn,:) = ZERO
            !nor(nn,abs(nffb)-1) = sign(-ONE,real(nffb))
            call F_gradient_method(i,k,j,lg,nn)
            !stop 'Nor sum is <= 0: Stop in set_surface_normal2'
!$omp end critical 
        endif
    else
        write(6,*) 'nf=',nff(nn)%f,nff(nn)%b,'nn=',nn,'icc=',icc
        stop
    endif   
enddo SURF_LOOP      
!$omp end do

! Now compute free surface normal for 
! air cells with F > 0 and fluid cells with F < 1
!$omp do 
do nnf = isfn+1,isfn2
    nn = sfn(nnf)
    icc =0 ; nor_tmp = ZERO
    i = in(nn) ; j = jn(nn) ; k = kn(nn) ; lg = [i, j, k] 
    ! Search surrounding nf = 2 cells and sum their normals
    do l = 0,2
        do ll = -1,1,2
            nnn = mn(i+ll*lf(l,0),k+ll*lf(l,2),j+ll*lf(l,1))
            if (nff(nnn)%f.eq.2) then
                icc = icc + 1
                nor_tmp(0:2) = nor_tmp(0:2) + nor(nnn,0:2)
            endif
        enddo
    enddo
    if (icc.ne.0) then
        nor_tmp = nor_tmp / dfloat(icc)
        nor_sum = sum(nor_tmp(0:2)**2)
        if (nor_sum.gt.ZERO) then
            nor(nn,0:2) = nor_tmp(0:2) / sqrt(nor_sum)
        else
!$omp critical 
            !write(6,*) n,'nn',nn,'nf=',nff(nn)%f,nff(nn)%b,'icc=',icc
            !write(6,*) 'nor_tmp=',nor_tmp
            !write(6,*) 'spn=',spn(nn,:)
            !write(6,*) 'Nor sum is <= 0: In surf loop'
            call F_gradient_method(i,k,j,lg,nn)
            !stop 'Nor sum is <= 0: Stop in set_surface_normal2'
!$omp end critical 
        endif
    else
        ! OK, no need to calculate normal
    endif   
enddo
!$omp end do
!$omp end parallel
end subroutine set_surface_normal2

function surface_point(nn,ml,mc,mh,l,lg)
    use arrays, only: nff, spn, x, in, jn, kn, dx, f
    use variables, only: sdum, ZERO, HALF
    implicit none
    real*8 :: surface_point(0:2)
    integer,intent(in) :: nn,ml,mc,mh,l,lg(0:2)
    real*8 :: dxc
    integer:: nfc,ll,i

    surface_point = ZERO
    nfc = nff(mc)%f
    ll = abs(l) - 1
    i = lg(ll)
    if (l.gt.0) dxc =   dx(ll,i)
    if (l.lt.0) dxc = - dx(ll,i-1)
                    !Skip when we have a new free surface
    if (nfc.eq.2.and.spn(mc,0).ne.sdum) then
            surface_point(:)  = spn(mc,:) - spn(nn,:)
            surface_point(ll) = surface_point(ll) + dxc
            return
    elseif (nfc.eq.-1) then 
        if (nff(mc)%b.eq.2) then
            surface_point(2)  = f(mc)*dx(2,kn(mc)) - spn(nn,2)
        endif
    elseif (nfc.eq.0) then   !Skip when we have a new free surface
        if (nff(ml)%f.eq.2.and.spn(ml,0).ne.sdum) then
            surface_point(:)  = spn(ml,:) - spn(nn,:)
            surface_point(ll) = surface_point(ll) + dxc
            ! And add a cell length/width in normal direction
            nfc = nff(nn)%b
            ll = abs(nfc) - 1
            i = lg(ll)
            if (nfc.gt.0) dxc =   dx(ll,i)
            if (nfc.lt.0) dxc = - dx(ll,i-1)
            surface_point(ll) = surface_point(ll) + dxc  
            return
        else
 
        endif
    elseif (nfc.eq.1) then    !Skip when we have a new free surface
        if (nff(mh)%f.eq.2.and.spn(mh,0).ne.sdum) then       
            surface_point(:)  = spn(mh,:) - spn(nn,:)
            surface_point(ll) = surface_point(ll) + dxc
            ! And add a cell length/width in normal direction
            nfc = nff(nn)%b
            ll = abs(nfc) - 1
            i = lg(ll)
            if (nfc.gt.0) dxc = - dx(ll,i-1)
            if (nfc.lt.0) dxc = dx(ll,i)
            surface_point(ll) = surface_point(ll) + dxc        
            return
        else

        endif
    endif
    surface_point(ll) = HALF*dx(ll,i+1) - spn(nn,ll) + dxc
end function

subroutine F_gradient_method(i,k,j,lg,nn)
    use variables, only: ONE, ZERO, HALF
    use arrays, only: FB, mn, F, nff, dx, nor, spn, lf
    implicit none
    integer,intent(in) :: i, k, j, lg(0:2), nn
    integer :: icc, l, ll, nnn, nn1, nn2
    real*8  :: nor_tmp(0:2), nor_sum
    ! Method based on F gradient  - William Pringle April 28 2015
    nor_tmp = ZERO
    do l = 0,2
        icc = 0
        do ll = 0,2
            if (l.eq.ll) cycle
            ! Calculate the F gradients at each eight nodes on the plane
            nnn = mn(i+lf(l,0),k+lf(l,2),j+lf(l,1))             
            if (nff(nnn)%f.ne.-1) then
                nn1 = mn(i+lf(ll,0),k+lf(ll,2),j+lf(ll,1))
                nn2 = mn(i+lf(l,0)+lf(ll,0),k+lf(l,2)+lf(ll,2),               &
                         j+lf(l,1)+lf(ll,1))
                if (nff(nn1)%f.ne.-1.and.nff(nn2)%f.ne.-1) then 
                    nor_tmp(l) = nor_tmp(l) + HALF * (FB(nnn)*(F(nnn) - ONE)  &
                               - FB(nn)*(F(nn) - ONE) + FB(nn2)*(F(nn2) - ONE)&
                               - FB(nn1)*(F(nn1) - ONE))
                    icc = icc + 1
                endif
                nn1 = mn(i-lf(ll,0),k-lf(ll,2),j-lf(ll,1))
                nn2 = mn(i+lf(l,0)-lf(ll,0),k+lf(l,2)-lf(ll,2),               &
                         j+lf(l,1)-lf(ll,1))
                if (nff(nn1)%f.ne.-1.and.nff(nn2)%f.ne.-1) then 
                    nor_tmp(l) = nor_tmp(l) + HALF * (FB(nnn)*(F(nnn) - ONE)  &
                               - FB(nn)*(F(nn) - ONE) + FB(nn2)*(F(nn2) - ONE)&
                               - FB(nn1)*(F(nn1) - ONE))
                    icc = icc + 1
                endif
            endif
            nnn = mn(i-lf(l,0),k-lf(l,2),j-lf(l,1))
            if (nff(nnn)%f.ne.-1) then
                nn1 = mn(i+lf(ll,0),k+lf(ll,2),j+lf(ll,1))
                nn2 = mn(i-lf(l,0)+lf(ll,0),k-lf(l,2)+lf(ll,2),               &
                         j-lf(l,1)+lf(ll,1))
                if (nff(nn1)%f.ne.-1.and.nff(nn2)%f.ne.-1) then 
                    nor_tmp(l) = nor_tmp(l) - HALF * (FB(nnn)*(F(nnn) - ONE)  &
                               - FB(nn)*(F(nn) - ONE) + FB(nn2)*(F(nn2) - ONE)&
                               - FB(nn1)*(F(nn1) - ONE))
                    icc = icc + 1
                endif
                nn1 = mn(i-lf(ll,0),k-lf(ll,2),j-lf(ll,1))
                nn2 = mn(i-lf(l,0)-lf(ll,0),k-lf(l,2)-lf(ll,2),               &
                         j-lf(l,1)-lf(ll,1))
                if (nff(nn1)%f.ne.-1.and.nff(nn2)%f.ne.-1) then 
                    nor_tmp(l) = nor_tmp(l) - HALF * (FB(nnn)*(F(nnn) - ONE)  &
                               - FB(nn)*(F(nn) - ONE) + FB(nn2)*(F(nn2) - ONE)&
                               - FB(nn1)*(F(nn1) - ONE))
                    icc = icc + 1
                endif
            endif
        enddo
        ! Now average to the cell center
        if (icc.ne.0) then
            nor_tmp(l) = - nor_tmp(l) / DX(l,lg(l)) / real(icc)
        endif
    enddo
    ! And normalize
    nor_sum = sum(nor_tmp(0:2)**2)
    if (nor_sum.gt.ZERO) then
        nor(nn,0:2) = nor_tmp(0:2) / sqrt(nor_sum)
    else
        write(6,*) 'nf=',nff(nn)%f,nff(nn)%b,'nn=',nn,'icc=',icc
        write(6,*) 'nor_tmp=',nor_tmp
        write(6,*) 'spn=',spn(nn,:)
        write(6,*) 'Nor sum is <= 0: In surf loop'
        nor(nn,:) = ZERO
        if (nff(nn)%b.eq.0) then
            ! Guess
            nor(nn,0:2) = [ZERO,ZERO,ONE]  
        else
            ! Take from nfb direction
            nor(nn,abs(nff(nn)%b)-1) = sign(-ONE,real(nff(nn)%b)) 
        endif
    !   stop 'Nor sum is <= 0: Stop in set_surface_normal2'
    endif
endsubroutine F_gradient_method
#endif    