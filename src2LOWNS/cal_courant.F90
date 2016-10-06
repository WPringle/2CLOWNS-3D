!%%%%%%%%%%%%%%%%%%%% FILE: cal_courant.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the largest local Courant number based on the maximum of
! velocities or the long wave speed
subroutine cal_courant
    use variables, only: icong, inns, inne, courant, g, hasoku_c, is, ie, js, &
                         je, ks, ke, dt, t, velmax, n, inn2d, ZERO, dtm, nu,  &
                         Uweight, Dweight, cr_ul, cr_ml, ce2i, vt_op, TWO,    &
			 HALF, SMALL
    use arrays, only: f, fb, in, jn, kn, mn, dx, u, un, fp, in2d, nff, lf,    &
                      r1, d, dim3d
    !use main_module, only:t_3d
    implicit none
    real*8  :: sumdv, uvv(0:2), uvmax(0:2), dep, turb_c, diff_c
    real*8,allocatable :: turb_ca(:), diff_ca(:)
    integer :: i, j, k, nn, l, lg(0:2), nnu(0:2)
      
    icong = 0 ; nnu = 0
    courant = ZERO  ; uvmax = ZERO  ; sumdv = ZERO
    hasoku_c = ZERO ; turb_c = ZERO ; diff_c = ZERO
!$omp parallel default(none)                                                  &
!$omp   shared(inns,inne,inn2d,in2d,lf,Uweight,Dweight,cr_ul,nu,d,r1,dim3d    &
!$omp          ,in,jn,kn,is,ie,js,je,ks,ke,mn,f,fb,g,un,u,dt,dx,nff,vt_op     &
#ifdef DRI
!$omp          ,nfd &
#endif
!$omp          )    &
!$omp   private  (i,k,j,uvv,l,lg,dep) &
!$omp   reduction(max:uvmax,hasoku_c) &
#ifdef RAN
!$omp   reduction(max:turb_c,diff_c)  &
#endif
!$omp   reduction(+:sumdv)

!$omp do
COURANT_LOOP: do nn = inns,inne
    i = in(nn) ; j = jn(nn) ; k = kn(nn); lg = [i,j,k]
    !Calc Courant in each direction
    do l = 0,2
        if (nff(nn)%f.eq.1.or.nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.1   &
#ifdef DRI
            .or.nfd(i,k,j).eq.1.or.nfd(i-lf(l,0),k-lf(l,2),j-lf(l,1)).eq.1    & 
#endif
            ) then
            uvv(l) = ( max( un(l,nn),ZERO) * dt / dx(l,lg(l))                 &
                   +   max(-un(l,nn),ZERO) * dt / dx(l,lg(l)-1) )
            if (uvv(l).gt.uvmax(l)) uvmax(l) = uvv(l)
        endif
    enddo
#ifdef RAN
    if (nff(nn)%f.ge.1) then
	if (vt_op.ne.2) then
        	! Get Courant number corresponding to turbulent time scale
        	if (r1(nn)%k.gt.SMALL) then
            		turb_c = max( turb_c, dt * r1(nn)%e / r1(nn)%k )
        	endif
	endif
        ! Get Courant number corresponding to diffusion
        if (dim3d) then
            diff_c = max(diff_c, TWO * dt * (d(nn) + nu) *                    &
                        (dx(0,lg(0))**2 * dx(1,lg(1))**2 +                    &
                         dx(0,lg(0))**2 * dx(2,lg(2))**2 +                    &
                         dx(1,lg(1))**2 * dx(2,lg(2))**2) /                   &
                        (dx(0,lg(0))**2 * dx(1,lg(1))**2 * dx(2,lg(2))**2))
        else
            if (js.eq.je-1) then
                diff_c = max(diff_c, TWO * dt * (d(nn) + nu) *                &
                         (dx(0,lg(0))**2 + dx(2,lg(2))**2) /                  &
                         (dx(0,lg(0))**2 * dx(2,lg(2))**2))
            elseif (is.eq.ie-1) then
                diff_c = max(diff_c, TWO * dt * (d(nn) + nu) *                &
                         (dx(1,lg(1))**2 + dx(2,lg(2))**2) /                  &
                         (dx(1,lg(1))**2 * dx(2,lg(2))**2))
            endif
        endif
    endif
#endif
    !Sum change in velocity
    if (nff(nn)%f.eq.1.and.nff(mn(i-1,k,j))%f.ne.0                            &
        .and.nff(mn(i,k-1,j))%f.ne.0.and.nff(mn(i,k,j-1))%f.ne.0 ) then
        sumdv = sumdv + sum(abs(un(:,nn)-u(:,nn)))
    endif
enddo COURANT_LOOP
!$omp end do

!$omp do
!Compute courant from the wave speed
!=== Modified Sep 10 using 2d cell counting system - Will ===========
do nn = 1,inn2d
    i = in2d(0,nn) ; j = in2d(1,nn)
    if (nff(mn(i,ke-1,j))%f.eq.-1) cycle
    dep = ZERO
    do k = ke-1,ks,-1
        if (nff(mn(i,k,j))%f.eq.-1) exit
        if (nff(mn(i,k,j))%f.ge.1) then
            dep = dep + f(mn(i,k,j)) * fb(mn(i,k,j)) * dx(2,k)
        endif
    enddo
    if (dep.gt.ZERO) then
        hasoku_c = max( hasoku_c, sqrt(dep*g) * dt / min(dx(0,i),dx(1,j)) )
    endif
enddo
!$omp end do

!$omp end parallel

if (maxval(uvmax).lt.1.0d-2) then 
    turb_c = ZERO; diff_c = ZERO
endif
courant = max(maxval(uvmax)*Uweight,                                          &
#ifdef RAN
              turb_c*ce2i,diff_c*Dweight,                                     &
#endif
              hasoku_c)
velmax  = sumdv  
if (courant.gt.cr_ul) then
    write(6,*) '(IN VOL2) t=',t,'n=',n,maxval(uvmax)*Uweight,                 &
#ifdef RAN
                turb_c*ce2i,diff_c*Dweight,                                   &
#endif
                'hasoku=', hasoku_c, 'uvmax=', uvmax
                uvmax = ZERO
#ifdef RAN
    allocate(turb_ca(inns:inne),diff_ca(inns:inne))
#endif
    uvmax = ZERO
    do nn = inns,inne
        i = in(nn) ; j = jn(nn) ; k = kn(nn) ; lg = [i,j,k]
#ifdef RAN
        if (nff(nn)%f.ge.1) then
            ! Get Courant number corresponding to turbulent time scale
            turb_ca(nn) = dt * r1(nn)%e / r1(nn)%k 
            ! Get Courant number corresponding to diffusion
            if (dim3d) then
                diff_ca(nn) = TWO * dt * (d(nn) + nu) *                       &
                             (dx(0,lg(0))**2 * dx(1,lg(1))**2 +               &
                              dx(0,lg(0))**2 * dx(2,lg(2))**2 +               &
                              dx(1,lg(1))**2 * dx(2,lg(2))**2) /              &
                             (dx(0,lg(0))**2 * dx(1,lg(1))**2 * dx(2,lg(2))**2)
            else
                if (js.eq.je-1) then
                    diff_ca(nn)  = TWO * dt * (d(nn) + nu) *                  &   
                                   (dx(0,lg(0))**2 + dx(2,lg(2))**2) /        &
                                   (dx(0,lg(0))**2 * dx(2,lg(2))**2)
                elseif (is.eq.ie-1) then
                    diff_ca(nn) = TWO * dt * (d(nn) + nu) *                   &
                                  (dx(1,lg(1))**2 + dx(2,lg(2))**2) /         &
                                  (dx(1,lg(1))**2 * dx(2,lg(2))**2)
                endif
            endif
        else
            turb_ca(nn) = ZERO  ; diff_ca(nn) = ZERO 
        endif
#endif
        do l = 0,2
            if(nff(nn)%f.eq.1.or.nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.1&
#ifdef DRI
                .or.nfd(i,k,j).eq.1.or.nfd(i-lf(l,0),k-lf(l,2),j-lf(l,1)).eq.1& 
#endif
                ) then
                uvv(l) = ( max( un(l,nn),ZERO) * dt / dx(l,lg(l))             &
                       +   max(-un(l,nn),ZERO) * dt / dx(l,lg(l)-1) )
                if (uvv(l).gt.uvmax(l)) then
                    nnu(l)   = nn
                    uvmax(l) = uvv(l)
                endif
            endif
        enddo
    enddo
    if (courant.gt.1.5d0) then
        icong = 6000
        i=in(nnu(2));k=kn(nnu(2));j=jn(nnu(2))
        if (courant.eq.uvmax(2).and.nff(nn)%f.eq.2.and.nff(nn)%b.eq.-3        &
            .and.nff(mn(i,k-1,j))%f.eq.1) then
            nff(nn)%f = 0
            f(nnu(2)) = 0 ;  fp(nnu(2))=0 ; 
            nff(mn(i,k-1,j))%f=2 ; nff(mn(i,k-1,j))%b=-3
        endif
	else
        icong = 600
    endif
        
    write(6,*) '(IN VOL) t=',t,'n=',n,maxval(uvmax),                          &
#ifdef RAN
                turb_c*ce2i,diff_c,                                           &
#endif
                hasoku_c
    write(6,900) courant,icong,uvmax
    write(6,*)  'nn high vel=',nnu
    write(6,*)  'high vel=',u(0,nnu(0)),u(1,nnu(1)),u(2,nnu(2))
    write(6,*)  'high vel_u(ijk)=',in(nnu(0)),jn(nnu(0)),kn(nnu(0)),          &
                 nff(nnu(0))%f,nff(mn(in(nnu(0))-1,kn(nnu(0)),jn(nnu(0))))%f
    if(je-1.ne.js) &
    write(6,*)  'high vel_v(ijk)=',in(nnu(1)),jn(nnu(1)),kn(nnu(1)),          &
                 nff(nnu(1))%f,nff(mn(in(nnu(1)),kn(nnu(1)),jn(nnu(1))-1))%f
    write(6,*)  'high vel_w(ijk)=',in(nnu(2)),jn(nnu(2)),kn(nnu(2)),          &
                 nff(nnu(2))%f,nff(mn(in(nnu(2)),kn(nnu(2))-1,jn(nnu(2))))%f
#ifdef RAN
    nn = maxloc(turb_ca,1)
    write(6,*)  'turb_c*ce2i=',turb_c*ce2i,maxval(r1%k),maxval(r1%e),         &
                in(nn),jn(nn),kn(nn),nn
    nn = maxloc(diff_ca,1)
    write(6,*)  'diff_c=',diff_c,maxval(d),in(nn),jn(nn),kn(nn),nn
    deallocate(turb_ca,diff_ca)
#endif
    if (max(abs(u(0,nnu(0))),abs(u(1,nnu(1))),abs(u(2,nnu(2)))).gt.100d0) then
        icong = 1000
    endif
!    ! Make dt half of what is now
!456 dt = dt * HALF ; dtm = dt
!    courant = courant * HALF
!    write(6,*)  'make dt half: dt =',dt,'Cran=',courant
!    if (courant.gt.cr_ml) goto 456
elseif (courant.gt.0.1d0) then
    icong = 10
endif

! write(6,*) t_3d,f(462357),nff(462357)%f,nff(462357)%b,un(:,462357),courant
    
    
900  format(1h ,'vmax=',e12.5,' icong=',i4,' uvmax=',e12.5          &
                             ,' vvmax=',e12.5,' wvmax=',e12.5)
901  format(1h ,'u=',e12.5,' nf=',i3,' nfb=',i3,' nfxm=',i3,        &
            ' nfbxm=',i3,' i,j,k=',i4,1x,i4,1x,i3,'dx=',f7.4,1x,f7.4)
902  format(1h ,'v=',e12.5,' nf=',i3,' nfb=',i3,' nfym=',i3,        &
            ' nfbym=',i3,' i,j,k=',i4,1x,i4,1x,i3,'dy=',f7.4,1x,f7.4)
903  format(1h ,'w=',e12.5,' nf=',i3,' nfb=',i3,' nfzm=',i3,        &
            ' nfbzm=',i3,' i,j,k=',i4,1x,i4,1x,i3,'dz=',f7.4,1x,f7.4)

endsubroutine
