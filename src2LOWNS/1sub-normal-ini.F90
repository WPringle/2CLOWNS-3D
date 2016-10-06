!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-normal-ini.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains subroutines to initialise the normal position



subroutine set_initial_surface_position
    use arrays, only: sp, dx, f, in, jn, kn, nff, x, nor, objno, fb0
    use variables, only: inns, inne, mms, sdum, ONE, ZERO, HALF, TINY
    implicit none
    integer :: nfbb, anfbb, i, k, j, nn, mms1
    real*8  :: fff, xx(0:2), dxx(0:2)

    mms1    = 0
    sp(0,:) = sdum
!$omp parallel do private(i,j,k,nfbb,anfbb) reduction(+:mms1)
    do nn = inns,inne
        if (nff(nn)%f.ne.2) then
            if (nff(nn)%f.eq.0.and.f(nn).ne.ZERO) then
                nor(nn,0:1) = ZERO ; nor(nn,2) = ONE
                cycle
            else
                sp(nn,:) = sdum
            endif
        else
            i = in(nn) ; j = jn(nn) ; k = kn(nn)
            nfbb = nff(nn)%b ; anfbb = abs(nfbb)
            dxx = [ dx(0,i), dx(1,j), dx(2,k) ]
            xx  = [ x(0,i),  x(1,j),  x(2,k)  ]
            if (nfbb.eq.-3.and.fb0(nn).lt.ONE) then
                fff = ( ONE - fb0(nn) * ( ONE - f(nn) ) )
            elseif (nfbb.lt.0) then
                fff = f(nn) 
            elseif (nfbb.gt.0) then
                fff = ONE - f(nn)
            endif
            nor(nn,:)         = ZERO 
            nor(nn,anfbb - 1) = ONE * sign(1,-nfbb)
            sp(nn,0:2)        = xx(0:2) + dxx(0:2) * HALF
            sp(nn,anfbb - 1)  = xx(anfbb - 1) + dxx(anfbb - 1) * fff
            where (abs(sp(nn,:)).lt.TINY) sp(nn,:) = ZERO
            nor(nn,3) = - sum( sp(nn,:) * nor(nn,0:2) )
            mms1      = mms1 + 1
        endif
    enddo
!$omp end parallel do
    mms=mms1
end subroutine set_initial_surface_position   
   
subroutine narabe4(L,il)
    real*8,dimension(0:8),intent(in)::L
    real*8,dimension(0:8)::L1
    Integer,dimension(0:8)::il
    integer::ix,ii
    L1=L
    do ix=1,4
        il(ix)=1
        do ii=2,4
            if(L1(ii).lt.L1(il(ix))) then
                il(ix)=ii
            endif
        enddo
        L1(il(ix)) = 1d4
    enddo
end subroutine narabe4
    
subroutine narabe(L,il)
    real*8,dimension(0:8),intent(in)::L
    real*8,dimension(0:8)::L1
    Integer,dimension(0:8)::il
    integer::ix,ii
    L1=L
    do ix=1,8
        il(ix)=1
        do ii=2,8
            if(L1(ii).lt.L1(il(ix))) then
                il(ix)=ii
            endif
        enddo
        L1(il(ix)) = 1d4
    enddo
end subroutine
    
subroutine get_vertex_values(pp,dxx1,nn)
      use arrays, only: in, jn, kn, dx, x, cpl, objno, fb
      use variables, only: ZERO, ONE, ncpls0max, ncpls1max, ncpls2max
      implicit none
      real*8, dimension(0:8,0:2), intent(out):: pp
      real*8, dimension(0:2), intent(out) :: dxx1
      integer, intent(in) :: nn
      real*8:: dx1,dy1,dz1,zb
      integer:: i,k,j

      i = in(nn) ; j = jn(nn) ; k = kn(nn)
      dx1 = dx(0,i)
      dy1 = dx(1,j)
      
      if (fb(nn).eq.ONE                                                 &
#ifdef OBJECT
          .or.OBJNO(nn).ne.0                                            &
#endif                
          ) then
         dz1 = dx(2,k)
         zb  = x(2,k)
      else
         dz1 = dx(2,k) * fb(nn)
         zb  = x(2,k+1) - dz1
      endif      
      !Added Aug 08 2014 for pp in local coordinates
      !Set the origin to transfer local to global coordinates
      pp(0,0:2) = [x(0,i),x(1,j),zb]
      !Make pp in local coordinates
      pp(1,0:2) = [ZERO,ZERO,ZERO]
      pp(2,0:2) = [dx1,ZERO,ZERO]
      pp(3,0:2) = [dx1,dy1,ZERO]
      pp(4,0:2) = [ZERO,dy1,ZERO]
      pp(5,0:2) = [ZERO,ZERO,dz1]
      pp(6,0:2) = [dx1,ZERO,dz1]
      pp(7,0:2) = [dx1,dy1,dz1]
      pp(8,0:2) = [ZERO,dy1,dz1]
      !Calculate the outward normals
      cpl(1,3)= -ONE*sum(cpl(1,0:2)*pp(1,0:2))
      cpl(2,3)= -ONE*sum(cpl(2,0:2)*pp(5,0:2))
      cpl(3,3)= -ONE*sum(cpl(3,0:2)*pp(1,0:2))
      cpl(4,3)= -ONE*sum(cpl(4,0:2)*pp(2,0:2))
      cpl(5,3)= -ONE*sum(cpl(5,0:2)*pp(4,0:2))
      cpl(6,3)= -ONE*sum(cpl(6,0:2)*pp(1,0:2))
      !Get the cell size array
      dxx1(0:2)=[dx1,dy1,dz1]

end subroutine get_vertex_values


