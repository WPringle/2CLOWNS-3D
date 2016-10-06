!*********************************
!  計算フラグの設定
!  vof.fの中でCALLされる
!*********************************   
!%%%%%%%%%%%%%%%%%%%% FILE: nfd2015.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assigns the flag to determine a cell type and orientation for surface ones
subroutine set_nflag
    use interface_list, only: fxx1
    use variables, only: is, js, ks, ie, je, ke, isfn, isfn2, ifln,           &
                         inns, inne, n, ZERO, VSMALL, Leps, LITTLE
    use arrays, only: nff, in, jn, kn, f, fb, mn, fln, a, sfn, objno, dx, lf
    implicit none
    integer :: k1, jjs, jje, n2, nloop
    integer :: ix(0:2), imax, kmax, kssu, igai, nnf, dir
    integer :: i, k, j, nn, l
    real*8  :: fxx(0:2)
#ifdef DRI       
    integer :: ib,ib0,ib1,ib2,icpl2,i1,i2,ii
    real*8,dimension(0:2) :: tmpV
    real*8 :: s1
#endif       
!==============================================
!   物体セル  = -1                            
!   気体セル  =  0                            
!   OTHERS   =  1 (FLOUD-SELL TO KATEI)       
!
!   1.物体セルは対象外
!   2.F <  1.0d-12のとき、 F=0,NF=0,NFB=0
!   3.F >= 1.0d-12のとき、 F=1,NF=1,NFB=0
!==============================================
    imax = 0
    kmax = 0
    jjs  = js
    jje  = je-1
!$omp parallel do private(k,i,j) reduction(max:imax,kmax) 
    do nn = inns,inne
        if(nff(nn)%f.eq.-1) cycle    !境界セル以外
        i = in(nn) ; j = jn(nn) ; k = kn(nn)          
        nff(nn)%fp = nff(nn)%f
        nff(nn)%bp = nff(nn)%b  
        if (F(nn).le.VSMALL.or.(f(nn)*fb(nn).lt.1.0d-3.and.          &
            nff(mn(i,k-1,j))%f.eq.-1.and.f(mn(i,k+1,j)).eq.ZERO)) then
            nff(mn(I,k,j))%f = 0                         
            nff(mn(I,k,j))%b = 0                         
            F(nn) = max( ZERO, f(nn) )
#ifdef DRI       
        elseif (nn.le.inne.and.nn.ge.1.and.fb(nn).eq.ZERO) then 
            nff(nn)%f = -1 ; nff(nn)%b = 10
#endif
        else
            nff(nn)%f = 1
            nff(nn)%b = 0
            imax = max(imax,i)
            kmax = max(kmax,k)
        endif                               
    enddo
!$omp end parallel do
    imax = min(ie-1,imax+1)
    kmax = min(ke-1,kmax+1)

!=========================================================
!     1.気体セルに接していれば,nf=2(水表面セル)
!=========================================================
    kssu = ke
3   nloop = 0
4   N2 = 0                                   
    isfn = 0 ; ifln = 0
      
    do k = ks, kmax     
        do j = jjs, jje
            do i = is, imax            
                nn = mn(i,k,j)
                if (nff(nn)%f.eq.-1) cycle   
                if (nff(nn)%f.ge.1) then     
                    !流体セルのとき
                    !Add a fluid cell (surface or full fluid cell)
                    ifln = ifln + 1
                    fln(ifln) = nn       
                    do l = 0,2
                        if (nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.0.and.&
                            a(l,mn(i+lf(l,0),k+lf(l,2),j+lf(l,1))).gt.ZERO)   &
                            goto 122
                        if (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.0.and.&
                            a(l,nn).gt.ZERO) goto 122
                    enddo                      
                    k1 = ke
                    ! Add a full fluid cell
                    nff(nn)%f = 1
                    goto 123
                    ! Add a surface cell
122                 nff(nn)%f = 2   
                    isfn = isfn + 1
                    sfn(isfn) = nn
                    nff(nn)%b = -3 
                    k1 = k
                else
                    k1 = ke
                endif
123             kssu = min(k1,kssu)
            enddo
        enddo
    enddo
!===========================================================
! 以下に該当しない場合 気体セルと見なす
! ただし、f値は,1を越えない限りそのまま。
!===========================================================
    n2 = 0
!$omp parallel do default(none) &
!$omp shared(isfn,sfn,in,jn,kn,nff,n2,f,mn,lf) private(nn,i,k,j,l)  
    do nnf = isfn,1,-1     
        nn = sfn(nnf)
        i = in(nn) ; j = jn(nn) ; k=kn(nn);
    !　 流体セルに接する水面セルを残す   ここから 
     !  底面に接する表面セルは可
        if (nff(mn(i,k-1,j))%f.eq.-1)  cycle 
        do l = 0,2
            if (nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.1) goto 124
            if (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.1) goto 124
            !　流体セルに接する水面セルを残す   ここまで
        enddo
        if (n2.eq.0) N2 = 1
        nff(nn)%f = 0
        nff(nn)%b = 0
124     continue
    enddo
!$omp end parallel do

    if (n2.eq.1) then
        nloop = nloop + 1
        goto 4
    endif
    
!$omp parallel do default(none) &
!$omp shared(isfn,sfn,in,jn,kn,nff,mn,n2,lf) private(nn,i,k,j,l)  
    do nnf = isfn,1,-1     
        nn = sfn(nnf)
        i = in(nn) ; j = jn(nn) ; k=kn(nn);
!     底面に接する表面セルは可
        if (nff(mn(i,k-1,j))%f.eq.-1.and.nff(mn(i,k+1,j))%f.eq.0) cycle
        do l = 0,2
!　          オリジナルの条件　! もっとも厳しい   ここから    
            if (nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.1.and.    &
                nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.0) goto 125
            if (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.1.and.    &
                nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.0) goto 125
!　         オリジナルの条件　! もっとも厳しい   ここまで
        enddo
        if (n2.eq.0) N2 = 1
        nff(nn)%f = 0
        nff(nn)%b = 0
125     continue
    enddo
!$omp end parallel do

    if (n2.eq.1) then
        nloop = nloop + 1
        goto 4
    endif

!=========================================================
!     nf=2のときのnfbの設定
!=========================================================
      igai = 0
!$omp parallel do default(none) &
!$omp  shared(isfn,sfn,in,jn,kn,nff,mn,objno,fb,n,lf,f,a) &
!$omp  private(nn,dir,i,k,j,ix,l,fxx) reduction(+:igai)
    NFB_loop: do nnf = isfn,1,-1
        nn = sfn(nnf)
        i = in(nn) ; j = jn(nn) ; k = kn(nn)      
        nff(nn)%b = 0
        ix = 0
        do l = 0,2
            if (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.1.and.       &
                nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.0.and.       &
                a(l,mn(i+lf(l,0),k+lf(l,2),j+lf(l,1))).gt.ZERO) then
                ix(l) = -1
            elseif (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.0.and.   &
                nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.1.and.       &
                a(l,nn).gt.ZERO) then
                ix(l) = 1
            endif
        enddo
        if (sum(abs(ix)).eq.0) then
!           底面に接する表面セルは可
            if (nff(mn(i,k-1,j))%f.eq.-1.and.nff(mn(i,k+1,j))%f.eq.0) then
#ifndef NORMAL
                ix(2) = -1   
#else
      !-------------------------------------------------------------------------  
      !Added 25, 26 Aug to improve performance of wave front inundation by Will
      !-------------------------------------------------------------------------
#ifdef OBJECT !We don't need to worry when there is an object (e.g. a slope)
                if (OBJNO(nn).ne.0.or.fb(nn).ne.ZERO) then
                    ix(2) = -1
                else
#endif      !If no object then we need to check if ix and iy could be non-zero...
            !(Otherwise the flow will spread too quickly on a flat surface)
                    do l = 0,1
                        if (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.gt.0.and. &
                            nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.0.and. &
                            a(l,mn(i+lf(l,0),k+lf(l,2),j+lf(l,1))).gt.ZERO) then
                            ix(l) = -1
                        elseif (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.0  &
                             .and.nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.gt.0&
                               .and.a(l,nn).gt.ZERO) then
                            ix(l) = 1
                        endif
                    enddo
                    if (sum(abs(ix(0:1))).eq.0) ix(2) = -1
#ifdef OBJECT
                endif
#endif    
#endif
      !--------------------- End addition Will ---------------------------------
            endif
        endif
        if (sum(abs(ix)).eq.1) then
            ! Only one direction has a value; set as direction
            dir = maxloc(abs(ix),1)
            nff(nn)%b = dir * ix(dir-1)
        elseif (sum(abs(ix)).eq.0) then
            ! No direction has a value; something strange
            write(6,*) 'souteigai111',i,k,j,'f=',f(nn),n
            igai = igai + 1
            nff(nn)%f = 0
        else
            ! Two or three directions have a value; find largest F orientation
            do l = 0,2
                if (ix(l).eq.0) then
                    fxx(l) = ZERO ; cycle
                endif
                fxx(l) = fxx1(ix(l),l,[i,j,k])
            enddo   
            dir = maxloc(fxx,1)
            nff(nn)%b = dir * ix(dir-1)
        endif
#ifdef BUGWIN    
        if(nff(nn)%b.eq.0) then
            write(6,*) sum(abs(ix)),ix,nn,nnf
            continue
        endif
#endif
    enddo NFB_loop
!$omp end parallel do

999 if (igai.ne.0) goto 3

#ifdef DRI
    do nnf=1,ifln
    nn=fln(nnf);i=in(nn);j=jn(nn);k=kn(nn);
    if(nff(mn(i,k,j))%f.eq.1.and.f(nn).lt.0.9999d0)then
        if(a(0,nn).gt.ZERO.and.mn(i-1,k,j).gt.0) then
            if(SFNO(mn(i-1,k,j)).gt.0.and.DRINO(mn(i-1,k,j)).gt.0) then
            mms1=SFNO(mn(i-1,k,j))
                do i1=1,ncpls(mms1,0,0)
                if(ncpls(mms1,i1,-1).eq.0.and.ncpls(mms1,i1,-2).eq.4) then
                continue
                tmpV=ZERO
                    ib0=ncpls(mms1,i1,1)
                    do ib=2,ncpls(mms1,i1,0)-1
                    ib1=ncpls(mms1,i1,ib)
                    ib2=ncpls(mms1,i1,ib+1)
                    tmpV=crossV(pp_s(mms1,ib1,0:2)-pp_s(mms1,ib0,0:2),pp_s(mms1,ib2,0:2)-pp_s(mms1,ib0,0:2))+tmpV
                    continue
                    enddo
                    s1=sqrt(sum(tmpV**2))*HALF/dx(1,j)/dx(2,k)  !流体部分の面積率
                    if(a(0,nn)-s1.gt.1.0d-9) then   
                    nff(nn)%f=2
                    nff(nn)%b=nff(mn(i-1,k,j))%b
                    isfn=isfn+1
                    sfn(isfn)=mn(i,k,j)
                    endif 
                endif
                enddo
            continue
            endif
        endif      
    endif
    enddo
#endif
!================================================================!
!     Counting full fluid cells that have a small F value and    !
!      air cells that have a non-trivial F value                 !
!================================================================!
    isfn2 = isfn
!$omp parallel do default(none) &
!$omp shared(isfn2,sfn,in,jn,kn,f,nff,mn,inns,inne,lf) private(nn,i,k,j)  
    do nn = inns, inne
        if ((nff(nn)%f.eq.1.and.f(nn).lt.Leps).or.              &!Fluid cell
            (nff(nn)%f.eq.0.and.f(nn).gt.LITTLE)) then           !Air   cell
            i = in(nn) ; k = kn(nn) ; j = jn(nn)
            do l = 0,2
                if (nff(mn(i-lf(l,0),k-lf(l,2),j-lf(l,1)))%f.eq.2.or.         &
                    nff(mn(i+lf(l,0),k+lf(l,2),j+lf(l,1)))%f.eq.2) then
!$omp critical
                    isfn2 = isfn2 + 1
                    sfn(isfn2) = nn
!$omp end critical
                    exit
                endif
            enddo
        endif
    enddo
!$omp end parallel do
    
endsubroutine
    
function fxx1(ix,ll,lg)
    use arrays, only: mn, dx, f, fb, lf
    use variables, only: ZERO, is, ie, js, je
    implicit none
    real*8 :: fxx1
    integer,intent(in) :: ix, ll, lg(0:2)
    integer :: nnc, l, lll, ii
    fxx1 = ZERO
    do l = 0,2
        if (l.eq.ll) cycle
        if (is.eq.ie-1.and.l.eq.0) cycle
        if (js.eq.je-1.and.l.eq.1) cycle
        do ii = -1,1
            nnc = mn(lg(0) + ix*lf(ll,0) + ii*lf(l,0),                        &
                     lg(2) + ix*lf(ll,2) + ii*lf(l,2),                        &
                     lg(1) + ix*lf(ll,1) + ii*lf(l,1))
            fxx1 = fxx1 + dx(l,lg(l)+ii) * f(nnc) * fb(nnc)
        enddo
        do lll = 0,2
            if (lll.eq.l.or.ll.eq.lll) cycle
            fxx1 = fxx1 * dx(lll,lg(lll))
        enddo
    enddo
    fxx1 = fxx1 * dx(ll,lg(ll)+ix)
endfunction