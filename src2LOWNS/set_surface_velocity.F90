!C*********************************
!C  水面での流速の境界条件設定
!C*********************************
!%%%%%%%%%%%%%%%%%%%% FILE: set_surface_velocity.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets the velocity at the interface between the free surface and air cell
! based on the surrounding values and continuity
! -> iflag = 0 or 1 : The velocity at the boundary normal to the orientation of 
!                     the surface calculated by no gradient condition
! -> iflag = 1      : Velocity is further capped to some specified Courant no.
! -> iflag = 2      : The velocity at the boundary normal to the orientation of 
!                     the surface calculated by continuity
subroutine set_surface_velocity(iflag)
    use interface_list, only: setV, cal_cap
    use variables,only: ks, inne2, inns, inne ,isfn, dt, g, iob, js, icon, is,&
                        ie, ke, je, ZERO, ONE, cr_ml, Uweight
    use arrays,only: un, sfn, in, jn, kn, mn, a, f, dx, fb, sfno, lf, nff
    implicit none
    integer,intent(in) :: iflag 
    integer :: nnf, nn, i, k, j, sf0, l, ll, kd, k1, k2
    integer :: mmm, mmp, mm1, nfm, nfp, lg(0:2)
    
    if (isfn.eq.0) return
!$omp parallel default(none)                                                  &
!$omp shared(isfn,sfn,inns,inne,in,jn,kn,mn,nff,lf,a,un,f,dx,iflag,           &
!$omp        dt,g,sfno,cr_ml,Uweight)                                         &
!$omp private(nnf,nn,i,k,j,l,ll,kd,k1,k2,mmm,mmp,mm1,nfm,nfp,lg,sf0) 
!$omp do schedule(static)
      do nn = inns,inne
        if (nff(nn)%f.eq.0) un(0:2,nn) = ZERO
      enddo
!$omp end do
!$omp do schedule(static)  ! 計算の順番が変わると収束しない。
    SURF_DO: do nnf = 1,isfn
        nn = sfn(nnf) ; i = in(nn) ; j = jn(nn) ; k = kn(nn)
        lg = [ i, j, k]
#ifdef NORMAL
        SF0 = SFNO(nn)
#endif
        ! Get the dimension related to the direction of the surface normal
        if (nff(nn)%b.eq.0) then
            write(6,*) 'i,j,k=',i,j,k
            write(6,*) 'nff,f=',nff(nn),f(nn)
            stop 'Stop in set_surface_velocity: surface cell but nfb = 0'
        endif
        l = abs(nff(nn)%b) - 1    !水面の方向　l=0,1,2 -> x,y,z
        KD = sign(1,nff(nn)%b); K1 = (1-KD)/2; K2 = (1+KD)/2   ! 水面の向き　nnf%b>0 kd=1,k1=0,k2=1  nnf%b<0 kd=-1,k1=1,k2=0
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !  Get all the surrounding velocities except those in the l dimension  !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        do ll = 0,2
            if (ll.eq.l) cycle ! 水面法線方向の流速は考えない
            nfm = nff(mn(i-lf(ll,0),k-lf(ll,2),j-lf(ll,1)))%f    ! 隣接するマイナス方向のセルフラグ
            mmp = mn(i+lf(ll,0),k+lf(ll,2),j+lf(ll,1))           ! 隣接するプラス方向のセル番号
            nfp = nff(mmp)%f                                     ! 隣接するプラス方向のセルフラグ
            ! Get the cell number of the adjacent cell in surface normal 
            ! direction plus one in the ll direct
            mmm = mn(i+kd*lf(l,0)+lf(ll,0),                                   &
                     k+kd*lf(l,2)+lf(ll,2),                                   & ! 隣接するプラス方向のセルの法線の逆方向セル番号
                     j+kd*lf(l,1)+lf(ll,1))
            if (a(ll,mmp).ne.ZERO) un(ll,mmp) = setV(nfp,nfm,un(ll,mmp),      &
                              un(ll,nn),un(ll,mmm),a(ll,mmp),a(ll,mmm),ll+1,nn) 
            ! Get the cell number of the adjacent cell in surface normal dir.
            mmm = mn(i+kd*lf(l,0),k+kd*lf(l,2),j+kd*lf(l,1))
            ! Velocity on the minus cell boundaries
            if (a(ll,nn).ne.ZERO) un(ll,nn) = setV(nfm,nfp,un(ll,nn),         &
                           un(ll,mmp),un(ll,mmm),a(ll,nn),a(ll,mmm),-(ll+1),nn)
        enddo
        !%%%%%%%%%%%%%%%%%%%%%%%% REQUIRED ??? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !            For nearby air cells that have a lot of fluid..          !
        !     Do for all directions or just the surface orientation??         !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        mmm = mn(i-kd*lf(l,0),k-kd*lf(l,2),j-kd*lf(l,1))
        if (nff(mmm)%f.eq.0.and.f(nn).gt.0.99d0.and.f(mmm).gt.ZERO) then
            do ll = 0,2
                if (ll.eq.l) cycle
                if (nff(mn(i-kd*lf(l,0)-lf(ll,0),k-kd*lf(l,2)-lf(ll,2),       &
                     j-kd*lf(l,1)-lf(ll,1)))%f.eq.0.and.a(ll,mmm).ne.ZERO) then
                    un(ll,mmm) = un(ll,nn)
                endif
                mmp = mn(i-kd*lf(l,0)+lf(ll,0),                               &
                         k-kd*lf(l,2)+lf(ll,2),                               &
                         j-kd*lf(l,1)+lf(ll,1))
                if (nff(mmp)%f.eq.0.and.a(ll,mmp).ne.ZERO) then
                    un(ll,mmp) = un(ll,mn(i+lf(ll,0),k+lf(ll,2),j+lf(ll,1)))
                endif
            enddo
        endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !         Calculate velocity in l direction from continuity if         ! 
        !      iflag = 2 or the openness of the boundary is equal to zero      !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mmm = mn(i+k1*lf(l,0),k+k1*lf(l,2),j+k1*lf(l,1))
        mmp = mn(i+k2*lf(l,0),k+k2*lf(l,2),j+k2*lf(l,1))
        if (a(l,mmp).eq.ZERO.or.iflag.eq.2) then
            if (a(l,mmm).eq.ZERO) then
                un(l,mmm) = ZERO
            else
                ! Calculation from continuity
                un(l,mmm) = ZERO
                do ll = 0,2
                    if (ll.eq.l) cycle
                    mm1 = mn(i+lf(ll,0),k+lf(ll,2),j+lf(ll,1))
                    un(l,mmm) = un(l,mmm) + ( a(ll,mm1) * un(ll,mm1)           &
                              - a(ll,nn) * un(ll,nn) ) / dx(ll,lg(ll))
                enddo
                un(l,mmm) = un(l,mmp) * a(l,mmp) + real(kd)                    &
                          * un(l,mmm) * dx(l,lg(l)) / a(l,mmm)
                !if(a(l,mmm).lt.0.1d0) then
           !     if(kd.lt.0.and.un(l,mmm).lt.0.0d0) then
           !         un(l,mmm)=max(un(l,mmm),un(l,mmm))
           !     else if(kd.gt.0.and.un(l,mmm).gt.0.0d0) then
           !         un(l,mmm)=min(un(l,mmm),un(l,mmm))
           !     endif
                !endif
            endif
        else
            ! Else if the openness is non-zero and iflag /= 2
            ! then just get from the nearby cell
            un(l,mmm) = un(l,mmp)
        endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !                 Cap the velocity if iflag = 1                        !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !   if (iflag.eq.1) then
            un(l,mmm) = cal_cap( un(l,mmm), un(l,mmp),                         &
                                 dx(l,lg(l)), dt, cr_ml/Uweight )
     !   endif
        
    enddo SURF_DO
!$omp end do
!$omp end parallel
endsubroutine

function setV(nffc,nffo,Vc,Vo,Vkm,axx,axkm,idir,nnn)
    use variables, only: g, honma, ZERO, ONE, TWO, dt
    use arrays, only: in, jn, kn, a, f, fb, mn, x, dx, nff,   &
                      sfno, pp_s, sp, ncpls, fc, objno, lf, spn
    implicit none
    real*8 :: setV
    real*8,intent(in)  :: Vc, Vo, Vkm, axx, axkm
    integer,intent(in) :: nnn
    integer,intent(in) :: nffc !流速設定境界面を挟むセルフラグ                     ↓流速設定境界面ax  
    integer,intent(in) :: nffo !流速設定境界面と反対側のセルのフラグ　　! 隣接セル nffc !! 対象セル　!! 反対セル nffo  !
    integer,intent(in) :: idir !設定流速の向き　　　　　　　　　　　　　　　         　Vc           Vo       ! (k)
    integer :: ib,ib1,nne,i,j,k,l,kd ,ll                                !        Vkm                   ! (k-1)
    integer,dimension(-2:2) :: imen = [3,6,0,4,5]                    !          ↑一つ下の境界面axkm
    real*8 :: za,zae,hfc,Vupw,Vda !,hn
!-------------------------------------------------------------------------------
! Updated by William Pringle 2014 15 Dec.                                      !
! -> Deleted need for FLATBED                                                  ! 
! -> Only calculate flow using OBJECT method if object or drift is             !
!    present under OBJECT preprocessor                                         ! 
! Updated by William Pringle 2015 29 July.                                     !
! -> Improved algorithm for non OBJ cell when next cell is a surface cell by   !
!    considering the wall height (axx) to not allow fluid to flow through      !
!-------------------------------------------------------------------------------
#ifdef OBJECT
integer :: isfconct_nn,isfconct_nn_o,sf0 !,isfconct_ne
#endif
! Return if next cell is 1 or -1 and 2,3,4,5
if (abs(nffc).eq.1) then
    setV = Vc
    return
endif
l = abs(idir) - 1 
NFF_IF: if (nffc.eq.2) then  !隣接セルが水面の場合．
!隣接セルnffが境界面axの水境界面を持つか調べる
    if (Vc.eq.ZERO.and.idir.lt.0) then
        ! If velocity zero between two fluid cells then set the velocity 
        ! equal to that below + influence due to gravity. idir < 0 is 
        ! just because for idir > 0 this velocity will already be 
        ! calculated at another cell in the main loop
        setV = Vkm - fc(l) * dt           
    else
        ! William 2015/1/2 - Works much better for runup and rundown on
        !                    slope with or without object when simply 
        !                    setting this to the value calculated via 
        !                    momentum eqn.
        !                    i.e. assume continuous flow always
        setV = Vc 
    endif 
!対象隣接セルが空の場合．
elseif (nffc.eq.0) then
    if (axkm.ne.ZERO) then !一つ下の境界面が閉じていない場合
        ! When there exists some fluid cell or such below set the  
        ! velocity equal to that below + influence due to gravity
        setV = Vkm - fc(l) * dt    
#ifdef SATO
        if(nffo.eq.-1 .and. abs(nff(nnn)%b)-1.ne.2 ) then
            setV = (Vkm*fb(nnn)+Vo)/(1.0d0+fb(nnn)) - fc(l) * dt  
    !        if(fb(nnn).lt.0.1d0) write(6,*) nff(nnn)%b ,vkm,vo,fb(nnn),setV
        endif
#endif
        return
    endif
    ! If in the z direction we can now just setV = Vo and return
    if (l.eq.2) then
        setV = Vo - fc(l) * dt   
        return
    endif
    ! If velocity in x or y direction and the next cell is an object or 
    ! surface cell we need to check whether flow is allowed into the next 
    ! cell and we calculate the velocity from homma's weir equation or from Vo
    i = in(nnn) ; j = jn(nnn) ; k = kn(nnn)
    kd = sign(1,idir)  
    !　段落ち式を使うs
#ifdef OBJECT
    sf0  = SFNO(nnn) 
    if (OBJNO(nnn).ne.0                                                        &
#if defined(HENDO) || defined(DRI)
        .or.DRINO(nnn).ne.0                                                    & 
#endif            
        ) then
        if (sf0.ne.0) then
            Vda = ZERO   
            isfconct_nn = 0
            do ib = 1,ncpls(SF0,0,0)
                if (ncpls(SF0,ib,-2).eq.imen(idir)) then
                    isfconct_nn = ib
                endif
            enddo
            if (isfconct_nn.gt.0) then
                hfc = x(2,k);za = x(2,k+1)
                do ib1 = 1,ncpls(SF0,isfconct_nn,0)
                    hfc = max(hfc,pp_s(sf0,ncpls(sf0,isfconct_nn,ib1),2))
                    za  = min(za,pp_s(sf0,ncpls(sf0,isfconct_nn,ib1),2))   
                enddo
                Vda = honma*sqrt(TWO*g*max(hfc-za,ZERO))
            endif
        else
            hfc = f(nnn)*dx(2,k)
            za  = (ONE-axx)*dx(2,k)
            Vda = honma*sqrt(TWO*g*max(hfc-za,ZERO))
        endif
        !反対側の流速を使う
        Vupw = ZERO
        if (real(idir)*Vo.gt.ZERO.and.& !20120326 gt.1 -> ge.1        
            nff(mn(i-kd*lf(l,0),k-kd*lf(l,2),j-kd*lf(l,1)))%f.ge.1) then 
            isfconct_nn_o = 0
            do ib = 1,ncpls(SF0,0,0)
                if (ncpls(SF0,ib,-2).eq.imen(-1*idir)) then
                    isfconct_nn_o = ib
                endif
            enddo
            if (isfconct_nn_o.ne.0) then
                Vupw = Vo
            endif
        endif
        setV = dsign(ONE,dfloat(idir))*max(abs(Vupw),Vda) 
        return
    endif
#endif   
    ! Get the height of the fluid
#ifdef NORMAL
    hfc = spn(nnn,2) - x(2,k)
#else
    hfc = ( f(nnn)*fb(nnn) + (ONE-fb(nnn)) ) * dx(2,k)
#endif
    ! Get the height of the wall
    za = ( ONE - axx ) * dx(2,k)
    if (hfc.gt.za) then
        ! When fluid is higher than the height of the wall
        ! 隣接するセルのセル番号nne
        nne = mn(i+kd*lf(l,0),k+kd*lf(l,2),j+kd*lf(l,1)) 
        zae = ( ONE - fb(nne) ) * dx(2,k)
        za  = ( ONE - fb(nnn) ) * dx(2,k)
        if (zae.lt.za) then
            ! When the adjacent floor level is less than current cell's
            ! get the fluid amount above the floor level
#ifdef NORMAL
            hfc = max(ZERO,spn(nnn,2) - za)
#else
            hfc = f(nnn)*fb(nnn)*dx(2,k)
#endif
            Vda = honma*sqrt(TWO*g*hfc)
        else
            ! When the adjacent floor level is greater than current cell's
            ! subtract the adjacent floor level from the fluid level
            Vda = honma*sqrt(TWO*g*max(hfc-zae,ZERO))    
        endif    
    else
        ! When fluid is lower than the height of the wall
        Vda = ZERO
    endif
    !反対側の流速を使う             !(If there is a wall blocking the flow)
    if (Vo*real(idir).gt.ZERO.and.Vda.gt.ZERO) then
        Vupw = abs(Vo) ! Set velocity equal to that on 
                       ! the opposite cell boundary
    else
        Vupw = ZERO   ! 向きが異なる場合使わない        
    endif   
    setV = dsign(ONE,real(idir))*max(Vupw,Vda)   
    return
else
    stop 'Stop in set_surface_velocity: Unexpected value of NFF'
endif NFF_IF

endfunction   
         
function cal_cap(wt,wb,dxx,dtt,cran)
    use variables, only: ONE
    implicit none
    real*8 :: cal_cap
    real*8,intent(in) :: wt, wb, dxx, cran, dtt
    
    if (abs(wt)*dtt.gt.dxx*cran) then
        cal_cap = wb !sign(ONE,wt) * dxx * cran / dtt 
    else
        cal_cap = wt
    endif  
endfunction
    