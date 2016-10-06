!%%%%%%%%%%%%%%%%%% MODULE: TYPE_LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Lists the special variable types
module type_list
    type bcond
        integer            :: North, South, East, West
    end type bcond
    type cellflag
        integer :: f, b, fp, bp
    end type cellflag 
    type rancomp
        real*8 :: k, e, d, x, a
    end type rancomp
    type ranke
        real*8 :: k, e
    end type ranke 
    type velocity
        real*8 :: v, x, am, ap
        integer :: np, npb, nm, nmb, nnp, nnm
    end type velocity
    ! For boundary conditions and wave maker
    type boundary
        character*4        :: North, South, East, West
    end type boundary
    type locinfo
        character*20       :: locname
        integer            :: x_int,y_int
        real*8             :: x_m,y_m
        integer            :: ks,ke,layer
    endtype locinfo
    type inputcond
        character*2        :: wavetype, order
        character*20       :: fmt1
        integer            :: size,width
        real*8             :: dt,T_end,h_wave,v_wave,m_wave,k,                &
                              h_shoki,period,angle,x0(0:1),t_adj
        real*8,allocatable :: hwave(:),vwave(:,:),mwave(:)
    endtype inputcond
    type header
        character*12 :: t, x, eta, depth, u, v, w, mom,                       &
                        tke, e, nu_t, tang, p, area, Mass, Ep, Ek
    endtype header
    
    ! For 2DH calculation    
    type Layer 
    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
    !&&    Each layer contains the relevant information for calculation     &&!
    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
    integer :: GEOC, & !GEOC = 0 : Cartesian coord., GEOC = 1 : Spherical coord.
               SWE , & !SWE  = 0 : Linear eqn.  ,    SWE  = 1 : Non-Linear eqn.
               DISP, & !DISP = 0 : No dispersion.  , DISP = 1 : Explicit disp.
                       !DISP = 2 : Implicit disp. correction
               CD  , & !CD   = 0 : Negative depths are allowed
                       !CD   = 1 : Adjust outgoing fluxes if neg. depth occurs
               BC  , & !BC   = 0 : No breaking condition
                       !BC   = 1 : Breaking condition in dispersion terms
      COUPLING_DIM , & ! = 1 : 1-way coupling , = 2: 2-way coupling
             WRITE_OUT ! = 1 : Write out data , = 0: Do not write out data 
    integer :: DIM,  & !Dimension (1 or 2)
               DEC , & !Number of decimal places for rounding
               INNE, & !Number of calculation cells
               INNE2,& !Number of calculation cells in next outer layer
               INWET,& !Number of wet cells
              INDISP,& !Number of initially wet cells (for matrix calc.)
               IS,IE,& !First and last cell number: X-direction 
               JS,JE,& !First and last cell number: Y-direction 
               KS,KE,& !First and last cell number: Z-direction 
             IS2,IE2,& !First and last cell number in next layer: X-direction 
             JS2,JE2   !First and last cell number in next layer: Y-direction 
    real*8 :: DT,             &     !Time step
              cl_set,         &     !Courant number set (to calculate dt)
              CRAN, VMAX,     &     !Courant number, Maximum velocity
              Dmin,           &     !Minimum Depth where momentum calc. occurs
              MANDEF,         &     !Mannings n (-1d4 for non-uniform value)
              DZ                    !Domain height
    character*2 :: IC  ! Initial condition (IC or GC or NO (none) etc)
    !Integer 3D arrays: order is i,k,j => Only used for visualization
    integer,allocatable,dimension(:,:,:) :: NF,   & !Cell flag (0 = dry,
                                                    !2 = wet, -1 = all object)
                                            NFB  ,& !Free surface orientation
                                                    !(= -3 if nf = 2)
                                              MN2,& !Returns nn inside inner layer
                                        INTER_MAP,& !Used for mapping from the
                                                    !next inner layer
                                        INTER_MAPb  !Used for mapping to the
                                                    !inner layer boundary
    real*8,allocatable,dimension(:,:,:) :: INTER_CF,& !Used for interpolation of
                                                      !next inner layer
                                           INTER_CFb,&!Used for interpolation to
                                                      !inner layer boundary
                                           ETA_IN     !The initial conditions on ETA (or ZK)
    !Integer 2D arrays: order is i,j => nn is used to denote the cell number
    integer,allocatable,dimension(:,:) :: IN,&!Returns i or j for given nn
                                         INW,&!Returns i or j for given nn(wet)
                                         IND,&!Returns i or j for given nn(disp)
                                         IN2,&!Returns i or j inside inner layer
                                          MN,&!Returns nn for given i and j
                                         MND,&!Returns nn(disp) for given i and j
                                          LL  !nn of non-zero non-diagonal
                                              !elements in a row for matrix
                                              !calc. (when DISP = 2)
    !Real 2D arrays: (1,:) => X-direction, (2,:) => Y-direction 
    real*8,allocatable,dimension(:,:) :: U,  &!Depth-averaged velocity 
                                       Umax, &!Max. U at a cell over time 
                                         Q,  &!Current momentum flux (t=n)
                                         Qn, &!New momentum flux (t=n+1)
                                         Qs, &!Interpolated momentum flux (t=n) 
                                         D,  &!Current water depth at cell bound
                                         Dn, &!New water depth at cell boundary
                                         H,  &!Initial water depth at cell bound
                                         DX, &!Cell size
                                         RX, &!Inverse of cell size
                                         RXc,&!Inverse of cell size at center
                                         APA,&!Alpha values (when DISP = 1)
                                         GMA,&!Gamma values (when DISP = 1)
                                         AL   !Non-zero non-diagonal elements
                                              !for matrix calc. (when DISP = 2)
    real*8,allocatable,dimension(:) :: X,Y,Z,&!X Y and Z coordinate vectors
                                       F,    &!Cell volume fraction
                                       ETA,  &!Current free surface (t=n)
                                       ETAn, &!New free surface (t=n+1)
                                      ETAmax,&!Max. ETA at a cell over time
                                       ZK,   &!Current Ground elevation
                                       ZKm,  &!Old Ground elevation
                                       MAN,  &!Mannings n value at each cell
                                       secY, &!sec(latitude)        (GEOC = 1)
                                       sinY, &!sin(latitude) * 2PSI (GEOC = 1)
                                      BREAK, &!Breaking array (for when DISP = 2)
                                       DIAG   !Diagonal elemets for matrix
                                              !calc. (when DISP = 2)
    endtype Layer
    
end module type_list   
!
!%%%%%%%%%%%%%%%%%% MODULE: INTERFACE_LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Lists the various subroutines and functions used 
module interface_list
interface
    function nodeofLineandPlain(A,B,Face)
        real*8,dimension(0:2)::nodeofLineandPlain
        real*8,dimension(0:2),intent(in)::A,B
        real*8,dimension(0:3),intent(in)::face
    end function
!
    subroutine PlainAdjust(nn,vv_1,surf,iz1,iz2,ncpl2,                        &
                           axyz_sd,iz3,pp,iiz1,iiz2,ncpl3)
        integer,intent(in) :: iz1,iz2,iz3,nn,iiz1,iiz2
        real*8,dimension(0:3),intent(in)::surf
        integer,dimension(0:iiz1,-2:iiz2),intent(out)::ncpl3
        integer,dimension(0:iz1,-2:iz2),intent(in)::ncpl2
        real*8,intent(out) :: vv_1
        real*8,dimension(6) :: axyz_sd
        real*8 :: pp(0:iz3,0:2)
    end subroutine PlainAdjust
                           
   
    subroutine PCGPME_WILL(AL,NL,NA,LL,D,N,N2,B,EPS,PARM,ITER)
        integer,intent(in)                 :: N,N2,NA,NL
        real*8,intent(in)                  :: PARM,EPS
        integer,intent(in)                 :: ITER
        real*8,dimension(NA,NL),intent(in) :: AL
        real*8,dimension(N),intent(in)     :: D
        integer,dimension(NA,NL),intent(in):: LL
        real*8,dimension(N),intent(inout)  :: B
    end subroutine PCGPME_WILL
                        
    function cal_UNV(vtc,ipoc)
        implicit none
        real*8::cal_UNV(1:4)
        integer,intent(in)::ipoc
        real*8,intent(in)::vtc(ipoc,1:3)
    end function
!
    function cal_area(pp,ipp)
        real*8::cal_area
        integer,intent(in)::ipp
        real*8,intent(in)::pp(1:ipp,0:2)
    end function cal_area

    subroutine zahyo_setting2(vt1,vt7,pp,dxx1,nn)
        integer,intent(in)::nn
        real*8,dimension(0:2)::vt1,vt7
        real*8,dimension(1:8,0:2),intent(out)::pp
        real*8,dimension(0:2),intent(out)::dxx1
    end subroutine zahyo_setting2
!
    function renzoku_n(nn)
        real*8::renzoku_n
        integer,intent(in)::nn
    end function
!
    function renzoku_p(nn)
        real*8::renzoku_p
        integer,intent(in)::nn
    end function
!
    function cal_drr(nn)
        real*8::cal_drr
        integer,intent(in)::nn
    end function
!
    function cal_drr2(nn)
        real*8::cal_drr2
        integer,intent(in)::nn
    end function
!
    function renzoku_fv(axp,axm,ayp,aym,azp,azm,up,um,vp,vm,wp,wm,dx,dy,dz,fb)
        real*8::renzoku_fv
        real*8,intent(in)::axp,axm,ayp,aym,azp,azm,up,um,vp,vm,wp,wm,dx,dy,dz,fb
    end function
!
    function momentum(uc,ox,pm,pc,dx,f,dt)
    real*8::momentum
        real*8::uc,ox,pm,pc,dx,f,dt
    end function
!
    function momentum2(uc,ox,ox2,pm,pc,dx,f,dt)
        real*8::momentum2
        real*8,intent(in)::uc,ox,ox2,pm,pc,dx,f,dt
    end function
!
    function momentum_d(uc,pm,pc,dx,dt)
        real*8::momentum_d
        real*8::uc,pm,pc,dx,dt
    end function
!
    function momentum_rho(uc,rc,rcn,ox,pm,pc,dx,f,dt)
        real*8::momentum_rho
        real*8,intent(in)::uc,rc,rcn,ox,pm,pc,dx,f,dt
    end function
!
    function momentum_rho2(uc,rc,rcn,ox,ox2,pm,pc,dx,f,dt)
        real*8::momentum_rho2
        real*8,intent(in)::uc,rc,rcn,ox,ox2,pm,pc,dx,f,dt
    end function
!
    function momentum_rho_d(uc,rcn,pm,pc,dx,dt)  !,rc,
        real*8::momentum_rho_d
        real*8,intent(in)::uc,rcn,pm,pc,dx,dt !,rc
    end function  
!
    function rbnd(rc,vvv,dxb,ks_in)
        use type_list,only:rancomp
        type(rancomp)::rbnd
        type(rancomp),intent(in)::rc !,az
        real*8,intent(in)::vvv,dxb,ks_in
    end function rbnd
!
    function rbndNY(rc,vvv,dxb)
        use type_list,only:rancomp
        type(rancomp)::rbndNY
        type(rancomp),intent(in)::rc !,az
        real*8,intent(in)::vvv,dxb
    end function rbndNY
!
    function ustar(uc,dxc,ks_r)
        real*8::ustar
        real*8,intent(in)::uc,dxc,ks_r
    end function 
!
    real*8 function center_diff2(Vp,vc,Vm,anup,anum,dt)
        use type_list,only:velocity
        type(velocity),intent(in)::Vp,vc,Vm
        real*8,intent(in)::anup,anum,dt
    end function
!
    function bibun(up,uc,um,dxp,dxc,dxm)
        real*8::bibun
        real*8,intent(in)::up,uc,um,dxp,dxc,dxm
    end function
!
    function center_diff_r(anu1,rc,rxp,rxm,ryp,rym,rzp,rzm,dt)
        use type_list,only:ranke,rancomp
        type(ranke)::center_diff_r
        real*8,intent(in)::anu1,dt
        type(rancomp),intent(in)::rc,rxp,rxm,ryp,rym,rzp,rzm
    end function
!
    function first_upwind(uc,vc,vp,Vm)
        use type_list,only:velocity
        real*8::first_upwind
        real*8::uc
        type(velocity)::Vp,vc,Vm
    end function
!
    function vbnd(nfm,nfbm,wc,wmm,wmp,dxc,dxm,az,axp,axm)
        real*8::vbnd
        real*8,intent(in)  :: wc,wmp,wmm,dxc,dxm,az,axm,axp
        integer,intent(in) :: nfm,nfbm
    end function vbnd
!
    function vbnd_1(v,va,ax,iflag)
        real*8::vbnd_1
        real*8,intent(in)  :: v,va,ax
        integer,intent(in) :: iflag
    end function vbnd_1
!
    function seikika(vec)
        implicit none
        real*8,dimension(0:2)::seikika   
        real*8,dimension(0:2),intent(in)::vec
    end function
!
    function df_vof(nn)
        integer,intent(in)::nn
        real*8::df_vof
    end function
!      
    function QFF(QX,DXX,nnd,nna,nde)
        real*8::qff
        real*8,intent(in)::QX,DXX
        integer,intent(in)::nnd,nna,nde
    end function
!
    function cal_idou2(nfdb,nfux,ux,xx,yy,zz,dxx,dyy,dzz,nn1)
        real*8::cal_idou2,ux,xx,yy,zz,dxx,dyy,dzz
        integer::nfdb,nfux,nn1
    end function
!
    function crossV(V1,V2)
        real*8,DIMENSION(0:2):: crossV,V1,V2
    end function
!
    function crossL(V1,V2)
        real*8::crossL
        real*8,DIMENSION(0:2):: V1,V2
    end function
!
    subroutine pritime2(date_time_start,date_time_old)
        integer::date_time_old(1:8),date_time_start(1:8)
    end subroutine
!
    function centroid2(PS,n)
        integer,intent(in)::n
        real*8,dimension(0:2)::centroid2
        real*8,dimension(1:n,0:2)::PS
    end function
!
    function surface_point(n1,n2,n3,n4,n5,n6)
        real*8,dimension(0:2) :: surface_point
        integer::n1,n2,n3,n4,n5,n6(0:2)
    end function
!
    function DISTA(vec)
        real*8::DISTA
        real*8,dimension(0:2)::vec
    end function
!
    function ubndS(ax,uc,dxc,dxb)
        real*8::ubndS
        real*8,intent(in)::ax,uc,dxc,dxb
    end function
!
    function ubndS2(anut,uc,dxc,dxb)
        real*8::ubndS2
        real*8,intent(in)::anut,uc,dxc,dxb
    end function
!
    function RHO0T2(Ten,cc,rhos)
        implicit none
        real*8::   RHO0T2
        real*8,intent(in)::   Ten,cc,rhos     
    end function 
!
    function PSHU(PP,PS,DXS,DXP,FS)
        implicit none
        real*8::PSHU
        real*8,intent(in)::PP,PS,DXS,DXP,FS
    end function    
    ! For matrix
    function cal_alu_f(dx,dxm,ax,fb)
        real*8::cal_alu_f
        real*8,intent(in)::dx,dxm,ax,fb    
    end function
    
    function cal_alu(fnn,dx,dxm)
        real*8::cal_alu
        real*8,intent(in)::fnn,dx,dxm
    end function cal_alu
    
    function cal_alu_f_rho(rho,rhom,dx,dxm,ax,fb)
        real*8::cal_alu_f_rho
        real*8,intent(in)::rho,rhom,dx,dxm,ax,fb
    end function  cal_alu_f_rho
    
    function fxx1(ix,ll,lg)
        implicit none
        real*8 :: fxx1
        integer,intent(in) :: ix, ll, lg(0:2)
    end function 
    
    function setV(nff,nffo,Vc,Vo,Vkm,axx,axkm,idir,nnn)
        implicit none
        real*8 :: setV
        real*8,intent(in)  :: Vc, Vo, Vkm, axx, axkm
        integer,intent(in) :: nnn, nff, nffo, idir
    endfunction
    
    function cal_cap(wt,wb,dxx,dtt,cran)
        implicit none
        real*8 :: cal_cap
        real*8,intent(in) :: wt,wb,dxx,cran,dtt
    endfunction   
    
    function hn_sanders_obc(h_wave,h_shoki,eta,mom,g)    
        implicit none
        real*8 :: hn_sanders_obc
        real*8,intent(in) :: h_shoki,h_wave,eta,mom,g
    end function
    
    subroutine sw_eta_vel(eta,vel,period,H,d,t,x,xx,g,kss,kee,z)
        implicit none
        integer,intent(in) :: kss,kee
        real*8,intent(in) :: H,period,t,x,xx,g,d,z(kss:kee)
        real*8,intent(out) :: eta, vel(min(kee-kss,2),kss:kee-1)
    endsubroutine
    
    function calc_sw_per(h_shoki,h_wave,g)
        implicit none
        real*8 :: calc_sw_per
        real*8,intent(in) :: h_shoki,h_wave,g
    end function
    
    function calc_ns_per(h_wave,size,dt)
        implicit none
        real*8 :: calc_ns_per
        integer,intent(in) :: size
        real*8,intent(in) :: dt, h_wave(size)
    end function
    
    function calc_k(period,h,g)
        implicit none
        real*8 :: calc_k
        real*8,intent(in) :: period, h, g
    end function
    
    function set_wl(t,dt,size,hwave,order) 
        implicit none
        real*8 :: set_wl
        integer,intent(in) :: size,order
        real*8,intent(in) :: t,dt
        real*8,intent(in) :: hwave(size)
    end function
    
    function round(num,dec)
        implicit none
        real*8 :: round
        real*8,intent(in) :: num
        integer,intent(in) :: dec
    endfunction round
    
    function LI_1(fm,fp,dxm,dxp)
        implicit none
        real*8 :: LI_1
        real*8,intent(in) :: fm,fp,dxm,dxp
    endfunction LI_1
    
    function UUN(VV,vs,vcc,xc,xcc)
        use type_list, only: velocity   
        implicit none
        real*8 :: UUN
        type(velocity),intent(in)::VV,vs
        real*8,intent(in)::vcc,xc,xcc
    end function UUN
    
    function u3ss1(VP,V0,VM,VMM)
        use type_list, only: velocity   
        implicit none
        real*8:: u3ss1
        type(velocity)::VP,V0,VM,VMM
    end function u3ss1
    
    function third_upwind(uc,Vp,Vc,Vm,Vmm)
        use type_list, only: velocity        
      implicit none
        real*8::third_upwind
	    real*8::uc
	    type(velocity)::Vp,Vc,Vm,Vmm
    end function third_upwind
    
    function iryu_term(uc,VS,VX1,VX2,VX3,VX4)
        use type_list, only: velocity   
        implicit none
        real*8::iryu_term
        real*8::uc
        type(velocity)::VS,VX1,VX2,VX3,VX4
    end function iryu_term

    function iryu_term_x(vuc,VS,VY1,VY2,VY3,VY4)
        use type_list, only: velocity   
        implicit none
        real*8:: iryu_term_x   
        real*8:: vuc
        type(velocity)::VS,VY1,VY2,VY3,VY4
    end function iryu_term_x 

    function BIBN2(Vp,vc,Vm) ! “ñŠK”÷•ª
        use type_list, only: velocity
        implicit none
        real*8::  BIBN2
        type(velocity)::Vp,vc,Vm
    end function BIBN2    

    function UPW_COR(xm,ym,dm,x,y,d,xp,yp,dp,dym,dyp,                         &
                     dymxm,dypxm,rxm,rx,DISP,h,g,dt)
        implicit none
        real*8             :: UPW_COR
        integer,intent(in) :: DISP 
        real*8,intent(in)  :: xm, ym, x, y, xp, yp, dm, d, dp, dym, dyp,      &
                              dymxm, dypxm, rxm, rx, h, g, dt
    endfunction
    
    function CDCUN_2_1(fm,rxm,f,rxp,fp,iflag)
        implicit none
        real*8             :: CDCUN_2_1, f
        real*8,intent(in)  :: fm, fp, rxm, rxp
        integer,intent(in) :: iflag
    endfunction CDCUN_2_1
    
    function wdiv(upper,lower)
        implicit none
        real*8            :: wdiv
        real*8,intent(in) :: upper,lower
    endfunction wdiv 
    
    subroutine find_layer_interp(L,L_out,LN)
        use type_list, only: Layer
        implicit none
        type(Layer),intent(inout) :: L
        type(Layer),intent(in)    :: L_out
        integer,intent(in)        :: LN
    endsubroutine find_layer_interp
    
end interface
end module interface_list
