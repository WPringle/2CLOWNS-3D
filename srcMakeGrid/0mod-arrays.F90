!%%%%%%%%%%%%%%%%%%%%% MODULE: ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Lists the all the arrays used in the program
module arrays
    use type_list, only: cellflag, ranke, Layer, locinfo, boundary
    use variables, only: ZERO
    logical :: one_dim_array= .true. ! 基本は一次元配列。falseのとき三次元配列。
    logical :: xyzn = .true. ! 基本はxyzn。falseのときxyz。
    logical :: dim3d = .true. ! 基本は三次元解析、false のとき二次元解析
    logical :: drift = .false. ! 漂流物の有無　基本はなし。
    logical :: nsdata = .false. ! nsdataの有無　基本はなし。
    logical :: cplusp = .false. ! 地形＋ポリゴンの場合
    integer :: date_time_old(1:8)=0,date_time_start(1:8)=0!　直前計測時間 および　初回計測時間
    integer :: ls(0:2),le(0:2)  ! ls={is,js,ks} le={ie,je,ke}
    integer :: lf(0:2,0:2)

    !From wave module
    type(boundary) :: BC(2:3)
    type(locinfo),allocatable,dimension(:) :: locout
    integer,allocatable,dimension(:) :: filenums
    real*8,allocatable,dimension(:) :: IWD, FSL, FSLmax
    real*8,allocatable,dimension(:,:) :: Flux_3D, U_3D, U_3Dmax
    real*8,allocatable :: eta_in(:,:), vel_in(:,:,:,:)
    
    ! 2DH array !
    type(Layer),allocatable,dimension(:) :: L !The array for each layer
    integer :: ix(2,2) !To add or substract from the correct dimension inside loop
    
    character(len=255),allocatable,dimension(:)::OBJECTshape,DRIFTshape !形状データ名

    real*8,allocatable:: dx(:,:) ! 格子間隔
    real*8,allocatable:: x(:,:) ! 格子座標
    real*8,allocatable,dimension(:):: dxi,dy,dz ! 格子間隔1D表示
    real*8,allocatable,dimension(:):: xi,y,z ! 格子座標1D表示
    integer,allocatable:: nf3(:,:,:)  ,nf3p(:,:,:) ! 格子フラグ,前ステップ格子フラグ
    integer,allocatable:: nfb3(:,:,:) ,nfb3p(:,:,:) ! 格子サブフラグ,前ステップ格子サブフラグ
    type(cellflag),allocatable,dimension(:):: nff !統合格子フラグ
    integer,allocatable:: mn(:,:,:)  ! 三次元位置(i,k,j) -> 一次元位置　nn
    integer,allocatable,dimension(:)::in,jn,kn  ! 一次元位置　nn  ->  三次元位置(i,k,j)
    integer,allocatable,dimension(:,:)::in2D, mn2D !2D cell mapping
    integer,allocatable,dimension(:)::fln  ! 内部セル＋水面セルのリスト　（内部+水面通し番号）ー＞　セル番号
    integer,allocatable,dimension(:)::sfn  ! 水面セルのリスト　（水面通し番号）ー＞　セル番号

    integer,allocatable::ipo(:,:) !各面の頂点数
    integer,allocatable::npo(:,:,:) !各面の頂点のならび
    integer,allocatable::ivt(:),ipl(:) ! 各kdcの総頂点数、総面数
    integer,allocatable::objnoi(:),objno(:) ! 地形セル通し番号　－＞ セル番号、  セル番号－＞地形セル通し番号

    integer,allocatable::sfnoi(:),sfno(:) ! 水面セル通し番号　－＞ セル番号、  セル番号－＞水面セル通し番号
    integer,allocatable::sfno2i(:),sfno2(:) ! 水面通し番号　－＞ セル番号、  セル番号－＞水面通し番号
    integer,allocatable,dimension(:,:,:)::ncpl0,info_0 ! 地形形状情報
    integer,allocatable,dimension(:,:,:)::ncpl1,info_1 ! 漂流物形状情報
    integer,allocatable,dimension(:,:)::ncplc  !セルの頂点の並び
    integer,allocatable,dimension(:,:)::ncpl  !各面内の頂点の並び
    integer,allocatable,dimension(:,:,:)::ncpls ! 水面形状情報

    integer,allocatable,dimension(:)::mmd ! セグメント表面の総数or通し番号　(kdc)
    integer,allocatable,dimension(:,:)::nd ! セグメント各面頂数　 (kdc,mmdmx)
    integer,allocatable,dimension(:,:,:)::IDP
                                        ! セグメントの各面の情報　 (kdc,mmd(kdc),1:2) -> (含まれるセル番号:全体での面番号）
    integer,allocatable,dimension(:)::mmd2 ! 二次元可視化時のセグメント面総数　 (kdc,mmd2max)
    integer,allocatable,dimension(:,:)::nd2 ! 二次元可視化時のセグメント各面頂点数　 (kdc,mmd2max)
    integer,allocatable,dimension(:,:,:)::IDP2 ! 二次元可視化時のセグメント各面の情報　 (kdc,mm2,1:2)

    integer,allocatable,dimension(:)::nds ! 各水面の頂点数　（水面通し番号）

    integer,allocatable,dimension(:)::idri_kdc ! (ndri_face)  漂流物表面通し番号 -> 漂流物番号の関係
    integer,allocatable,dimension(:)::idri_fc_n ! (ndri_face)  漂流物表面通し番号 -> セル内面番号
    integer,allocatable,dimension(:,:)::idri_fc_all_n !  漂流物番号,セル内面番号 ->漂流物表面通し番号
    integer,allocatable,dimension(:,:)::men ! セルの面情報
    integer,allocatable,dimension(:,:)::path0 ! 辺のつながり
    integer,allocatable,dimension(:,:)::side     ! 格子の辺の定義

    real*8,allocatable,dimension(:,:):: u,un ! 流速、次ステップ流速
    real*8,allocatable,dimension(:,:,:):: ui,v,w ! 流速、次ステップ流速(三次元配列)
    real*8,allocatable,dimension(:,:):: qf ! 移動水塊体積[m/s]
    real*8,allocatable,dimension(:,:,:):: ox ! adv+vis
    real*8,allocatable,dimension(:):: f,fp,p ! 充填率、次ステップ充填率、圧力
    real*8,allocatable,dimension(:,:,:):: f3 ! 充填率３次元配置
    real*8,allocatable,dimension(:,:,:):: p3 ! 圧力３次元配置
    real*8,allocatable,dimension(:,:):: A,A0,AD ! 開口率、開口率読み込み値
    real*8,allocatable,dimension(:,:,:):: ax,ay,az ! 開口率３次元
    real*8,allocatable,dimension(:):: fb,fb0,fbd ! 空隙率、空隙率読み込み値
    real*8,allocatable,dimension(:,:,:):: fb3 ! 空隙率３次元
    real*8,allocatable,dimension(:,:,:)::vt !  頂点座標
    real*8,allocatable,dimension(:,:)::vtplus !  頂点の平行移動
    real*8,allocatable,dimension(:)::vtplus_s !  頂点の平行移動
    type(ranke),allocatable,dimension(:)::r,r1 !  乱流エネルギー

    real*8,allocatable,dimension(:) :: d,pro,prop,SS !  渦動粘性、乱流エネルギー生成項, 乱流エネルギー生成項(previous),strain-rate
    real*8,allocatable,dimension(:) :: AsUs 
    real*8,allocatable,dimension(:,:)::sp,spn ! 水面の重心,次ステップ
    real*8,allocatable,dimension(:,:)::nor ! 水面の方程式

    real*8,allocatable,dimension(:)::rhod ! 物体密度
    real*8,allocatable,dimension(:)::md ! 物体質量
    real*8,allocatable,dimension(:,:,:)::pp_0 ! 地形形状の頂点
    real*8,allocatable,dimension(:,:,:)::pp_1 ! 漂流部形状の頂点
    real*8,allocatable,dimension(:,:,:)::pp_s ! 水面形状の頂点
    real*8,allocatable,dimension(:,:,:)::xxmax,xxmin !物体表面の(xmax,ymax,zmax)および(xmin,ymin,zmin)
    real*8,allocatable,dimension(:,:)::vtmx,vtmn !物体の(xmax,ymax,zmax)および(xmin,ymin,zmin)

    real*8,allocatable,dimension(:,:,:,:)::dp  ! セグメント各面の頂点情報
    real*8,allocatable,dimension(:,:,:,:)::dp2 ! 二次元可視化時のセグメント各面の頂点情報
    real*8,allocatable,dimension(:,:,:)::surfp !　水面の頂点　（水面番号、並び、座標）
    real*8,allocatable,dimension(:,:)::cplc     !(1:6,0:3) 格子面の平面方程式
    real*8,allocatable,dimension(:,:)::cpl      !(1:mmp,0:3) 全ての平面方程式
    real*8,allocatable,dimension(:,:,:)::upl !(kdc,面番号,1:4) 物体各面の方程式（外向き正）
    real*8,allocatable,dimension(:,:)::upl0      !  セル面の単位法線ベクトル

    !allocate(pp_0(1:mm0mx,0:mcipmx,0:2),ncpl0(mm0mx,0:ncpl00max,-2:mcfvtmx),info_0(mm0mx,1:2,1:mcipmx))
    real*8,allocatable,DIMENSION(:,:)::c !,cfx,cfy,cfz
    real*8,allocatable,DIMENSION(:,:,:)::cf
    real*8,allocatable,DIMENSION(:)::dia
    real*8,allocatable,DIMENSION(:)::te
    real*8,allocatable,DIMENSION(:,:)::CN
    real*8,allocatable,DIMENSION(:)::TEN
    real*8,allocatable,DIMENSION(:)::cmax,cmin

    real*8,allocatable,DIMENSION(:)::rho,rhon  !,rhow
    real*8 :: fc(0:2) ! 外力

    real*8,allocatable,dimension(:,:)::GD,GDP !(mdri,3)
    real*8,allocatable,dimension(:,:)::vvd,vvdp
    real*8,allocatable,dimension(:,:)::acc !(mdri,3)
    real*8,allocatable,dimension(:,:)::FFD,FFDp,rot,rotp !(mdri,3)
    real*8,allocatable,dimension(:,:)::roll_radius,omega,omegap,domega,domegap,theta,thetap

    !Matrix ones
    real*8,allocatable,dimension(:)   :: b, de, de_c, xe
    real*8,allocatable,dimension(:,:) :: alu, alu_c
    integer,allocatable :: mno(:,:,:), mno2mn(:)
    integer,allocatable,dimension(:,:) :: llu, LL_C

    real*8,allocatable,dimension(:,:)::Iner  ! Iner(kdc,xyz方向):漂流物の慣性モーメント(mdri,1:3)
    real*8,allocatable,dimension(:,:)::Inerbym  ! Iner(kdc,xyz方向):漂流物の慣性モーメント(mdri,1:3)
    !======================= 2014/05/19 =======================================================================
    real*8,allocatable,dimension(:,:,:)::Iner2,Iner2p !Iner2(kdc,3*3):漂流物の慣性モーメント(mdri,1:3,1:3)
    real*8,allocatable,dimension(:,:,:)::INV_iner2 !Iner2(kdc,3*3):漂流物の慣性モーメント(mdri,1:3,1:3) inverse matrix
    !===========================================================================================================
    real*8,allocatable,dimension(:,:)::dgd   !(mdri,1:3)
    integer,allocatable,dimension(:)::idm  ! 漂流物が陸上に静止しているとき0，動いているとき1 (mdri)
    real*8,allocatable,dimension(:,:)::uplb  ! 解析領域の表面(1:mp,1:4)

    !real*8,allocatable,dimension(:,:)::mftan,mftanp  !係留策による接線方向の力
    real*8,allocatable,dimension(:,:)::dtheta  ! theta(i):i軸回りの回転角 (mdri,3)

    real*8,allocatable,dimension(:,:,:)::DVECp  ! DVEC(kdc,i,1～3):i軸の方向ベクトル !(mdri,1:3,1:3)
    real*8,allocatable::fix(:),fixp(:,:)  !回転軸
    !integer,allocatable::drinoi(:),drino(:) ! 漂流物を部分的に含むセル通し番号　－＞ セル番号、  セル番号－＞漂流物を部分的にセル通し番号
    !integer,allocatable::drino2i(:),drino2(:,:) ! 漂流物内部セル通し番号　－＞ セル番号、  （セル番号,漂流物番号）－＞漂流物内部セル通し番号
    integer,allocatable,dimension(:)::DRINO,DRINOP,DRINOK !漂流物を含むセルの通し番号
    integer,allocatable,dimension(:,:)::DRINO2,DRINO2K,DRINO2p !漂流物内部セルの通し番号
    integer,allocatable,dimension(:)::DRINOI,DRINO2I,DRINOIK,DRINO2IK

        !!!!!!!!!!!
        !real*8,allocatable::fix(:),fixp(:,:)  !回転軸
        !!!!!!!!!!!    
!#ifdef DRI2        
!        integer,dimension(mnh)::DRI2PNOI,DRI2PNOIP  !(nn)  分断セルの通し番号からセル通し番号の逆引き
!        integer,dimension(100,3)::DRI2PNO !(mm2p,1) 当該セルのセルの通し番号nn, (mm2p,2) 分断数
!        real*8,dimension(100,10,6)::axyz,axyzd,axyz2  !分断された各部分の開口率
!        real*8,dimension(100,10,6)::Vxyz,Vxyzn  !分断された各部分の流速        
!        real*8,dimension(100,10)::P2P  !分断された各部分の圧力        
!        real*8,dimension(100,10)::FB2P,FB2PD,FB2P2,F2P !(mm2p,idvd) 分断された各部分の空隙率  
!        integer,dimension(100,10)::NF2P !(mm2p,idvd) 分断された各部分のセルフラグ              
!#endif

    contains

!%%%%%%%%%%%%%%%%%%%%% SUBROUTINES: ALLOCATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> The following subroutines allocate the sizes of the arrays listed above
    subroutine alloc_drift(ndri)
        implicit none
        integer,intent(in)::ndri
        allocate(md(1:ndri))
        allocate(Iner(ndri,3),Iner2(ndri,3,3),Iner2p(ndri,3,3), inv_iner2(ndri, 3, 3))
#ifdef FLAPGATE    
        allocate(fix(ndri))                            !!!!!!!!!!!!回転中心判定
        allocate(fixp(ndri,3))                       !!!!!!!!!!!!回転中心となる点   
#endif    
        allocate(dtheta(ndri,3),omega(ndri,3),vvd(ndri,3),domega(ndri,3),domegap(ndri,3));dtheta=ZERO;omega=ZERO;vvd=ZERO;domega=ZERO
        allocate(omegap(ndri,3),vvdp(ndri,3),roll_radius(ndri,3)) ;vvdp=ZERO;omegap=0;roll_radius=0
        allocate(rot(ndri,3),rotp(ndri,3),ffdp(ndri,3),ffd(ndri,3),acc(ndri,3));rot=ZERO;rotp=ZERO;ffd=ZERO;ffdp=ZERO;acc=ZERO
        allocate(GDp(ndri,3),thetap(ndri,3),dGD(ndri,3))   ;GDp=ZERO  ;thetap=0 ;dGD=ZERO
        allocate(DVECp(1:ndri,1:3,1:3)) ;DVECp=ZERO
        !

        !
    end subroutine

    subroutine alloc_quvwf(mnh1)
        implicit none
        integer::mnh1
        allocate(qf(0:2,0:mnh1+1))
    end subroutine   

    subroutine alloc_dxdydz(mx0,my0,mz0,mx1,my1,mz1)
        implicit none
!        logical::oned
        integer::mx0,my0,mz0,mx1,my1,mz1,mnh1,mx2
        mx2=max(mx1,my1,mz1)
        mnh1=(mx1-mx0+1)*(my1-my0+1)*(mz1-mz0+1)
    
    if(one_dim_array) then
        allocate(dx(0:2,mx2));dx=ZERO
        allocate(x(0:2,mx2));x=ZERO
    else
        allocate(dxi(mx0-2:mx1),dy(my0-2:my1),dz(mz0-2:mz1));dxi=ZERO;dy=ZERO;dz=ZERO
        allocate(xi(mx0-2:mx1),y(my0-2:my1),z(mz0-2:mz1));xi=ZERO;y=ZERO;z=ZERO
    endif
    !    allocate(nf(mx0-2:MX1,mz0-2:mz1,my0-2:my1),nfb(mx0-2:MX1,mz0-2:mz1,my0-2:my1))
    !    allocate(nfbp(mx0-2:MX1,mz0-2:mz1,my0-2:my1),nfp(mx0-2:MX1,mz0-2:mz1,my0-2:my1))
        allocate(mn(mx0-2:MX1,mz0-2:mz1,my0-2:my1));mn=0
        allocate(nff(0:mnh1));nff(:)%f=-1;nff(:)%b=0
        allocate(in(Mnh1),jn(Mnh1),kn(Mnh1))
    end subroutine

    subroutine alloc_uvwfp(mnh1)
        use variables, only: TOME
        implicit none
        integer :: mnh1
        allocate(u(0:2,0:mnh1+1),un(0:2,0:mnh1+1));u=ZERO
        allocate(f(0:mnh1+1),fp(0:mnh1+1),p(0:mnh1));f=ZERO;p=ZERO    
        allocate(ox(TOME,0:2,0:mnh1));ox = ZERO !;ox1=0.0d0ox1(0:2,0:mnh1));
        allocate(sfn(Mnh1),fln(Mnh1))
    end subroutine alloc_uvwfp

    subroutine alloc_a0fb0(mnh1,mnh2)
        implicit none
        integer::mnh1,mnh2
        allocate(a0(0:2,0:mnh1+1),fb0(0:mnh2+1));a0=ZERO;fb0=ZERO
    end subroutine alloc_a0fb0

    subroutine alloc_afb(mnh1,mnh2)
        implicit none
        integer::mnh1,mnh2
        allocate(A(0:2,0:mnh1+1),fb(0:mnh2+1));a=ZERO;fb=ZERO
    end subroutine alloc_afb

    subroutine alloc_adfbd(mnh1,mnh2)
        implicit none
        integer::mnh1,mnh2
        allocate(Ad(0:2,0:mnh1+1),fbd(0:mnh2+1));ad=ZERO;fbd=ZERO
    end subroutine alloc_adfbd

    subroutine alloc_aked(mnh1)
    use variables, only: vt_op
    implicit none
    integer::mnh1
    ! Turbulent viscosity and production
    allocate(d(0:mnh1+1),pro(0:mnh1+1),prop(0:mnh1+1))
    allocate(r1(0:mnh1+1),r(0:mnh1+1)) 
    if (vt_op.eq.2) then
        ! Smagorinsky
        allocate(ss(0:mnh1+1))
    elseif (vt_op.eq.3) then
        ! Realizable k-eps
        allocate(ss(0:mnh1+1),AsUs(0:mnh1+1)) 
    endif
#ifdef FURYOKU
    allocate(frate(0:mnh1))
#endif
#ifdef REY
    allocate(ru(0:mnh1));allocate(ru(1,0:mnh1));allocate(ru(2,0:mnh1))
#endif
    end subroutine

    subroutine alloc_for_segment(ndri,mm0mx,mmdmx,mmd2mx,ncpl02max,inne2)
      !  use variables
        implicit none
        integer,intent(in)::ndri,mm0mx,mmdmx,mmd2mx,ncpl02max,inne2
     !   if(.not.allocated(OBJNO)) then
     !   allocate(OBJNOi(1:mm0mx))
     !   allocate(OBJNO(0:inne2))
     !   endif
        allocate(vtmx(1:ndri,0:2),vtmn(1:ndri,0:2))
        allocate(mmd(1:ndri),mmd2(1:ndri))
     !   allocate(IDP(1:ndri,1:mmdmx,1:2))
     !   allocate(dp(1:ndri,1:mmdmx,1:ncpl02max,1:3))
     !   allocate(nd(1:ndri,1:mmdmx))
     !   allocate(IDP2(1:ndri,1:mmd2mx,1:2))
     !   allocate(dp2(1:ndri,1:mmd2mx,1:ncpl02max,1:3))
     !   allocate(nd2(1:ndri,1:mmd2mx))
    
    end subroutine alloc_for_segment

    subroutine alloc_norsp(mnh1)
          implicit none
          integer::mnh1
          allocate(sp(0:mnh1+1,0:2));allocate(spn(0:mnh1+1,0:2))
          allocate(nor(0:mnh1,0:3))
    end subroutine

    subroutine alloc_uivwp3(ie1,je1,ke1)
    implicit none
    integer::ie1,ke1,je1
    allocate(ui(1:ie1+1,1:ke1+1,1:je1+1))
    allocate(v(1:ie1+1,1:ke1+1,1:je1+1))
    allocate(w(1:ie1+1,1:ke1+1,1:je1+1))
    allocate(p3(1:ie1+1,1:ke1+1,1:je1+1))
    allocate(f3(1:ie1+1,1:ke1+1,1:je1+1))
    ui=ZERO;v=ZERO;w=ZERO;p3=ZERO;f3=ZERO
    end subroutine alloc_uivwp3

    subroutine alloc_fnfnfb(ie1,je1,ke1)
    implicit none
    integer::ie1,ke1,je1
    if(.not.allocated(f3)) then
        allocate(f3(1:ie1+1,1:ke1+1,1:je1+1))
        f3=ZERO
    endif
    allocate(nf3(1:ie1+1,1:ke1+1,1:je1+1),nfb3(1:ie1+1,1:ke1+1,1:je1+1))
    end subroutine

    subroutine alloc_axyzfb3(ie1,je1,ke1)
    implicit none
    integer::ie1,ke1,je1
    allocate(ax(1:ie1+1,1:ke1+1,1:je1+1),ay(1:ie1+1,1:ke1+1,1:je1+1),az(1:ie1+1,1:ke1+1,1:je1+1),fb3(1:ie1+1,1:ke1+1,1:je1+1))
    ax=ZERO;ay=ZERO;az=ZERO;fb3=ZERO
    end subroutine
    subroutine alloc_c(mfk1,mnh1)
    integer::mnh1,mfk1
    allocate(c(0:mfk1,0:mnh1),cf(0:2,0:mfk1,0:mnh1))  !,cfy(0:mfk1,0:mnh1),cfz(0:mfk1,0:mnh1))
    allocate(cn(0:mfk1,0:mnh1))
    c=ZERO;cn=ZERO;cf=ZERO  !;cfy=ZERO;cfz=ZERO
    end subroutine alloc_c

    subroutine alloc_te(mnh1)
    implicit none
    integer::mnh1
    allocate(te(0:mnh1))
    allocate(ten(0:mnh1))
    te=ZERO;ten=ZERO
    end subroutine alloc_te

    subroutine alloc_rho(mnh1)
    implicit none
    integer::mnh1
    allocate(rho(0:mnh1))
    allocate(rhon(0:mnh1))
    rho=ZERO;rhon=ZERO
    end subroutine
    ! For matrix elements
    subroutine dealloc_matrix
    deallocate(b)
    deallocate(de)
    deallocate(xe)
!    deallocate(dd)
!    deallocate(wk)
    deallocate(alu)
    deallocate(mno)
    deallocate(llu)
    deallocate(mno2mn)
    end subroutine

    subroutine alloc_matrix(mm0,mx1,my1,mz1)
    use variables, only: mm2, nl
    implicit none
    integer,intent(in)::mm0,mx1,my1,mz1
    mm2 = mm0
    allocate(b(mm0))
    allocate(de(mm0))
    allocate(xe(0:mm0))
!    allocate(dd(0:mm0))
!    allocate(wk(0:mm0,5))
    allocate(alu(mm0,nl))
    allocate(mno(mx1,mz1,my1))
    allocate(llu(mm0,nl))
     allocate(mno2mn(1:mm0))
    end subroutine
    
    ! 2DH allocation    
    subroutine alloc_XY(LN,mn2Dx,mn2Dy,ks2D,ke2D)
        integer,intent(in) :: LN,mn2Dy,mn2Dx,ks2D,ke2D
        !Allocate the arrays with multi-dimensions
        allocate(L(LN)%Z(ks2D-2:ke2D+1),L(LN)%X(mn2Dx),L(LN)%Y(mn2Dy))
        allocate(L(LN)%MN(mn2Dx,mn2Dy),L(LN)%MND(mn2Dx,mn2Dy))
        allocate(L(LN)%NF(mn2Dx,ks2D-1:ke2D,mn2Dy),&
                 L(LN)%NFB(mn2Dx,ks2D-1:ke2D,mn2Dy))
        ! The spherical coordinate arrays
        if (L(LN)%GEOC.eq.1) allocate(L(LN)%secY(mn2Dy),L(LN)%sinY(mn2Dy))
        IX(1,:) = [ 1, 0 ]
        IX(2,:) = [ 0, 1 ]
    endsubroutine alloc_XY
!
    subroutine alloc_F(LN,DIM,mn2Dh)
    integer,intent(in) :: LN,DIM,mn2Dh
        !Allocate arrays of single dimension (or with 2 for x and y directions)
        allocate(L(LN)%MAN(0:mn2Dh),L(LN)%ZK(0:mn2Dh),L(LN)%ZKm(0:mn2Dh))
        allocate(L(LN)%ETA(0:mn2Dh),L(LN)%ETAn(0:mn2Dh),L(LN)%F(0:mn2Dh))
        allocate(L(LN)%INW(2,0:mn2Dh),L(LN)%IND(2,0:mn2Dh))
        allocate(L(LN)%RX(DIM,0:mn2Dh),L(LN)%RXc(DIM,0:mn2Dh))
        allocate(L(LN)%DX(DIM,0:mn2Dh),L(LN)%H(DIM,0:mn2Dh))
        allocate(L(LN)%Q(DIM,0:mn2Dh),L(LN)%Qn(DIM,0:mn2Dh))
        allocate(L(LN)%Dn(DIM,0:mn2Dh),L(LN)%D(DIM,0:mn2Dh))
        allocate(L(LN)%U(L(LN)%DIM,0:L(LN)%inne),L(LN)%Qs(DIM,0:mn2Dh))
        L(LN)%U(:,0)= ZERO ; L(LN)%MAN(0) = ZERO ; L(LN)%INW = 0
        L(LN)%ZK(0) = L(LN)%Z(L(LN)%ke) ; L(LN)%ZKm(0) = L(LN)%Z(L(LN)%ke)
        L(LN)%ETAn(0) = L(LN)%Z(L(LN)%ke) ; L(LN)%F(0) = ZERO
        L(LN)%DX(:,0) = ZERO ; L(LN)%Q(:,0) = ZERO ; L(LN)%RX(:,0) = ZERO
        L(LN)%Dn(:,0) = ZERO;L(LN)%Qs(:,0) = ZERO ; L(LN)%Qn(:,0) = ZERO
        L(LN)%D(:,0) = ZERO ; L(LN)%H(:,0) = -L(LN)%Z(L(LN)%ke)
        L(LN)%ETA(0) = L(LN)%Z(L(LN)%ke) ; L(LN)%IND = 0 ; L(LN)%MND = 0
        !Maximum Value arrays
        allocate(L(LN)%ETAmax(0:mn2Dh)) ; L(LN)%ETAmax(0) = ZERO
        allocate(L(LN)%Umax(DIM,0:mn2Dh)) ; L(LN)%Umax(:,0) = ZERO
    endsubroutine alloc_F
!%%%%%%%%%%%%%%%%%%%%% END OF MODULE: ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module arrays    
