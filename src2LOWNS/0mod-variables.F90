!%%%%%%%%%%%%%%%%%%%%% MODULE: VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> Lists the all the variables used in the program
module variables
    use type_list, only: bcond, ranke, boundary, inputcond

    type(bcond) :: openbound
    character*3 :: order
    character*10 :: fmt,fmt1
    real :: refval
    
    ! From wave module
    type(inputcond) :: inflow
    character*255 :: inflow_file
    character*255 :: output_file
    integer :: locnum    
    real*8 :: h_in
    real*8 :: GDR = 1.0d0 ! Ground Deformation Rate
    
    !For 2CLOWNS
    integer :: nstart
    real*8  :: Limit_3D = 1d-8
    
    ! 2DH Variables
    integer :: LNUM = 0                       !Number of layers
    real*8  :: tide_level                     !Tide adjustment   
    character*255 :: inputxyz2D
    character*255 :: inputf2D
    character*255 :: input2D 
    character*255 :: title2D
    
    !Parameters
    real*8,parameter :: ZERO = 0.0d0, HALF = 0.5d0, ONE = 1.0d0, TWO = 2.0d0, &
                        THREE = 3d0, THREEQ = 0.75d0, PI_8 = 4d0*atan(1d0),   &                           
                        TWFTH = 1d0/12d0, FTW7 = 4d0/27d0, THIRD = 1d0/3d0,   &  
                        TWOTHIRD = 2d0/3d0, SIXTH = 1d0/6d0, EIGHTH = 1d0/8d0,&
                        FFTH = 1d0/15d0, TW4 = 1d0/24d0, SQRSIX = sqrt(6d0),  &
                        TW7 = 27d0, FIVE = 5d0,                               &
                        ONEHALF = 1.5d0, SVNTHIRD = 7d0/3d0, V_lim = 20d0,    &
                        TINY = 1d-15, VSMALL = 1d-12, SMALL = 1d-9,           &
                        LITTLE = 1d-3, Leps = 1d0 - 1d-5, F_lim = 1d0 + 1d-5, &
                        DT_lim = 1d-6, TEN = 10d0, TENTH = 0.1d0,             &
                        A = 0.4d0, B = 1d0/15d0, karman = 0.41d0,             &
                        COMFLOW = 0.35d0, SUBFLOW = 0.91d0, HCIRC = 180d0,    &
                        R = 1d0/6378137d0, PSI = 14.58423d-5 !2*psi rad/s
                        
    character(len=255) :: INPUT   ! 入力メインファイル名
    character(len=255) :: INPUTXYZ ! 入力XYZファイル名
    character(len=255) :: INPUTKE ! 入力rdataファイル名

    character(len=255) :: INPUTF='non'  ! 入力fdataファイル名
    character(len=255) :: INPUTf0='non' ! 入力fdataファイル名(地形)
    character(len=255) :: INPUTSHAPE ! 入力odataファイル名
    character(len=255) :: DRIINIT ! 漂流物情報設定ファイル名
    character(len=255) :: INPUTtc      

    character(len=255) :: title,title_temp ! 出力ファイル名
    
    character(len=255) :: TESTA ! 出力ファイル名(t?.???より前の部分）
    character(len=8)   :: idate  !時刻取得用引数
    character(len=255) :: num

    real*8 :: cl_set  !Desired setting of courant number to automatically get dt
    real*8 :: courant !クーラン数
    real*8 :: hasoku_c !波速クーラン数    
    real*8 :: velmax !流速変動の総和
    integer:: ivmax  !　クーラン数に関する情報受け渡し
    integer:: icon,icong ! エラーコード
    integer,parameter::iob=-1 !物体セルのnf値
    integer,parameter::ifu=1 !流体セルのnf値
    integer::n=0 !現在の計算ステップ数
    integer::nkai !計算の繰り返し回数
    integer::mwr ! 情報出力間隔（mwrステップごと）
    integer::dwr ! local data output frequency
    real*8 ::twl ! local data output time
    integer::fg_xyz ! =0 -> .xyz ; =1 -> .xyzn
    integer::ko_tmp ! tempデータを出力する回数（ko_t=ko_tmpのとき，時間付きデータを出力）
    integer::ko_t  ! tempデータを出力した回数
    integer::ko_out  ! sfdataを maindataやrdataよりも細かく出力する場合。（２D解析の場合無視される）
    integer::is,js,ks,ie,je,ke
    integer::inns,inne,inne2,inn2d !セル番号開始と終了,境界セルを含めた終了
    integer::ndri=0 !　漂流物数(初期値0)
    integer::nx1,nx2
    integer:: NTHR = 1 !Number of OMP threads
    integer:: IM = 1   !Model type: = 1 - Hybrid 2D-3D, = 2 - 2DH only, = 3 - 3D only
    

    integer::iplmax
    integer::npln=6  !総平面数　6+sum(ipl(1:ndri))
    integer::ndri_face  ! 平面通し番号
    integer::mm0=0,mm0mx=10000  ! 地形を含むセル数, 地形を含むセル数の最大値
    integer::mm1=0,mm1mx=10000  ! 漂流物を部分的に含むセル数, 漂流物を部分的に含むセル数の最大値
    integer::mm1f=0,mm1fmx=10000 !漂流物に完全に含まれるセルの総数
    integer::mms,mmsmx=10000 ! 水面を含むセル数,水面を含むセル数の最大値
    integer::mms2 ! 水面数（一セルに複数水面を含む。）
    integer::ncpl01max=60 !25!vertex number in a cell ex. ncpl*(0,1) ! 
    integer::ncpl00max=60 !plain number in a cell  ex. ncpl*(0,0) 
    integer::ncpl02max=80 !25!vertex number of each face in a cell ex ncpl(0,2) !mcfvtmx
    integer::minfopath=300 !path数の上限
    integer::ic1_0,ic2_0,ic3_0
    integer::ipl0=6 ! セルの面数
    integer::ncpls1max=50 !pp_s　第二引数の最大値
    integer::ncpls0max=50 !ncpls　第二引数の最大値
    integer::ncpls2max=100 !ncpls　第三引数の最大値
    !integer::ub_surfp_2=30 !surfp　第二引数の最大値

    integer::mmdmx=30000    !　セグメント面の総数の上限
    integer::mmd2mx=30000    !　二次元可視化時セグメント面の総数の上限

    integer::ifln  ! 内部セル+水面セルの総数
    integer::isfn  ! 水面セルの総数
    integer::isfn2  ! 水面セル＋（f<1の内部セル,f>0の空気セル）の総数

    integer :: TOME = 1 !Temporal order of momentum explicit integration 
    
    real*8::t=ZERO,t0 ! 時間,読み込んだデータの時間
    real*8::dt=0.0d0,dtmax,dt0=0.0d0 ! 時間刻み、時間刻み最大値,読み込んだデータの時間刻み
    real*8::dtm !時間刻み，読み込み時間刻み,Previous timestep
    real*8::g  = 9.80665d0 ! 重力加速度 Standard gravity
    real*8::nu = 1.0d-6 ! 動粘性係数 Viscosity of water at 20C
    real*8::d0 ! 連続式打切り誤差
    real*8::eps0 ! 行列計算打切り誤差

    real*8::t_end ! 計算終了時刻
    real*8::tw0  ! 初回データ出力時刻
    real*8::tw   ! 次回データ出力時刻
    real*8::dtw  ! データ出力時間間隔
    real*8::g_read   ! 読み込みデータの重力加速度
    real*8::nu_read   ! 読み込みデータの動粘性係数
    real*8::PI = 4d0*atan(1d0) !円周率
    real*8::sdum= -1d4  ! 水面がない場合のダミー値

    integer :: ike !乱流量の読み込み調整
    integer :: n_ak=1 !繰り返し回数
    integer :: ieee !
    real*8::akrate,erate,drate !乱流量の壁面低減率
    real*8::ratez=ZERO,ratex=ZERO,ratey=ZERO !乱流量の調整
    real*8::aratez=1.0d0,aratex=1.0d0,aratey=1.0d0 !乱流量の調整
    real*8::akrateS !乱流量の水面低減率
    real*8,allocatable,dimension(:)::ak,ak1 !乱流エネ，次ステップ乱流エネ
    real*8,allocatable,dimension(:)::e,e1 !乱流エネ散逸，次ステップ乱流エネ散逸
    real*8,allocatable,dimension(:)::pro !生成項
    real*8,allocatable,dimension(:)::d !渦動粘性係数
    real*8,parameter:: c_smag = 0.2d0
    ! k - eps coefficients                        ! Realizable coefficients
    real*8,parameter::ce1 = 1.44d0,ce2i = 1.92d0, ce2r = 1.9d0, c1r = 0.43d0, A0 = 4.04 
    real*8,parameter::cm = 0.09d0, Sk = 1.0d0, Se = 1.3d0, Ser = 1.2d0 
    real*8 :: mu ! 粘性係数
    real*8 :: ks_r(0:10) !Equivalent sand roughness (m) !<0 = no wall function (no gradient), 0 = smooth wall, >0 rough wall
    integer :: pro_lim !Limiter type for production: = 0 no limiter; = 1 Kato-Launder Mod; = 2 Menter limiter
    integer :: vt_op   !Option for varying nu_t for low RE
    type(ranke) :: r_lim !Limit for k and e
    real*8,parameter :: akmin = 1.92d-15, emin = 1.0d-15, D_rlim = 1d8  
    real*8 :: akini = 1.92d-15, eini = 1.0d-15
    real*8 :: TI_b  = 0.01d0, nut_b = 1d0
    real*8 :: honma = 0.35d0 !本間式
    real*8 :: svmax !　クーラン数に関する情報受け渡し
    ! Courant:  lowest limit   small limit     medium limit
    real*8 :: cr_ll = 0.1d0, cr_sl = 0.2d0, cr_ml = 0.4d0
    !          big limit        upper limit
    real*8 :: cr_bl = 0.6d0, cr_ul = 0.9d0
    
    real*8 :: Uweight = 4d0    !Weight given to Courant number using fluid velocity in 
                               !comparison to the wave speed and turbulent time scale  
    real*8 :: Dweight = 1.1d0  !Weight given to Courant number using diffusion in 
                               !comparison to the wave speed and turbulent time scale  
    real*8::prt=1.0d0/1.6d0
    real*8::prc=1.0d0/1.2d0
    real*8::rhos,dakuin
    real*8::c0,c00,c1,c10
    integer :: ifk0,mfk1
    real*8::rho0=1000.d0
    
    real*8::myu=HALF,myu1=0.3d0  ! myu:最大静止摩擦係数，myu1:動摩擦係数
    real*8::cr !=0.1d  ! cr:反発係数      
    integer :: mm1k, mm11k

    !Matrix ones
    integer::nl
    integer::mm2
    real*8 ::eps   ! 収束計算の打切り誤差(結果)
    integer::iter   ! 収束計算の繰り返し回数
    integer::mnoe,m_alloc=0
    integer::ik
    real*8::DMAX
    integer::igs,ige,ricon
    real*8 :: mp ! a multiplier of the matrix for normalisation purposes 
    
end module variables 
