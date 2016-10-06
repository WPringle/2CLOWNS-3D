!%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-normal-PlainAdjust.F90 %%%%%%%%%%%%%%%%%%%%%%
! Calculates the volume and planes of the segment in the cell 
subroutine PlainAdjust(nn,vv_1,surf,iz1,iz2,ncpl2,axyz_sd,iz3,pp,iiz1,iiz2,ncpl3)
!　ポリゴンで囲まれた形状A　を　面S　により切断し、新たな形状　Bを作成する。
    !nn,int:対象セル番号(in)
    !vv_1,dbl:形状Bの体積(out)
    !surf,dbl(0:3):面Sの方程式a,b,c,d (in)   
    !pp,dbl(:,0:3):形状Aの頂点 兼形状Bの頂点(in out)
    !ncpl2,int(:,:):形状Aの形状情報(in)
    !ncpl3,int(:,:):形状Bの形状情報(out)
    !axyz_sd,dbl(:):形状B各面の開口率(out) %底面の場合

    ! ipp           ：頂点番号（総数）
    ! vt_info(ipp)  : =2 水面に含まれる頂点，=1 水面下の頂点，=0 水面上の頂点
    ! men_info(icpl)：=0　すべての頂点が水面下,　=-1　すべての頂点が水面上，=1 水面を含む
    ! info(ipp)     : =2 水面に含まれる頂点，=1 水面下の頂点，=0 水面上の頂点
#ifdef CK_SEGMENT
    use param2
    use interface_list, only: dista, crossV
#else
    use interface_list, only: dista, crossV
    use arrays, only: cpl, dx, in, jn, kn
    use variables,only: pi, ZERO, SIXTH, VSMALL, ncpl02max, ncpl00max
    use omp_lib
#endif   
    implicit none
!------------------------- In/Out Variables ------------------------------------
    integer,intent(in)               :: nn, iz1, iz2, iz3, iiz1, iiz2
    real*8,dimension(0:3),intent(in) :: surf
    real*8,intent(inout)             :: pp(0:iz3,0:2), axyz_sd(6)  
    real*8,intent(out)               :: vv_1
    integer,dimension(0:iiz1,-2:iiz2),intent(out) :: ncpl3
    integer,dimension(0:iz1,-2:iz2),intent(in) :: ncpl2   
!------------------------- Temporary variables ---------------------------------
#ifdef BEDMOVE
    integer                      :: ipath3,kk,mm
    integer,dimension(0:50,-2:3) :: path3
#endif
    integer :: iis,iie,irmin,icpl2,ia,icc,icpl3,ipath,iap,ib,ib2,ib1,nnp,ia1,ic
    integer :: spath,ncmax,ncmin,icc1,ii,iam,ipath2,ic2,ic22,icc3,ib4 !,icpath
    integer :: ipp,myhr,i_1,ncmin1,ncmax1,ic3,ic1,icca,ic1a !,loop,i1,ibb,nnp1
    integer :: ic1b,onevt,icpl20,impt !,iib,idd,ico,nnpall,loop3,loop4,ppsize
    integer :: ic4,ic7,ic8 !,ibw,ia2,ib3,ib3p,ic5,ic6,ic9,iiaa,iibb,iicc
    integer :: ii_2,ii_1,i_2,i_3,i_4,i_5,i_6,ibp,ib5,i_4p !,j1,j2,j3,i_7,i_8
    integer :: icc3v(10), nmpt(10) !,igpath(11:13), icv3(20,20)
    integer,dimension(20) :: icv6,icv7,icv8
    integer,dimension(0:300) :: icv,icv1,icv2,info,men_info,vt_info !,icvm,icvn,icvnt
    integer,dimension(0:50,-2:3) :: path,path2 
    integer,dimension(0:50,-2:4) :: paths
    integer,dimension(0:50,0:50) :: ivv
    integer,allocatable,dimension(:)   :: imin,icv4,icv5 !,ictemp
    integer,allocatable,dimension(:,:) :: route,routem1,routem2 !,route2,route2m1,route2m2
    real*8 :: dd1,dd2,r1,r2,ttemp !,dd3,bm1,bm2,aa,bb,cc,round,pp_round
    real*8 :: dda(ncpl02max+5), tmpVV(10,0:2), rr(10), vv(0:30) !, gpathV(11:13,0:2)  
    real*8,dimension(0:2) :: tmpV,tmpV1,tmpV2,ppc !,tmpV3,tmpV4
    real*8,dimension(0:3) :: surf1,cpltemp    
!    real*8,allocatable,dimension(:) :: stemp  
!    real*8,allocatable,dimension(:,:) :: pp_temp
#ifdef BUGWIN
!    real*8,allocatable,dimension(:) :: s,leng
!    real*8 :: len
!===============================================================================
    info=0 ; myhr = 0
!$  myhr=omp_get_thread_num()
#else
    myhr=0
#endif
    !cpl2とncpl2には手を加えない．
    ncpl3=ncpl2 ! 形状Aの情報をBに写す。
    icpl2=ncpl2(0,0)  !Aの面の数
    ncpl3(0,0)=icpl2+1 !Bの面の数(水面考慮)：仮定
    ipp=ncpl2(0,1)    !頂点の数(水面考慮以前）
    icpl3=0;icca=0;ic1a=0;ic1b=0;ppc=0
 
! 1. 形状Aと面Sとの関係を確認し、分割の必要性を判断
! 1-1. 形状Aの頂点と面S（以下水面とする）との関係を登録、分類
    do i_1=1,ipp    
    ! vt_info(頂点番号)=2 ;   水面に含まれる頂点
    ! vt_info(頂点番号)=1 ;   水面より下の頂点
    ! vt_info(頂点番号)=0 ;   水面より上の頂点
    !水面に含まれる頂点の数 icca++
    !水面下の頂点の数 ic1a++
        dda(i_1) = sum( surf(0:2) * pp(i_1,0:2) ) + surf(3)  
        if (abs(dda(i_1)).lt.1d-14) then !水面に含まれる頂点の数 icca++
            vt_info(i_1)=2
            icca=icca+1     
            onevt=i_1
            ppc=ppc+pp(i_1,0:2)
        elseif (dda(i_1).lt.ZERO) then !水面下の頂点の数 ic1a++
            vt_info(i_1)=1
            ic1a=ic1a+1 
        else
            vt_info(i_1)=0 !水面上の頂点
        endif
    enddo
! 1-2.水面に含まれる頂点の数 icca 、水面下の頂点の数 ic1a　により形状Aと水面関係を判定
    if(icca+ic1a.eq.ipp) then !すべて水面or水面下の場合 
        ncpl3(icpl2+1,-1)=2   !icpl2+1を水面(=2)として登録
        if(icca.eq.0) then          !水面に含まれる頂点がない場合
            ncpl3(icpl2+1,0)=0            !水面なし
            goto 346
        else if(icca.lt.3) then     !水面に含まれる頂点で面が形成できない場合
            ncpl3(icpl2+1,0)=1          !頂点が１点の面として登録
            ncpl3(icpl2+1,1)=onevt      !１点の頂点番号を登録
            goto 346
        else                    ! icca >=3 && icca>0 ; 水面に一致する面もあり、水面下にも頂点あり 
            icpl20=icpl2
            do i_1=1,icpl2      ! 形状Aの各面i_1について 
                icc=0
                do ia=1,ncpl2(i_1,0)    ! 面i_1の頂点iaが水面に含まれる場合、icc++ 
                    if(vt_info(ncpl2(i_1,ia)).eq.2) icc=icc+1   
                enddo
                if(icc.eq.ncpl2(i_1,0)) then   ! 面i_1のすべての頂点が水面に含まれる場合、それも水面として登録 
                    icpl20=icpl20+1     ! 面の数を増やす
                    ncpl3(icpl20,0:ncpl2(i_1,0))=ncpl2(i_1,0:ncpl2(i_1,0))  !面i_1の情報をコピー
                    ncpl3(icpl20,-1)=2   !icpl2+1を水面として登録    
                    continue
                else if(icc.ge.3) then   ! 面i_1の三つの頂点が水面に含まれる場合、残りが含まれない場合・おかしい 
                    ncpl3(icpl2+1,0)=1
                    ncpl3(icpl2+1,1)=onevt
                    !cpl22(icpl2+1,:)=surf(0:3)  !水面の登録
                    ncpl3(icpl2+1,-1)=2   !icpl2+1を水面として登録                
                    goto 346
                else 
                endif
            enddo        
            if(icpl20.gt.icpl2) then ! 面が増えていれば、形状Bの面の数を新しく設定
                ncpl3(0,0)=icpl20
                goto 345
            else ! 面が増えてなかった場合　一直線上にある可能性大（おかしい）
                ncpl3(icpl2+1,0)=1
                ncpl3(icpl2+1,1)=onevt    
                goto 346     
            endif 

        endif    
!すべて水面下の頂点が０の場合
    else if(ic1a.eq.0) then !水面下の頂点数が0 すなわち、水面下に形状なし
        !Switch 水面下の頂点 and 水面に含まれる頂点
        vv_1=ZERO; return
    endif
    
! 2.形状Aを構成する各ポリゴンを水面が横切るか調べる。（横切る場合には、水面を含む新しいポリゴンを作成）    
ipath=0 ! 水面を横切るポリゴンi_1のパス数を初期化
every_men: do i_1=1,icpl2    ! 2-1.ポリゴンi_1と水面の関係を調べ、i_1を切断するか調べる。
876     icc=0;ic1=0
        do ia=1,ncpl2(i_1,0)    !形状Aの面i_1の頂点iaが・・・
            ic1b=ic1b+1
            if(vt_info(ncpl2(i_1,ia)).eq.2) then ! 水面に含まれる場合　icc++
                icc=icc+1
            else if(vt_info(ncpl2(i_1,ia)).eq.1) then ! 水面下の場合   ic1++
                ic1=ic1+1
            else
            endif
        enddo
        
        if(ic1.eq.ncpl2(i_1,0)) then !面i_1の頂点がすべて水面下ならば・・・
            men_info(i_1)=0;cycle             !切断の必要はないが、面i_1は形状Bには含まれる。
        else if(ic1.eq.0) then       !面i_1の頂点に水面下がないならば・・・
            men_info(i_1)=-1;cycle            !切断の必要はなく、面i_1は形状Bにも含まれない。
        else if(icc.eq.1.and.icc+ic1.eq.ncpl2(i_1,0)) then !面i_1の頂点の一つが水面に含まれ、その他が水面下ならば・・・
            men_info(i_1)=0;cycle                                           !切断の必要はないが、面i_1は形状Bには含まれる。
        else if(icc.eq.ncpl2(i_1,0)) then !面i_1の頂点がすべて水面に含まれるとき・・・
            if(ncpl2(i_1,-1).eq.0) then             !面i_1の面種がセル境界面(=0)のとき・・・
                if(dot_product(surf(0:2),cpl(ncpl2(i_1,-2),0:2)).lt.ZERO) then     !面i_1の法線と水面の法線が逆向きの場合・・・
                    write(6,*) 'mondai ari'
                    stop
                else if(dot_product(surf(0:2),cpl(ncpl2(i_1,-2),0:2)).gt.ZERO) then     !面i_1の法線と水面の法線が同じ向きの場合・・・
                    ncpl3(0,0)=icpl2        ! 形状Bの面の数を形状Aの面の数にする。(面は増えない。)
                    ncpl3(i_1,-1)=2         ! 面i_1の面種を水面(=2)に変更
                    goto 345
                else
                    write(6,*) 'mondai ari'
                    stop
                endif
            else 
                continue
                men_info(i_1)=-1;cycle        
            endif
        else ! それ以外の場合(水面が横切る可能性あり)
            men_info(i_1)=1
        endif

! 面i_1を切断する。　
 
! まず、面_1について、水面を横切る辺ABをpathに登録　水面との交点Sの頂点番号をpath(3)に登録。
        ic2=0
    every_vt: do ia=1,ncpl2(i_1,0)
!線分ABについて・・・・
        iap=ia+1; if(iap.gt.ncpl2(i_1,0)) iap=1  
        
        !すでに登録されているか調べる。
        icc=0
        do ib=1,ipath
            if(path(ib,1)==ncpl2(i_1,ia).and.path(ib,2)==ncpl2(i_1,iap)) icc=1
            if(path(ib,2)==ncpl2(i_1,ia).and.path(ib,1)==ncpl2(i_1,iap)) icc=1
        enddo
        if(icc.eq.1)  cycle 
        !まだ登録されていない場合
        
! (Aが水面下 かつ Bが水面上) または　(Aが水面上 かつ Bが水面下)　の場合、その間の水面で切断する。　
        if ((vt_info(ncpl2(i_1,ia)).eq.1.and.vt_info(ncpl2(i_1,iap)).eq.0).or.&
           (vt_info(ncpl2(i_1,ia)).eq.0.and.vt_info(ncpl2(i_1,iap)).eq.1)) then
            ipath=ipath+1       !   線分番号を設定
            ic2=ic2+1
            path(ipath,1)=ncpl2(i_1,ia);path(ipath,2)=ncpl2(i_1,iap)  !線分ABをpathに登録  
            tmpV = ( pp(ncpl2(i_1,ia),:) * abs(dda(ncpl2(i_1,iap))) + pp(ncpl2(i_1,iap),:)           &
                   * abs(dda(ncpl2(i_1,ia))) ) / (abs(dda(ncpl2(i_1,iap)) ) + abs(dda(ncpl2(i_1,ia))))  !線分ABと水面の交点S(の座標)を求める。   
            dd1=dista(tmpV-pp(ncpl2(i_1,ia),:))  !AとSの距離
            dd2=dista(tmpV-pp(ncpl2(i_1,iap),:)) !BとSの距離
            ic3=0
            do ic=1,ipp !交点Sがすでに頂点として登録されているか検索。
                if(ic.eq.ncpl2(i_1,ia)) cycle
                if(ic.eq.ncpl2(i_1,iap)) cycle    
                if(dista(tmpV-pp(ic,:)).lt.1.0d-9) ic3=ic ! 交点Sが既存の頂点に一致する場合、既存の頂点番号をic3に記録。
            enddo
            if(ic3.eq.0) then !交点Sが登録されていない(ic3=0)ので、交点を水面として登録
                ipp=ipp+1 !セルの登録頂点総数を増やす。交点Sの頂点番号をippとする。
                path(ipath,3)=ipp ! 線分AB内に交点S(頂点番号ipp)が存在することを登録                   
                pp(ipp,:)=tmpV !交点S(頂点番号ipp)の座標を登録。                                
                vt_info(ipp)=2 !交点S(頂点番号ipp)が水面に含まれることを登録。 
            else !水面との交点が既存の点（頂点番号ic3）の場合
                path(ipath,3)=ic3 ! 線分AB内に交点S(頂点番号ic3)が存在することを登録 
                vt_info(ic3)=2 !交点S(頂点番号ic3)が水面に含まれることを登録。      
            endif
        endif

        enddo every_vt

enddo every_men

!CC　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　この時点で、ipathは水面と交点を持つ形状Aの辺の総数。
! 3.各ポリゴンから水上部分を取り除いたポリゴンを作成する。また、水面となるポリゴンを構成するパスを登録。
ncpl3(0,1)=ipp !水面を考慮したセル全体の頂点数
spath=0         !水面を構成するバス数を初期化。

every_men2:do i_1=1,icpl2   
    icpl3=0 !平面i_1の水中部分の頂点数icpl3を初期化
#ifdef CK_SEGMENT
     write(6,1223) i_1,men_info(i_1)
1223 format(1h ,'memNo=',i3,',men_info=',i3)
#endif
    if(men_info(i_1).ne.1) then !ポリゴンi_1が水面を横切らない場合の処理   
        if(men_info(i_1).eq.-1) then ! 面_1がすべて水面上の場合、
            ncpl3(i_1,0)=icpl3 !            面の頂点数を0として、以降考慮しない。->次の面に 
        endif                        ! 面_1がすべて水面下の場合、そのまま、形状Bに利用する。
        cycle                     
    endif

!CC　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　水面を挟む（水面上の頂点、水面下の頂点をともに含む）場合、
! 3-1. 各ポリゴン頂点から、水面上のものを除き、水面との交点を加える。
    every_vt2: do ia=1,ncpl2(i_1,0)  !頂点(i_1,ia)を頂点Aとする。
        iap=ia+1; if(iap.gt.ncpl2(i_1,0)) iap=1   !頂点(i_1,iap)を頂点APとする。
        
        if(vt_info(ncpl2(i_1,ia)).eq.2) then !頂点Aが水面内の場合・・・・・
            iam=ia-1;if(iam.eq.0) iam=ncpl2(i_1,0)   !頂点(i_1,iam)を頂点AMとする。
            if(vt_info(ncpl2(i_1,iap)).eq.0.and.vt_info(ncpl2(i_1,iam)).eq.0) then ! 頂点AM,AP が水面上の場合・・・・
                                                                                        ! ncpl3に加えない。
            elseif (ncpl2(i_1,-2).eq.2.and.ncpl2(i_1,-1).eq.0.and.                    &
                    vt_info(ncpl2(i_1,iap)).eq.1.and.vt_info(ncpl2(i_1,iam)).eq.1) then   ! ポリゴンi_1がセル境界面　かつ　Z＋の面の場合で、頂点AM,APが水面下の場合    
                icpl3=icpl3+1; ncpl3(i_1,icpl3)=ncpl2(i_1,ia)   !頂点番号を登録．
                info(ncpl2(i_1,ia))=1            !水面下頂点として登録
#ifdef CK_SEGMENT
                write(6,1224) ncpl2(i_1,ia),info(ncpl2(i_1,ia))
1224            format(1h ,'----vtNo=',i3,',vtInfo=',i3)
#endif
            else
                icpl3=icpl3+1; ncpl3(i_1,icpl3)=ncpl2(i_1,ia)   !頂点番号を登録．
                info(ncpl2(i_1,ia))=2            !水面内頂点として登録
#ifdef CK_SEGMENT
                write(6,1224) ncpl2(i_1,ia),info(ncpl2(i_1,ia))
#endif
            endif
        elseif (vt_info(ncpl2(i_1,ia)).eq.1.or.vt_info(ncpl2(i_1,ia)).eq.0) then   !頂点Aが水面下または水面上の場合．
            icpl3=icpl3+1
            ncpl3(i_1,icpl3)=ncpl2(i_1,ia)   !頂点番号を登録．
            info(ncpl2(i_1,ia))=vt_info(ncpl2(i_1,ia))   !頂点情報をコピー
#ifdef CK_SEGMENT
            write(6,1224) ncpl2(i_1,ia),info(ncpl2(i_1,ia))
#endif  
        endif
        
        icc=0
        do ib=1,ipath   !頂点Aと次の頂点APとの間に水面との交点Sがあるか調べる　(ある場合には、pathに登録されている。)   
            if(path(ib,1)==ncpl2(i_1,ia).and.path(ib,2)==ncpl2(i_1,iap)) icc=ib
            if(path(ib,2)==ncpl2(i_1,ia).and.path(ib,1)==ncpl2(i_1,iap)) icc=ib
        enddo
        if(icc.ne.0) then !ある場合
            icpl3=icpl3+1 ; ncpl3(i_1,icpl3)=path(icc,3) !交点Sの頂点番号を登録．
            info(path(icc,3))=2          !水面内頂点として登録
#ifdef CK_SEGMENT
            write(6,1225) path(icc,3),info(path(icc,3)),path(icc,1),path(icc,2)
1225        format(1h ,'----vtNo=',i3,',vtInfo=',i3,' btw=',i3,' and',i3)
#endif
        endif
    enddo every_vt2

    ncpl3(i_1,0)=icpl3  !ポリゴンi_1の水中部分の頂点数を登録
    
    if(ipath.eq.0.and.men_info(i_1).ne.1) goto 214    !平面!が水面に絡まない場合はここまで．

    ipath2=0  !A,APの少なくともどちらかが水面下となるパス（水中部分の候補）を見つけて、path2に登録。総数ipath2
    do ia=1,ncpl3(i_1,0)                                     !ncpl3(i_1,ia)を頂点Aとする。
        iap=ia+1; if(iap.gt.ncpl3(i_1,0)) iap=1             !ncpl3(i_1,iap)を頂点APとする。
        if(info(ncpl3(i_1,ia)).eq.1.or.info(ncpl3(i_1,iap)).eq.1) then
            ipath2=ipath2+1
            path2(ipath2,1)=ncpl3(i_1,ia)   ! 頂点Aを登録
            path2(ipath2,2)=ncpl3(i_1,iap)  ! 頂点APを登録
#ifdef CK_SEGMENT
            write(6,1226) ncpl3(i_1,ia),ncpl3(i_1,iap)
1226        format(1h ,'---------どちらか一端が水面下となる線分=',i3,i3)
#endif
        endif
#ifdef CK_SEGMENT
        !       write(6,12261) ncpl3(i_1,ia),ncpl3(i_1,iap)
        !12261 format(1h ,'---------両端が水面となる線分=',i3,i3)
#endif
        !endif
    enddo

#ifdef BEDMOVE
    if (ground.eq.1) then
        if(ncpl2(i_1,-2).eq.1.or.ncpl2(i_1,-2).eq.3.or.ncpl2(i_1,-2).eq.6) then
            ipath3=0
            do ia=1,ncpl3(i_1,0)    !頂点iaについて
                iap=ia+1; if(iap.gt.ncpl3(i_1,0)) iap=1
                if(info(ncpl3(i_1,ia)).eq.0.or.info(ncpl3(i_1,iap)).eq.0) then
                    ipath3=ipath3+1
                    path3(ipath3,1)=ncpl3(i_1,ia)
                    path3(ipath3,2)=ncpl3(i_1,iap)
                    path3(ipath3,3)=1
                endif
            enddo

            if(ipath3.gt.1) then

            !水面から水面上に向かうパスを検索し、個数をic2とする。ic2は水面の分割数となる。
            ic2=0
            do ia=1,ipath3
                if(info(path3(ia,1)).ne.2) cycle
                ic2=ic2+1
                ivv(ic2,1:2)=path3(ia,1:2)
                path3(ia,3)=0
            enddo 
            if(ic2.gt.2) then
                continue
            else
                ic=1
                ic22=2
                ib1=1
                do while(ic22>0)
                    continue
                    if(path3(ib1,3).ne.0) then
                        if(path3(ib1,1).eq.ivv(ic,ic22)) then
                            ic22=ic22+1
                            ivv(ic,ic22)=path3(ib1,2)
                            if(info(ivv(ic,ic22)).eq.2.and.ic22.gt.2) then
                                ivv(ic,0)=ic22
                                exit
                            endif  
                        endif
                    endif
                    ib1=ib1+1
                    if(ib1.gt.ipath3) ib1=ib1-ipath3
                enddo
                ia=ncpl2(i_1,-2)    

                select case(ia)
                case(1);ib=3
                case(3);ib=2 
                case(6);ib=1             
                end select
                nd3nn(ib)=ivv(1,0)
                do kk=1,nd3nn(ib)
                    dp3nn(ib,kk,0:2)=pp(ivv(1,kk),0:2)
                enddo    
                continue            
            endif            
            endif 
        endif   

    endif
#endif

if(ipath2.lt.ncpl3(i_1,0)) then

! 3-2. 水面内頂点から始まり再び水面に到達する頂点の並びをivv(ic2)とする。該当する並びがポリゴンi_1に複数ある場合ic2>1. 
        ic2=0; path2(1:ipath2,3)=1 ! ic2,path2(:,3)を初期化。
        do ia=1,ipath2 !水面から水中に向かうパスを検索しivvの先頭（ivv(ic2,1:2)）に登録
            if(info(path2(ia,1)).ne.2) cycle ! 頂点Aが水面内の場合のみ検討。（path2はどちらかが水面下なので、頂点Aが水面内であればBは水面下となる。）
            ic2=ic2+1 ; ivv(ic2,1:2)=path2(ia,1:2) !path2(ia)をivvに登録。
            path2(ia,3)=0 !ivvに登録されていることを記録。
#ifdef CK_SEGMENT
            write(6,1227) path2(ia,1:2)
1227        format(1h ,'---------水面から水中に向かうパスを検索=',i3,i3)
#endif
        enddo       
        do ic=1,ic2  ! ivv(ic,1:2)に接続し、水面まで到達する頂点の並びを検索する。 
            ic22=2   ! ivv(ic)の頂点数、最初は２
            ib1=1    ! path2(1)から検索開始
            do while(ic22>0)
                continue
                if(path2(ib1,3).ne.0) then ! path2(ib1)がivvにまだ選択されていないとき，
                    if(path2(ib1,1).eq.ivv(ic,ic22)) then       ! path2(ib1)の始点が，ivvの終点(ivv(ic,2))と等しい場合
                        ic22=ic22+1;ivv(ic,ic22)=path2(ib1,2)   ! path2(ib1)の終点をivv(ic)の並びに追加する．
                        if(info(ivv(ic,ic22)).eq.2.and.ic22.gt.2) then ! 追加したivv(ic)が，水面内頂点であるとき(終了)。
                            ivv(ic,0)=ic22 ! ivv(ic)の並び数をic22として、新しいivvに移る。
                            exit
                        endif  
                    endif
                endif
                ib1=ib1+1 !次のpath2を検討する。
                if(ib1.gt.ipath2) ib1=ib1-ipath2
            enddo
        enddo

349     if(ic2.eq.1) then  ! 水面から水中に向かうpath2が１つの場合(ic2=1)，ivv(1)の並びは分割されたポリゴンの並びに等しい。
            ncpl3(i_1,0:ivv(1,0))=ivv(1,0:ivv(1,0)) 
#ifdef CK_SEGMENT
            write(6,1228) ncpl3(i_1,0:ncpl3(i_1,0))
1228        format(1h ,'-------------水中部の並び=',20(',',i3))
#endif
            goto 214
        else if(ic2.eq.0) then ! ic2のときエラーの可能性。
            do i_4=1,ncpl3(i_1,0)
             write(6,*) ncpl3(i_1,i_4),vt_info(ncpl3(i_1,i_4)),dda(ncpl3(i_1,i_4)),'ic2=',ic2
            enddo
            stop
        endif

!CC　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　 ic2>1の場合、切断後のポリゴンi_1の分割が必要となる場合もあり。
        icc=0         
        do ia1=2,ic2!ivv(ia)の始点から終点に向かうベクトルと，vv(1)の始点から終点に向かうベクトルとの向きを調べる．
            dd1 = sum( (pp(ivv(ia1,ivv(ia1,0)),0:2) - pp(ivv(ia1,1),0:2) ) &
                 * ( pp(ivv(1,ivv(1,0)),0:2) - pp(ivv(1,1),0:2) ) )
            if(dd1.lt.ZERO) icc=icc+1   !逆方向の場合icc=1とする．
        enddo

        if(icc.eq.0) then  !すべて向きが同じ ic2の数だけ面が分割される。
            do ia=0,ivv(1,0) !一つ目の面を面i_1とする。
                ncpl3(i_1,ia)=ivv(1,ia)
            enddo
            do ib=2,ic2 !二つ目以降の面は、面を新たに追加
                icpl2=icpl2+1                       !面icpl2を追加する
                ncpl3(icpl2,-2:-1)=ncpl2(i_1,-2:-1) !面i_1の情報を面icpl2にコピー
                ncpl3(0,0)=icpl2                    !形状Bの面の総数を増加させる。
                do ia=0,ivv(ib,0)                   ! 頂点の数、並びを登録。
                    ncpl3(icpl2,ia)=ivv(ib,ia)
                enddo   
            enddo
        else if(ic2.eq.2.and.icc.eq.1) then !水面下パスが二本で、逆向きの場合。一つのポリゴンで表すことができる。
            ncpl3(i_1,0:ivv(1,0))=ivv(1,0:ivv(1,0))  ! ic2=1の情報を登録。
            if(ncpl3(i_1,ivv(1,0)).eq.ivv(2,1)) then ! ivv(1)の終点と、ivv(2)の始点が同一点の場合。
                ncpl3(i_1,ivv(1,0)+1:ivv(1,0)+ivv(2,0)-1)=ivv(2,2:ivv(2,0))  
                ncpl3(i_1,0)=ivv(1,0)+ivv(2,0)-1
            else! ivv(1)の終点と、ivv(2)の始点が異なる場合。
                ncpl3(i_1,ivv(1,0)+1:ivv(1,0)+ivv(2,0))=ivv(2,1:ivv(2,0))  
                ncpl3(i_1,0)=ivv(1,0)+ivv(2,0)
            endif
            if(ncpl3(i_1,1).eq.ncpl3(i_1,ncpl3(i_1,0))) ncpl3(i_1,0)=ncpl3(i_1,0)-1 !頂点並びの始点と終点が統一点の場合
        else !想定外
            write(6,*) 'souteigai;;ic3,icc=',ic2,icc;stop
        endif
    endif

214 continue
#ifdef BUGWIN
    if(ipath.eq.0) then
        continue
    endif
#endif

!CC　　　　　　     　　　　　　　　　　　　　　　　　　　　　　　　　CCCCCCC 水面の頂点が他の水面に含まれる辺と接触する場合。
!CC　　　　　     　　　　　　　　　　　　　　　　　　　　　　　　　　CCCCCCC 水面の頂点で前後の頂点が水面下の場合は、接触する可能性がある。
    
! 
    impt=0
    do ia=1,ncpl3(i_1,0)   !　接触する可能性がある頂点番号をnmptに格納。総数はimpt
        if(info(ncpl3(i_1,ia)).ne.2) cycle     
        iap=ia+1; if(iap.gt.ncpl3(i_1,0)) iap=1
        iam=ia-1; if(iam.eq.0) iam=ncpl3(i_1,0)
        if(info(ncpl3(i_1,iap)).ne.2.and.info(ncpl3(i_1,iam)).ne.2) then     
            impt=impt+1
            nmpt(impt)=ncpl3(i_1,ia)   
        endif
    enddo

    ncmin1=0;ncmax1=0
! 3-3. 切断後のポリゴンパスのうち、水面内のものをpathsに登録（総数spath）。水面内のpathが他の頂点と接触する場合、その点でパスを分割する。
    do ia=1,ncpl3(i_1,0)  
        if(info(ncpl3(i_1,ia)).ne.2) cycle ! 頂点が水面に含れる場合のみ対象 
        iap=ia+1; if(iap.gt.ncpl3(i_1,0)) iap=1
        if(ncpl3(i_1,ia).eq.ncpl3(i_1,iap)) cycle                         !始点と終点が同じものをのぞく
        if(info(ncpl3(i_1,iap)).eq.2) then   ! A、APが水面内頂点の場合 
            ncmin=min(ncpl3(i_1,ia),ncpl3(i_1,iap))
            ncmax=max(ncpl3(i_1,ia),ncpl3(i_1,iap))
            icc=0
            if(icc.eq.0) then
                dd1=1.0d0
                do ic1=1,impt !辺A-AP内にnmptが含まれるか調べる。
                    dd1=sum( (pp(ncmax,:)-pp(nmpt(impt),:))*(pp(ncmin,:)-pp(nmpt(impt),:)) ) !nmptから辺の両端へのベクトルの
                    dd2=dista(pp(ncmax,:)-pp(nmpt(impt),:))*dista(pp(ncmin,:)-pp(nmpt(impt),:))
                    if(dd2.ne.ZERO) dd1=dd1/dd2
                    continue
                enddo
                ncmin1=0;ncmax1=0
                if((dd1).eq.-1.d0) then     ! 両端が水面内頂点のパスに両側が水面でない頂点が含まれた場合，その点で水面内パスを二つに切断
                    spath=spath+1
                    paths(spath,1)=ncmin;ncmin1=ncmin
                    paths(spath,2)=nmpt(impt)
                    paths(spath,3)=1
                    paths(spath,-2:-1)=ncpl2(i_1,-2:-1)
                    spath=spath+1
                    paths(spath,1)=nmpt(impt)
                    paths(spath,2)=ncmax;ncmax1=ncmax
                    paths(spath,3)=1
                    paths(spath,-2:-1)=ncpl2(i_1,-2:-1)            
                else                      ! 含まれない場合、水面内パスとして登録　　
                    spath=spath+1
                    paths(spath,1)=ncmin
                    paths(spath,2)=ncmax
                    paths(spath,3)=1
                    paths(spath,-2:-1)=ncpl2(i_1,-2:-1)


#ifdef CK_SEGMENT
                    write(6,1229) spath,paths(spath,1:2)
                    !1229   format(1h ,'-------------水面パスNo=',i3,'水面パス始点終点=',i3,i3)
#endif
                    !         endif
                endif   
            endif
        endif
    enddo

    if(impt.eq.0) then
        cycle
    else if(impt.ge.2) then
        write(6,*) 'souteiga impt=',impt
        stop
    endif
    
! 3-4, 水面内パスがnmpt(impt)で切断された場合、ポリゴンi_1もその点で分割する
    if(ncmax1.ne.ncmin1) then
        do ia=1,ncpl3(i_1,0) ! ポリゴンi_1の頂点並びncpl3(i_1,:)において、nmpt(impt)を切断した水面内パスの間に挿入。
            if(ncpl3(i_1,ia).eq.ncmax1.or.ncpl3(i_1,ia).eq.ncmin1) then
                ncpl3(i_1,ia+2:ncpl3(i_1,0)+1)=ncpl3(i_1,ia+1:ncpl3(i_1,0))
                ncpl3(i_1,ia+1)=nmpt(impt);exit
                continue
            endif
        enddo
        ncpl3(i_1,0)=ncpl3(i_1,0)+1 !ポリゴンの頂点数が１増える。
        continue
        ncmin1=0
        do ia=1,ncpl3(i_1,0) ! nmpt(impt)が、ncpl3(i_1,:)の中で何番目と何番目かを調べる。
            if(ncpl3(i_1,ia).eq.nmpt(impt).and.ncmin1.eq.0) then
                ncmin1=ia
            endif
            if(ncpl3(i_1,ia).eq.nmpt(impt).and.ncmin1.ne.0) then
                ncmax1=ia
            endif
        enddo 
        do ia=ncmin1,ncmax1-1 !! nmpt(impt)が、ncpl3(i_1,:)の中ではじめに出てから、次に出る直前までをivv(1)とする。
            ivv(1,ia-ncmin1+1)=ncpl3(i_1,ia)
        enddo 
        ivv(1,0)=ncmax1-1-ncmin1+1
        continue
        ncpl3(i_1,1:ivv(1,0))=ivv(1,1:ivv(1,0))
        ncpl3(i_1,0)=ivv(1,0)        
        !残りでもう一つivv(2)を作成する。
        !       並びを設定する。
        ivv(2,1:ncmin1)=ncpl3(i_1,1:ncmin1) 
        ivv(2,ncmin1+1:ncmin1+ncpl3(i_1,0)-ncmax1)=ncpl3(i_1,ncmax1+1:ncpl3(i_1,0))
        ivv(2,0)=ncmin1+ncpl3(i_1,0)-ncmax1
        continue
        !       ポリゴンを追加する。
        icpl2=icpl2+1;ncpl3(0,0)=icpl2
        ncpl3(icpl2,-2:-1)=ncpl3(i_1,-2:-1)
        ncpl3(icpl2,1:ivv(2,0))=ivv(2,1:ivv(2,0))
        ncpl3(icpl2,0)=ivv(2,0)
    endif
    
enddo every_men2

! 3-5. ポリゴン面が形状Aより増えている場合、増えた面のに水面内パスを追加。ただし、頂点とパスの接触についてcheckしていない。
if(icpl2.gt.ncpl2(0,0)) then 
        do i_1=ncpl2(0,0)+1,icpl2
            do ia=1,ncpl3(i_1,0)
                iap=ia+1; if(iap.gt.ncpl3(i_1,0)) iap=1
                if(ncpl3(i_1,ia).eq.ncpl3(i_1,iap)) cycle  !始点と終点が同じものをのぞく 
                if(info(ncpl3(i_1,ia)).eq.2.and.info(ncpl3(i_1,iap)).eq.2) then 
                    ncmin=min(ncpl3(i_1,ia),ncpl3(i_1,iap))
                    ncmax=max(ncpl3(i_1,ia),ncpl3(i_1,iap))
                    icc=0
                    if(icc.eq.0) then
                        spath=spath+1
                        paths(spath,1)=ncmin
                        paths(spath,2)=ncmax
                        paths(spath,3)=1
                        paths(spath,-2:-1)=ncpl3(i_1,-2:-1) ! ポリゴンi_1の面種(-1)、面番号(-2)をコピー
#ifdef CK_SEGMENT
                        write(6,1229) spath,paths(spath,1:2)
1229                    format(1h ,'-------------水面パスNo=',i3,'水面パス始点終点=',i3,i3)
#endif
                    endif
                endif
            enddo
        enddo
endif

ncpl3(icpl2+1,-1)=2   !icpl2+1を水面として登録

    if(spath.eq.0) then !水面パス数spathが０のときの対処
        if(ipath.eq.1) then !水面を横切るパスが１の場合(おかしい)
             write(6,*) 'ipath_1=',ipath; goto 912
#ifdef SKIPNORMAL
             VV_1=-10000d0 !William 23 Aug, in order to continue calculation where possible
             return
#endif
        else if(ipath.ne.0) then!水面を横切るパスが１の場合(おかしい)
             write(6,*) 'ipath_2=',ipath; goto 912
#ifdef SKIPNORMAL
             VV_1=-10000d0 !William 23 Aug, in order to continue calculation where possible
             return
#endif
            tmpV=ZERO
            do ib=1,ipath
                tmpV=tmpV+pp(path(ib,1),:)
            enddo
            ncpl3(icpl2+1,0)=1
            ncpl3(icpl2+1,1)=ncpl3(0,1)+1
            pp(ncpl3(0,1)+1,:)=tmpV/dfloat(ipath)
            ncpl3(0,1)=ncpl3(0,1)+1       
            continue    
        else
            ncpl3(icpl2+1,0)=0
            continue
        endif
    else if(spath.eq.1) then !水面パス数spathが１のときの対処
        ncpl3(icpl2+1,0)=spath
        ncpl3(icpl2+1,1:2)=paths(1,1:2)
        continue
    else if(spath.eq.2) then !水面パス数spathが２のときの対処
     !   if(ipath.eq.3) then
            ncpl3(icpl2+1,0)=1
            ipp=ipp+1;ib1=0;pp(ipp,:)=ZERO
             do ib=1,ipp-1
                 if(vt_info(ib).eq.2) then
                 ib1=ib1+1
                 pp(ipp,:)=pp(ipp,:)+pp(ib,:)
                 endif
             enddo
             if(ib1.ne.0) then
                pp(ipp,:)=pp(ipp,:)/dfloat(ib1)
            ncpl3(icpl2+1,1)=ipp
            ncpl3(0,1)=ncpl3(0,1)+1
            continue
            else
            write(6,*) 'spath=',spath,icca,ic1a,surf,'nn=',nn,'ib1=',ib1
            call checkdataA 
            call checkdataB
            goto 912
           endif
    else
        ncpl3(icpl2+1,0)=spath
    endif

    if(spath.le.2) goto 345 ! 水面を構成できない場合、終了。
       
    paths(1:spath,3)=1 !　水面パスの登録状況を初期化。
356 continue
    icv(1)=paths(1,1) !　paths(1)の始点を水面ポリゴン並びの始点とする。
    icc1=0;icc3=0
359 do ia=1,ipp ! 各頂点が何本のパスに属するか調べる。
        icc=0
        do ib=1,spath
            if(paths(ib,3).eq.0) cycle        
            if(paths(ib,1).eq.ia) icc=icc+1
            if(paths(ib,2).eq.ia) icc=icc+1
        enddo
        if(icc.gt.2) then !三つ以上のパスに属する場合。頂点の順番をicv1に、何本に属するかをicv2に登録。該当する頂点数をicc1。
            icc1=icc1+1
            icv1(icc1)=ia 
            icv2(icc1)=icc
        else if(icc.eq.1) then !一つしか属さない場合、頂点の順番をicc3vに登録。該当する頂点の総数をicc3とする。
            icc3=icc3+1
            icc3v(icc3)=ia
        endif
    enddo

    if(icc3.gt.0) then  !一つしか属さない頂点がある場合、パスを削除する。
        do ib=1,spath
            if(paths(ib,3).eq.0)  cycle       
            if(paths(ib,1).eq.icc3v(1)) then ! 該当するパスを削除（利用済みとする）
                paths(ib,3)=0
            else if(paths(ib,2).eq.icc3v(1)) then
                paths(ib,3)=0
            endif
        enddo
        icc1=0;icc3=0
        goto 359 !もう一度、パスをチェック。
    endif

    ic4=10  
    allocate(route(1:ic4,0:spath+1),routem1(1:ic4,0:spath+1),routem2(1:ic4,0:spath+1))
    allocate(icv4(1:spath),icv5(1:spath))
    icv4(1:spath)=paths(1:spath,3)  !paths(:,3)の情報をicv4(:)にコピー
    ii_1=0
3245 continue
          
    if(sum(paths(1:spath,3)).lt.3) goto 6785 ! 水面パスが２以下の場合、6785へ
    
    paths(1:spath,4)=paths(1:spath,3)  !paths(:,3)の情報をpaths(:,4)にコピー

!水面パスrouteの始点を設定する。
    if(icc1.eq.0) then !　icc1が0のとき、paths(*,3)!=0のpathsを始点とする。
        do i_1=1,spath
            if(paths(i_1,3).ne.0) then
                ii_1=ii_1+1;
                routem1(ii_1,1)=paths(i_1,-1);routem2(ii_1,1)=paths(i_1,-2) ! パスが属する面の情報をroutem1,routem2にコピー。
                route(ii_1,1:2)=paths(i_1,1:2);paths(i_1,3:4)=0;goto 2344 ! パスの両端をrouteにコピー。paths(i_1)を使用済みにする。
            endif
        enddo
    else !　icc1が1以上のとき、一端がicv1(i_2)と一致するpaths(i_1)を始点とする。
        do i_1=1,spath
            i_3=0
            if(paths(i_1,3).eq.0) cycle        
            do i_2=1,icc1 ! 
                if(paths(i_1,1).eq.icv1(i_2)) then ! pathsの始点が一致する場合、そのままコピー。
                    ii_1=ii_1+1
                    routem1(ii_1,1)=paths(i_1,-1);routem2(ii_1,1)=paths(i_1,-2) 
                    route(ii_1,1:2)=paths(i_1,1:2);paths(i_1,3:4)=0;goto 2344
                else if(paths(i_1,2).eq.icv1(i_2)) then ! pathsの終点が一致する場合、始点終点を逆転してコピー。
                    ii_1=ii_1+1
                    routem1(ii_1,1)=paths(i_1,-1);routem2(ii_1,1)=paths(i_1,-2) 
                    route(ii_1,1)=paths(i_1,2);route(ii_1,2)=paths(i_1,1);paths(i_1,3:4)=0;goto 2344
                endif
            enddo         
        enddo
       ! 重複頂点がなくなったのでicc1=0
        icc1=0
        goto 3245
    endif        

    goto 6785

2344 i_6=2
     icv5(1:spath)=paths(1:spath,3) !paths(:,3)の情報をicv5(:)にコピー。
2341 i_4=0   

    if(route(ii_1,1).eq.route(ii_1,i_6)) then ! route(ii_1)の始点と終点が一致する場合（一周した場合）、route(ii_1)の頂点数を一減して、次のrouteを探す。
        route(ii_1,0)=i_6-1
        goto 2343
    endif
    do i_1=1,spath! route(ii_1)の終点が、path(i_1)の始or終点と一致する場合、追加パス番号(総数）をi_4、paths番号をicv7(i_4),終or始点の頂点番号をicv6(i_4)に設定。
        if(paths(i_1,4).eq.0) cycle
        if(route(ii_1,i_6).eq.paths(i_1,1)) then 
            i_4=i_4+1;icv7(i_4)=i_1;icv6(i_4)=paths(i_1,2)
        else if(route(ii_1,i_6).eq.paths(i_1,2)) then
            i_4=i_4+1;icv7(i_4)=i_1;icv6(i_4)=paths(i_1,1)
        endif
    enddo
    if(i_4.eq.1) then !追加パスが一つだけのとき、見つけた頂点をicvに加える。paths(icv7(i_4))は登録済みパスに。これに接続する次のpathsを探す。
        icv(i_6)=icv6(i_4);paths(icv7(i_4),3:4)=0  !;paths(icv7(i_4),3)=0
        routem1(ii_1,i_6)=paths(icv7(i_4),-1);routem2(ii_1,i_6)=paths(icv7(i_4),-2) !面情報をroute*にコピー。
        i_6=i_6+1; route(ii_1,i_6)=icv6(i_4) ! route(ii_!)に頂点を追加。
        goto 2341 !次の点orパスを探す。
    else if(i_4.eq.2) then !追加パスの候補が二つあるとき、
        tmpV=crossV(pp(route(ii_1,i_6-1),:)-pp(route(ii_1,i_6),:),cpl(routem2(ii_1,i_6-1),0:2) ) !route(ii_1,i_6)から一つ前の点route(ii_1,i_6-1)へのベクトルと!route(ii_1）を含む水面法線との外積
        i_5=0
        do i_3=1,i_4
            tmpV1=crossV(pp(icv6(i_3),:)-pp(route(ii_1,i_6),:),cpl(paths(icv7(i_3),-2),0:2) )!route(ii_1,i_6)から候補の点icv6(i_3)へのベクトルとicv6(i_3)を含む水面法線との外積
            if(sum(tmpV*tmpV1).lt.ZERO) then !外積の向きが逆になる場合。
                i_5=i_5+1
                icv8(i_5)=i_3
                continue
            endif
        enddo

        if(i_5.eq.1) then ! ひとつに絞れたとき。
            icv(i_6)=icv6(icv8(i_5));;paths(icv7(icv8(i_5)),4)=0 ;paths(icv7(icv8(i_5)),3)=0
            routem1(ii_1,i_6)=paths(icv7(icv8(i_5)),-1);routem2(ii_1,i_6)=paths(icv7(icv8(i_5)),-2) 
            i_6=i_6+1   
            route(ii_1,i_6)=icv6(icv8(i_5)) 
            goto 2341
        endif
! route(ii_1,i_6)を基点として
        tmpV=pp(route(ii_1,i_6-1),:)-pp(route(ii_1,i_6),:);tmpV=tmpV/sqrt(sum(tmpV**2)) ! 一つ前の点route(ii_1,i_6-1)への単位ベクトルをtmpV
        tmpV1=pp(icv6(icv8(1)),:)-pp(route(ii_1,i_6),:);tmpV1=tmpV1/sqrt(sum(tmpV1**2)) ! 候補1：icv6(icv8(1))への単位ベクトルをtmpV1
        tmpV2=pp(icv6(icv8(2)),:)-pp(route(ii_1,i_6),:);tmpV2=tmpV2/sqrt(sum(tmpV2**2)) ! 候補2：icv6(icv8(2))への単位ベクトルをtmpV2
        dd1=sum(tmpV*tmpV1)
        dd2=sum(crossV(tmpV,tmpV1)*surf(0:2))
        r1=cal_theta(dd1,dd2) ! tmpVとtmpV1のなす角をr1

        dd1=sum(tmpV*tmpV2)         
        dd2=sum(crossV(tmpV,tmpV2)*surf(0:2))
        r2=cal_theta(dd1,dd2)  ! tmpVとtmpV2のなす角をr2

!　なす角の小さいほうを次の頂点とする。
        dd1=r1;dd2=r2
        if(dd1.lt.dd2) then
            icv(i_6)=icv6(1);;paths(icv7(1),4)=0 ;paths(icv7(1),3)=0
            routem1(ii_1,i_6)=paths(icv7(1),-1);routem2(ii_1,i_6)=paths(icv7(1),-2) 
            i_6=i_6+1   
            route(ii_1,i_6)=icv6(1)            
            goto 2341
        else
            icv(i_6)=icv6(2);;paths(icv7(2),4)=0 ;paths(icv7(2),3)=0
            routem1(ii_1,i_6)=paths(icv7(2),-1);routem2(ii_1,i_6)=paths(icv7(2),-2) 
            i_6=i_6+1   
            route(ii_1,i_6)=icv6(2) 
            goto 2341
        endif
        continue
    else if(i_4.ge.3) then
        ttemp=sum(surf(0:2)*cpl(routem2(ii_1,i_6-1),0:2))
        cpltemp(0:2)=cpl(routem2(ii_1,i_6-1),0:2)-surf(0:2)*ttemp  ! cpltemp：route(ii_1,i_6-1)を含む面の法線ベクトルを水面に射影。
        cpltemp=cpltemp/sqrt(sum(cpltemp(0:2)**2)) !　cpltempの単位ベクトル化
        tmpV=crossV(pp(route(ii_1,i_6-1),:)-pp(route(ii_1,i_6),:),cpltemp(0:2)) ! 一つ前の点route(ii_1,i_6-1)へのベクトルとcpltempの外積
        !      dd1=dista(pp(icv6(1),:)-pp(route(ii_1,i_6-1),:))
        i_5=0
        do i_3=1,i_4
            ttemp=sum(surf(0:2)*cpl(paths(icv7(i_3),-2),0:2))
            cpltemp(0:2)=cpl(paths(icv7(i_3),-2),0:2)-surf(0:2)*ttemp ! cpltemp：paths(icv7(i_3))を含む面の法線ベクトルを水面に射影。
            cpltemp=cpltemp/sqrt(sum(cpltemp(0:2)**2))
            tmpV1=crossV(pp(icv6(i_3),:)-pp(route(ii_1,i_6),:),cpltemp(0:2) )  !候補i_3：icv6(icv6(i_3))へのベクトルとcpltempの外積
            if(sum(tmpV*tmpV1).lt.0) then
                i_5=i_5+1
                icv8(i_5)=i_3
                continue
            endif
        enddo

        if(i_5.eq.1) then ! tmpVと逆方向のベクトルが１つだけであれば、それを採用。
            icv(i_6)=icv6(icv8(1));;paths(icv7(icv8(1)),4)=0 ;paths(icv7(icv8(1)),3)=0
            routem1(ii_1,i_6)=paths(icv7(icv8(1)),-1);routem2(ii_1,i_6)=paths(icv7(icv8(1)),-2) 
            i_6=i_6+1;route(ii_1,i_6)=icv6(icv8(1)) 
            goto 2341                
        else if(i_5.gt.10) then !　10以上あるとき　912へ（エラー）
            write(6,*) 'i_5=',i_5
            goto 912
        else if(i_5.eq.0) then  !   0のとき、icv8を初期化して、次に。
            i_5=i_4
            icv8(1:10)=(/1,2,3,4,5,6,7,8,9,10/)
        else

        
        ! tmpVと逆方向のベクトルが2-9この場合、route(ii_1,i_6)から頂点候補へのベクトルとroute(ii_1,i_6)から一つ前の点route(ii_1,i_6-1)へのベクトルのなす角を計算
        rr=100000.0d0
        tmpV=pp(route(ii_1,i_6-1),:)-pp(route(ii_1,i_6),:);tmpV=tmpV/sqrt(sum(tmpV**2))
        do ii_2=1,i_5
        tmpVV(ii_2,:)=pp(icv6(icv8(ii_2)),:)-pp(route(ii_1,i_6),:);tmpVV(ii_2,:)=tmpVV(ii_2,:)/sqrt(sum(tmpVV(ii_2,:)**2))
        dd1=sum(tmpV*tmpVV(ii_2,:))
        dd2=sum(crossV(tmpV,tmpVV(ii_2,:))*surf(0:2))
        rr(ii_2)=cal_theta(dd1,dd2)
        enddo

        ! 最小を与える頂点を採用する。
        allocate(imin(1:i_5))
        imin(1:i_5)=MINLOC(rr(1:i_5))
        irmin=icv8(imin(1));deallocate(imin)
  !      if(dd1.le.dd2.and.dd1.le.dd3) then
            icv(i_6)=icv6(irmin);;paths(icv7(irmin),4)=0 ;paths(icv7(irmin),3)=0
            routem1(ii_1,i_6)=paths(icv7(irmin),-1);routem2(ii_1,i_6)=paths(icv7(irmin),-2) 
            i_6=i_6+1   
            route(ii_1,i_6)=icv6(irmin) 
            goto 2341  
        endif

    else 
        write(6,*) 'i_4=',i_4,'i_55=',i_5
        goto 912
        ii_1=ii_1-1
        paths(1:spath,3)=icv5(1:spath)

    endif
2343 continue    
    goto 3245
6785 continue
     
    ib4=ii_1;
927 do i_1=1,ii_1 !ひとつのrouteに同じ頂点が二回出現する場合、その間にある点から新しいrouteを作成するとともに、抜けた分はつなぎ合わせる。
        nnp=route(i_1,0) 
        do iie=1,nnp
        do iis=iie-1,1,-1
            if(route(i_1,iis).eq.route(i_1,iie)) then
                 ib4=ib4+1
                 route(ib4,1:iie-iis)=route(i_1,iis+1:iie);route(ib4,0)=iie-iis !　新しいrouteを作成する
                 routem1(ib4,1:iie-iis)=routem1(i_1,iis+1:iie)  
                 routem2(ib4,1:iie-iis)=routem2(i_1,iis+1:iie) 
                 
                 route(i_1,iis+1:nnp-iie+iis)=route(i_1,iie+1:nnp);route(i_1,0)=nnp-iie+iis !　抜けた分はつなぎ合わせる
                 routem1(i_1,iis+1:nnp-iie+iis)=routem1(i_1,iie+1:nnp)
                 routem2(i_1,iis+1:nnp-iie+iis)=routem2(i_1,iie+1:nnp)
                 goto 936
            endif
        enddo
        enddo
936 continue        
    enddo
    if(ii_1.ne.ib4) then !　該当がなくなるまで繰り返す。
    ii_1=ib4
    goto 927
    endif
     
! 見つけたrouteについて、法線の向きを確定して、ncpl3に登録する。     
    ii=icpl2
    do i_1=1,ii_1 
        nnp=route(i_1,0)
        if(nnp.lt.3) cycle
        ii=ii+1 ! ii_1に応じて面を増やす。
        call surface_dir3(nnp,route(i_1,1:nnp),surf(:),ncpl3(ii,1:nnp)) !法線ベクトルが水面と同じ向きになるように並べ替える。ncpl3(ii)を作成
        ib=1;
591     ibp=ib+1;if(ib.eq.nnp) ibp=1
      !         do i_4=1,nnp! route(i_1,ib:ib3)がncpl3(ii,i_4:i_2)と一致する場合（反転も可）、それをncpl3(ii,ib1:ib2)に。
      !              i_4p=i_4+1; if(i_4.eq.nnp) i_4p=1
      !              if(ncpl3(ii,i_4).eq.route(i_1,ib).and.ncpl3(ii,i_4p).eq.route(i_1,ib3)) then 
      !                  ib1=i_4;ib2=i_4p
      !          else if(ncpl3(ii,i_4p).eq.route(i_1,ib).and.ncpl3(ii,i_4).eq.route(i_1,ib3)) then !
      !                  ib1=i_4;ib2=i_4p
      !              endif
      !          enddo
      !  ib1=ib;ib2=ib3
                icc=0
                do i_3=1,icpl2 !ncpl3(ii,ib1:ib2)が ncpl3(i_3,i_4:i_2)と一致する場合（反転も可、それをncpl3(ib3,ib4:ib4+1)とする。
                    do i_4=1,ncpl3(i_3,0)
                        i_4p=i_4+1; if(i_4.eq.ncpl3(i_3,0)) i_4p=1
                        if(ncpl3(i_3,i_4).eq.ncpl3(ii,ib).and.ncpl3(i_3,i_4p).eq.ncpl3(ii,ibp)) then ! (ib3,ib4:i_2)に。
                            icc=icc+1;ib5=-1    !ib3=i_3;ib4=i_4;
                        else if(ncpl3(i_3,i_4p).eq.ncpl3(ii,ib).and.ncpl3(i_3,i_4).eq.ncpl3(ii,ibp)) then
                            icc=icc+1;ib5=1      !;ib3=i_3;ib4=i_4
                        endif
                    enddo
                enddo
                if(icc.ne.1) then !　icc=1となるまで繰り返す。
                    if(ib.eq.nnp) then
                        write(6,*) 'decision of surface dimection'
                        goto 912
                    endif
                    ib=ib+1
                    goto 591
                endif
    !           i_3=ib3;i_4=ib4 ! ncpl3(ib3,ib4:ib4+1)とncpl3(ii,ib1:ib2)が同じ向きの場合、水面ではないので、法線ベクトルを反転する。
     !                   i_2=i_4+1; if(i_4.eq.ncpl3(i_3,0)) i_2=1
         !               if(ncpl3(i_3,i_4).eq.ncpl3(ii,ib1).and.ncpl3(i_3,i_2).eq.ncpl3(ii,ib2)) then
                if(ib5.eq.-1) then
                            surf1(:)=surf(:)*(-1.0d0)
                            call surface_dir3(nnp,route(i_1,1:nnp), surf1(:),ncpl3(ii,1:nnp));ncpl3(ii,-2)=-2
                endif            
          !              endif
        
 
437 ncpl3(ii,0)=nnp ;ncpl3(ii,-1)=2 ! 面の頂点数、面種を設定。
    enddo
    ncpl3(0,0)=ii; !セルの面数を確定。
  
! ポリゴン面が重複している（完全に一致した面が二つある）場合削除する。
345 do i_1=1,icpl2
        do ii=icpl2+1,ncpl3(0,0)
            if(ncpl3(ii,-1).ne.2) cycle
            if(ncpl3(i_1,0).gt.ncpl3(ii,0)) cycle
            icc=0
            ia1=1
            do ib1=1,ncpl3(i_1,0)
                if(ncpl3(i_1,ib1).eq.ncpl3(ii,ia1)) goto 433
            enddo
            cycle
433         icc=1
            do while(icc < ncpl3(ii,0)) 
                ib1=ib1+1; if(ib1.gt.ncpl3(i_1,0)) ib1=1
                ia1=ia1+1; if(ia1.gt.ncpl3(ii,0)) ia1=1
                if(ncpl3(i_1,ib1).ne.ncpl3(ii,ia1)) exit
                icc=icc+1
            end do
            if(icc.eq.ncpl3(i_1,0)) then
                !ncpl3(0,0)=icpl2
                ncpl3(i_1,0)=0
                !write(6,*) 'suimen chofuku'
            endif
        enddo
    enddo

!　水面ポリゴン以外で、水面頂点のみで構成される面があれば削除。
    do i_1=1,ncpl3(0,0)
        if(ncpl3(i_1,-1).eq.2) cycle
        if(ncpl3(i_1,0).eq.0) cycle
        do ia=1,ncpl3(i_1,0)
            if(vt_info(ncpl3(i_1,ia)).ne.2) goto 342
        enddo
           ncpl3(i_1,0)=0
342     continue
    enddo

! 形状Bの体積を計算する。
346 VV(myhr)=ZERO
#ifdef BEDMOVE
    if(ground.eq.1) then
        axmm(nn)=ZERO;aymm(nn)=ZERO;azmm(nn)=ZERO
        axpp(nn)=ZERO;aypp(nn)=ZERO;azpp(nn)=ZERO
    endif
#elif HENDO
    axyz_sd=ZERO
#endif
    do i_1=1,ncpl3(0,0) ! 体積を計算
        nnp=ncpl3(i_1,0)
        if(nnp.eq.0) cycle
        if(nnp.ne.1.and.ncpl3(i_1,1).eq.ncpl3(i_1,nnp)) then
            nnp=nnp-1;ncpl3(i_1,0)=nnp
            continue    
        endif
        tmpV=ZERO
        do ib=2,nnp-1
            ib1=ncpl3(i_1,ib)
            ib2=ncpl3(i_1,ib+1)
            tmpV=crossV(pp(ib1,0:2)-pp(ncpl3(i_1,1),0:2),pp(ib2,0:2)-pp(ncpl3(i_1,1),0:2))+tmpV
        enddo

#ifdef BEDMOVE
        if(ground.eq.1) then
            if(ncpl3(i_1,-1).eq.0) then
                if(ncpl3(i_1,-2).eq.1) then
                    azmm(nn)=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(1,jn(nn))
                else if(ncpl3(i_1,-2).eq.2) then
                    azpp(nn)=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(1,jn(nn))
                else if(ncpl3(i_1,-2).eq.3) then
                    aymm(nn)=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(2,kn(nn))
                else if(ncpl3(i_1,-2).eq.4) then
                    axpp(nn)=sqrt(sum(tmpv(0:2)**2))*HALF/dx(1,jn(nn))/dx(2,kn(nn))
                else if(ncpl3(i_1,-2).eq.5) then
                    aypp(nn)=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(2,kn(nn))
                else if(ncpl3(i_1,-2).eq.6) then
                    axmm(nn)=sqrt(sum(tmpv(0:2)**2))*HALF/dx(1,jn(nn))/dx(2,kn(nn))
                endif
                continue
            endif
        endif
#elif HENDO
        if(ncpl3(i_1,-1).eq.0) then
            !axyz_sd(ncpl3(i_1,-2))=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(1,jn(nn))
            if(ncpl3(i_1,-2).eq.1.or.ncpl3(i_1,-2).eq.2) then
                axyz_sd(ncpl3(i_1,-2))=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(1,jn(nn))
                !else if(ncpl3(i_1,-2).eq.2) then
                !axyz_sd(ncpl3(i_1,-2))=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(1,jn(nn))
            else if(ncpl3(i_1,-2).eq.3.or.ncpl3(i_1,-2).eq.5) then
                axyz_sd(ncpl3(i_1,-2))=sqrt(sum(tmpv(0:2)**2))*HALF/dx(0,in(nn))/dx(2,kn(nn))
            else
                axyz_sd(ncpl3(i_1,-2))=sqrt(sum(tmpv(0:2)**2))*HALF/dx(1,jn(nn))/dx(2,kn(nn))
            endif
        endif
#endif
        VV(myhr)=VV(myhr)+sum(tmpV*pp(ncpl3(i_1,1),0:2))*SIXTH
    enddo
    ic7=0
    do i_1=1,ncpl3(0,0)!すべての辺についてある面の辺が、他の面の辺と逆向きに接するかを確認。
        if(ncpl3(i_1,-1).eq.2.and.ncpl3(i_1,0).le.2) cycle
        if(ncpl3(i_1,0).eq.0) cycle
        ic8=0
        do ia=1,ncpl3(i_1,0)
            ib=ia+1;if(ib.gt.ncpl3(i_1,0)) ib=1

            do ii=1,ncpl3(0,0)
                if(i_1.eq.ii) cycle
                do ic=1,ncpl3(ii,0)
                    ic1=ic+1;if(ic1.gt.ncpl3(ii,0)) ic1=1
                    if(ncpl3(i_1,ib).eq.ncpl3(ii,ic).and.ncpl3(i_1,ia).eq.ncpl3(ii,ic1)) then !面i_1の辺が、他の面の辺と逆向きに接するかを確認。
                        goto 324
                    endif
                enddo                    
            enddo    
            ic7=1;ic8=ic8+1
            !       write(6,*) 'mondaiari(DriSeg3 nn,n=',nn,ia,ib
#ifndef CK_SEGMENT 
!            write(6,*) 'ijk=',in(nn),jn(nn),kn(nn)
#endif
#ifndef SKIPNORMAL
            write(6,*) 'i_1,ncpl3(i_1,ia),ncpl3(i_1,ia)=',i_1,ncpl3(i_1,ia),ncpl3(i_1,ib)
#endif
!    write(6,*) ncpl3(1:ncpl3(0,0),0)

324         continue        
        enddo
#ifdef BUGWIN
        if(ic8.eq.ncpl3(i_1,0)) then
        continue
        endif
#endif
    enddo
    if(ic7.eq.1) then ! ある面の辺が他の面の辺と逆向きに接しないとき、エラー(912へ）
#ifdef SKIPNORMAL
        VV_1=-10000d0 !William 23 Aug, in order to continue calculation where possible
        return
#endif
        write(6,*) 'ic7=',ic7,'nn=',nn
        goto 912
    endif
3241 continue
    !#endif
    VV_1=VV(myhr)

#ifdef CK_SEGMENT
    call checkdataA
    call checkdataB
    write(6,*) 'mondai nashi'
#endif
    return

912 continue
    call checkdataA 
    call checkdataB
#ifdef OBJECT
#ifndef CK_SEGMENT  
    call output_bindata
#endif    
#endif
    stop
    contains
#ifdef OBJECT    
#ifndef CK_SEGMENT  
 subroutine output_bindata
    use arrays, only: idri_kdc, idri_fc_n, upl, upl0, f, fb, in, jn, kn
    use variables, only: n, npln, t, ncpl00max, ncpl01max, ncpl02max
    implicit none
!    integer::i1,i2
    open(29,file='mondaiari.bin',form='unformatted')
    write(29) ncpl00max,ncpl01max,ncpl02max
    write(29) nn,vv_1,surf,pp,ncpl2
    write(29) npln
    write(29) cpl
    i_1=ubound(idri_kdc,1);i_2=ubound(idri_fc_n,1)
    write(29) i_1,i_2
    write(29) idri_kdc(1:i_1),idri_fc_n(1:i_2)
    i_1=ubound(upl0,1)
    write(29) i_1
 !   write(29) upl0(1:i_1,1:3)

    i_1=ubound(upl,1);i_2=ubound(upl,2)
    write(29) i_1,i_2
   ! if(i_2.ne.0) then        
   !     write(29) upl(1:i_1,1:i_2,1:3)
   ! endif
    write(6,*) 'stop in PlainAdjust',t,n,spath
    write(6,*) ncpl3(1:ncpl3(0,0),0),VV(myhr)
 !   write(6,*) 'icc0=',icc0,'n=',n
    close(29)
 endsubroutine
#endif
#endif
    doubleprecision function cal_theta(d1,d2)
    use variables, only: HALF
    implicit none
    doubleprecision::d1,d2

    if(abs(d1).lt.1.0d-12) then
        if(d2.gt.ZERO) then
            cal_theta=pi*HALF
        else if(d2.lt.ZERO) then
            cal_theta=pi*1.5d0
        else
            continue
        endif
    else if(d1.gt.ZERO) then
        if(d2.gt.ZERO) then
            cal_theta=atan(d2/d1)
        else
            cal_theta=2.0*pi-atan(d2/d1)
        endif
    else
        if(d2.gt.ZERO) then
            cal_theta=pi+atan(d2/d1)
        else
            cal_theta=pi-atan(d2/d1)
        endif
    endif

    end function

    subroutine checkdataA
    write(611,*) nn,ipp  !,ncpl2(0,1),objno(nn)
    do ib=1,ipp  ! ncpl2(0,1)
        write(611,11) pp(ib,:)
11      format(1h ,e12.5,10(' ',e12.5))
        continue
    enddo
    write(611,*) ncpl2(0,0)
    do ib=1, ncpl2(0,0)
        write(611,12) ncpl2(ib,-2:ncpl2(ib,0))
    enddo
12  format(1h ,i7,12(' ',i2))
    !write(611,*) 'nn=',nn,'i,k,j=',in(nn),kn(nn),jn(nn)
    end subroutine
    subroutine checkdataB
    write(613,*) nn,ncpl3(0,1)
    do ib=1, ncpl3(0,1)
        write(613,11) pp(ib,:)
11      format(1h ,e12.5,10(' ',e12.5))
        continue
    enddo
    write(613,*) ncpl3(0,0)
    do ib=1, ncpl3(0,0)
        write(613,12) ncpl3(ib,-2:ncpl3(ib,0))
12      format(1h ,i7,12(' ',i2))
    enddo
    end subroutine
    !#endif
    doubleprecision function surface_dir23(nnp,icv,cpls)
    implicit none
    doubleprecision,dimension(0:2)::tmpV
    doubleprecision,dimension(0:2),intent(in)::cpls
    !integer,dimension(1:nnp)::ncp
    integer,dimension(1:nnp)::icv
!    doubleprecision::dd2,dis1,dis2,dista
    integer,intent(in)::nnp
    integer::ib !,ic
    tmpV=ZERO
    do ib=2,nnp-1
        tmpV=crossV(pp(icv(ib),0:2)-pp(icv(1),0:2),pp(icv(ib+1),0:2)-pp(icv(1),0:2))+tmpV
    enddo
    surface_dir23=dot_product(tmpV,cpls(0:2))
    end function

    subroutine surface_dir3(nnp,icv,cpls,ncp)
    implicit none
    doubleprecision,dimension(0:2)::tmpV
    doubleprecision,dimension(0:2),intent(in)::cpls
    integer,dimension(1:nnp)::ncp
    integer,dimension(1:nnp)::icv
    doubleprecision::dd2 !,dis1,dis2,dista
    integer,intent(in)::nnp
    integer::ib,ic

    tmpV=0.0

    tmpV=ZERO
    do ib=2,nnp-1
        tmpV=crossV(pp(icv(ib),0:2)-pp(icv(1),0:2),pp(icv(ib+1),0:2)-pp(icv(1),0:2))+tmpV
    enddo


    dd2=dot_product(tmpV,cpls(0:2))
    if(dd2.gt.ZERO) then
        ncp(1:nnp)=icv(1:nnp)
    else
        do ic=1,nnp
            ncp(ic)=icv(nnp-ic+1)
        enddo
    endif   
    end subroutine
    
    
end subroutine PlainAdjust
