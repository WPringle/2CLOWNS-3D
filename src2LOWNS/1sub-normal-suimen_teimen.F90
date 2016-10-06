subroutine suimen_teimen(nn,x0,dx,nrm,spp,fnn,iposi)
    use variables, only: VSMALL, LEPS, SIXTH, THIRD, FTW7, TWFTH, HALF, ZERO
    implicit none
    Integer,intent(in) :: iposi,nn
    real*8,intent(in) :: fnn
    real*8,intent(in),dimension(0:2) :: x0,dx
    real*8,dimension(0:2),intent(out) :: spp
    Integer,dimension(0:8) :: il
    real*8,dimension(0:8,0:2) :: XL
    real*8,dimension(0:8) :: L
    real*8,dimension(0:2) :: nrm,oa,ob,oc,p1,p2,p3,p4,p5,p6 !,od
    real*8 :: sc,tc,uc,ss,uci,uc1,uc2,ucc,uccp,vvd !,sci,tci
    real*8 :: dm1,dm2,dmc,dmcp,axn,bxn,cxn,alpha,beta
    real*8 :: a2,b2,d2,c1,vv,x1,x2,a01,a02,a03 
    real*8 :: Ln,S1,S2,S3,S4,S5,S6,pp1,qq,r3a,r3b
    integer :: iflag,iconk,icc,Nosg,ic,i1,i2

dmcp = ZERO
sc = ZERO; tc = ZERO ; uc = ZERO; Nosg = 0
!Will Sep 2 -> Keep x0 definition and change xp to dx
if (abs(nrm(0)).lt.VSMALL.and.abs(nrm(1)).lt.VSMALL) then
    Nosg = 1; ic = 2; i1 = 0; i2 = 1
elseif (abs(nrm(2)).lt.VSMALL.and.abs(nrm(1)).lt.VSMALL) then
    Nosg = 1; ic = 0; i1 = 1; i2 = 2
elseif (abs(nrm(0)).lt.VSMALL.and.abs(nrm(2)).lt.VSMALL) then
    Nosg = 1; ic = 1; i1 = 0; i2 = 2
endif
!<<<<<<<<<<<<<<<<<<<<<<<1D Calculation<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (Nosg.eq.1) then
    spp(i1) = x0(i1) + dx(i1) * HALF ; spp(i2) = x0(i2) + dx(i2) * HALF
    if (nrm(ic).gt.ZERO) then
        spp(ic) = x0(ic) + dx(ic) * fnn  
    else
        spp(ic) = x0(ic) + dx(ic) * (1.0d0-fnn)
    endif
    return
endif
!<<<<<<<<<<<<<<<<<<<<<<<2D Calculation<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (abs(nrm(0)).lt.VSMALL.or.abs(nrm(1)).lt.VSMALL.or.abs(nrm(2)).lt.VSMALL) then 
    if (abs(nrm(0)).lt.VSMALL) ic = 0 !No surface gradient in x direction
    if (abs(nrm(1)).lt.VSMALL) ic = 1 !No surface gradient in y direction
    if (abs(nrm(2)).lt.VSMALL) ic = 2 !No surface gradient in z direction
    nrm = nrm_nrm(nrm,VSMALL)
    call make_XL_ob(ic)
    Ln = sqrt(sum(nrm(:)**2))
    L(1:8) = (nrm(2)*(XL(1:8,2)-XL(0,2))+nrm(0)*(XL(1:8,0)-XL(0,0))            &
           +  nrm(1)*(XL(1:8,1)-XL(0,1)))/Ln;    
    call narabe4(L,il)
    oa(:)=XL(il(2),:)-XL(0,:)
    oc(:)=XL(il(3),:)-XL(0,:) 
    axn = sum(oa(:)*nrm(:)) ; bxn = sum(ob(:)*nrm(:));cxn=sum(oc(:)*nrm(:))
    S1=HALF*axn/cxn
    S2=1.0d0-S1
    if (fnn.lt.S1) then
        sc=sqrt(2.0d0*fnn*cxn/axn) 
        uc=sc*axn/cxn
        vv=taiseki2(sc,uc)
        if(iposi.eq.1) then
            p1(:)=uc*oc(:)
            p2(:)=sc*oa(:)
            p3(:)=uc*oc(:)+ob(:)
            p4(:)=sc*oa(:)+ob(:)
            spp(:)=0.25d0*(p1(:)+p2(:)+p3(:)+p4(:))+XL(0,:)
        endif        
    elseif (fnn.lt.S2) then
        sc=fnn*(cxn/axn)+HALF 
        uc=sc*axn/cxn
#ifdef BUGWIN        
        vv=taiseki2(sc,uc)
#endif 
        if (iposi.eq.1) then
            p1(:)=uc*oc(:)
            p2(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)
            p3(:)=uc*oc(:)+ob(:)
            p4(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)+ob(:)
            spp(:)=0.25d0*(p1(:)+p2(:)+p3(:)+p4(:))+XL(0,:)
        endif        
    elseif (fnn.lt.Leps) then
        a01=1.0d0;a02=-2.0d0*(1.d0+cxn/axn)
        a03=2.0d0*fnn*cxn/axn+(cxn/axn)**2+1.0d0
        call kai2(a01,a02,a03,x1,x2,iconk)
        if(x1.lt.1.0d0+cxn/axn) then
            sc=x1
        elseif (x2.lt.1.0d0+cxn/axn) then
            sc=x2
        elseif (x1.eq.x2) then
            sc=x1
        endif
        uc=sc*axn/cxn
#ifdef BUGWIN        
        vv=taiseki2(sc,uc)
#endif        
        if(iposi.eq.1) then
            p1(:)=oc(:)+sc*(1.0d0-1.0d0/uc)*oa(:)
            p2(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)
            p3(:)=p1(:)+ob(:)
            p4(:)=p2(:)+ob(:)
            spp(:)=0.25d0*(p1(:)+p2(:)+p3(:)+p4(:))+XL(0,:)
        endif        
    else
        !write(6,*) 'Bad IF statement: in suimen_teimen line 107'   
        !write(6,*) fnn,nrm  
        return 
    endif   
!<<<<<<<<<<<<<<<<<<<<<<<Full 3D Calculation<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
else 
    call make_XL_ob(-1)
    Ln=sqrt(sum(nrm(:)**2))
    L(1:8)=(nrm(2)*(XL(1:8,2)-xl(0,2))+nrm(0)*(XL(1:8,0)-xl(0,0))              &
          +nrm(1)*(XL(1:8,1)-xl(0,1)))/Ln; 
    call narabe(L,il)
    oa(:)=xl(il(2),:)-xl(0,:)
    ob(:)=xl(il(3),:)-xl(0,:)
    oc(:)=xl(il(4),:)-xl(0,:)
    c1=abs(sum(oa(:)*oc(:)))+abs(sum(ob(:)*oc(:)))
    iflag=0
        if (c1.ne.ZERO) then
            oc(:)=xl(il(5),:)-xl(0,:)
            iflag=1
        endif
    axn=sum(oa(:)*nrm(:));bxn=sum(ob(:)*nrm(:));cxn=sum(oc(:)*nrm(:))
    alpha=cxn/axn;beta=cxn/bxn
    s1=axn**2/(bxn*cxn*6.0d0)
    a2=3.d0*bxn**2;b2=-3.0d0*bxn*axn;d2=axn**2
    S2=(axn**2+3.0d0*bxn**2-3.0d0*bxn*axn)/(bxn*cxn*6.0d0)! ’¸“_B‚ð‚Æ”„‚é‚Æ‚« t=1
    S5=1.0d0-S2
    S6=1.0d0-S1
    if(iflag.eq.1) then
        S3=HALF*(bxn+axn)/cxn
        S4=1.0d0-S3
    else
        S3=1.0d0-HALF*(1.0d0/alpha+1.0d0/beta)                                &
          -SIXTH*alpha*beta*(1.0d0-1.0d0/alpha-1.0d0/beta)**3
        uc=1.0d0
        sc=uc*alpha
        tc=uc*beta
        S4=1.0d0-S3
        uc=1.0d0/alpha+1.0d0/beta
        sc=uc*alpha
        tc=uc*beta
    endif
    if(fnn.le.S1) then
        uc=(6.0d0*fnn/alpha/beta)**(THIRD)
        sc=uc*cxn/axn
        tc=uc*cxn/bxn
        if(iposi.eq.1) then
            spp(:)=THIRD*(sc*oa(:)+tc*ob(:)+uc*oc(:))+XL(0,:)
        endif
    else if(fnn.le.S2) then
        uc=HALF/alpha+sqrt(2.0d0*fnn/beta-1.0d0/12.0d0/alpha**2)
        sc=uc*cxn/axn
        tc=uc*cxn/bxn
        if(iposi.eq.1) then
            p1(:)=uc*oc(:)
            p2(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)
            p3(:)=tc*ob(:)
            p4(:)=oa(:)+tc*(1.0d0-1.0d0/sc)*ob(:)
            spp(:)=0.25d0*(p1(:)+p2(:)+p3(:)+p4(:))+XL(0,:)
        endif
    else if(fnn.le.S3) then
        pp1=-6.0d0/alpha/beta
        qq=3.0d0*(2.0d0*fnn-1.0d0/alpha-1.0d0/beta)/alpha/beta
        ss=qq**2+FTW7*pp1*pp1*pp1
        if(ss.ge.ZERO) then
            r3a=HALF*(qq+sqrt(ss))
            r3b=HALF*(qq-sqrt(ss))
            uc=-r3a**(THIRD)-sign(1.0d0,r3b)*abs(r3b)**(THIRD)+1.0d0/alpha+1.0d0/beta
        else
            ucc=bxn/cxn;dmc=fnn-S2
            vvd=-HALF*alpha*beta*( ucc**2-(ucc-1.0d0/alpha)**2-(ucc-1.0d0/beta)**2 )
            a02=-2.0d0*(1.0d0/alpha+1.0d0/beta);a03=1.0d0/alpha**2+1.0d0/beta**2
            ucc=ucc-dmc/vvd
            dmc=fnn-SIXTH*alpha*beta*( ucc**3-(ucc-1.0d0/alpha)**3-(ucc-1.0d0/beta)**3 )
            dmcp=ZERO;icc=0
            231 vvd=-HALF*alpha*beta*( ucc**2-(ucc-1.0d0/alpha)**2-(ucc-1.0d0/beta)**2 )
            ucc=ucc-dmc/vvd
            dmc=fnn-SIXTH*alpha*beta*( ucc**3-(ucc-1.0d0/alpha)**3-(ucc-1.0d0/beta)**3 )
            if(abs(dmc).lt.1.0d-15) then
                uc=ucc
            else if(icc.gt.10) then
                if(abs(dmc).lt.1.0d-10) then
                    uc=uccp
                else
                    write(61,*) 'nashi',dmc,dmcp
                    write(61,*) x0,dx,nrm,fnn,iposi                    
                    write(62) x0,dx,nrm,fnn,iposi
                    stop
                endif
            else
                icc=icc+1
                dmcp=dmc
                uccp=ucc
                goto 231
            endif
        endif
        sc=uc*alpha
        tc=uc*beta
        vv=SIXTH*alpha*beta*uc**3*(1.0d0-(1.0d0-1.0d0/tc)**3-(1.0d0-1.0d0/sc)**3)
        if(iposi.eq.1) then
            if(uc.ne.ZERO) then
                p1(:)=uc*oc(:)
                p2(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)
                p3(:)=ob(:)+uc*(1.0d0-1.0d0/tc)*oc(:)
                p4(:)=oa(:)+tc*(1.0d0-1.0d0/sc)*ob(:)
                p5(:)=ob(:)+sc*(1.0d0-1.0d0/tc)*oa(:)
            else
                p1(:)=uc*oc(:)
                p2(:)=oa(:)+(uc-1.0d0/alpha)*oc(:)
                p3(:)=ob(:)+(uc-1.0d0/beta)*oc(:)
                p4(:)=oa(:)+(tc-beta/alpha)*ob(:)
                p5(:)=ob(:)+(sc-alpha/beta)*oa(:)
            endif
            spp(:)=0.2d0*(p1(:)+p2(:)+p3(:)+p4(:)+p5(:))+XL(0,:)
        endif
    elseif(fnn.le.S4) then
        if(iflag.eq.1) then
            uc=fnn+HALF*(bxn+axn)/cxn
            sc=uc*alpha
            tc=uc*beta
            if(iposi.eq.1) then
                p1(:)=uc*oc(:)
                p2(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)
                p3(:)=ob(:)+uc*(1.0d0-1.0d0/tc)*oc(:)
                p4(:)=oa(:)+ob(:)+uc*(1.0d0-1.0d0/tc-1.0d0/sc)*oc(:)
                spp(:)=0.25d0*(p1(:)+p2(:)+p3(:)+p4(:))+XL(0,:)
            endif
        else
            uc1=1.0d0
            dm1=fnn-SIXTH*alpha*beta*( 3.0d0*uc1**2-3.0d0*uc1+1.0d0-(uc1-1.0d0/alpha)**3-(uc1-1.0d0/beta)**3 )
            uc2=1.0d0/alpha+1.0d0/beta
            dm2=fnn-SIXTH*alpha*beta*( 3.0d0*uc2**2-3.0d0*uc2+1.0d0-(uc2-1.0d0/alpha)**3-(uc2-1.0d0/beta)**3 )
            if (dm1*dm2.gt.ZERO) then
                continue
            else
                232 ucc=(uc1*dm2-uc2*dm1)/(dm2-dm1)
                dmc=fnn-( alpha*beta-2.0d0*alpha*beta*ucc**3+3.0d0*(alpha*beta+alpha+beta)*ucc**2&
                                    -3.0d0*(alpha*beta+beta/alpha+alpha/beta)*ucc+beta/alpha**2+alpha/beta**2 )/6.0d0
                if (abs(dmc).lt.1.0d-15) then
                    uc=ucc
                    sc=uc*cxn/axn
                    tc=uc*cxn/bxn
                    if(iposi.eq.1) then
                        p1(:)=oc(:)+sc*(1.0d0-1.0d0/uc)*oa(:)
                        p2(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)
                        p3(:)=ob(:)+uc*(1.0d0-1.0d0/tc)*oc(:)
                        p4(:)=oa(:)+tc*(1.0d0-1.0d0/sc)*ob(:)
                        p5(:)=ob(:)+sc*(1.0d0-1.0d0/tc)*oa(:)
                        p6(:)=oc(:)+tc*(1.0d0-1.0d0/uc)*ob(:)
                        spp(:)=SIXTH*(p1(:)+p2(:)+p3(:)+p4(:)+p5(:)+p6(:))+XL(0,:)
                    endif
                elseif(dmc.eq.dmcp) then
                    if(abs(dmc).lt.1.0d-13) then
                        uc=ucc
                        sc=uc*cxn/axn
                        tc=uc*cxn/bxn
                        if(iposi.eq.1) then
                            p1(:)=oc(:)+sc*(1.0d0-1.0d0/uc)*oa(:)
                            p2(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:)
                            p3(:)=ob(:)+uc*(1.0d0-1.0d0/tc)*oc(:)
                            p4(:)=oa(:)+tc*(1.0d0-1.0d0/sc)*ob(:)
                            p5(:)=ob(:)+sc*(1.0d0-1.0d0/tc)*oa(:)
                            p6(:)=oc(:)+tc*(1.0d0-1.0d0/uc)*ob(:)
                            spp(:)=SIXTH*(p1(:)+p2(:)+p3(:)+p4(:)+p5(:)+p6(:))+XL(0,:)
                        endif
                    endif
                else
                    if(dmc*dm1.lt.ZERO) then
                        dm2=dmc
                        uc2=ucc
                    else if(dmc*dm2.lt.ZERO) then
                        uc1=ucc
                        dm1=dmc
                    endif
                    dmcp=dmc
                    goto 232
                endif
            endif
        endif
    elseif (fnn.le.S5) then
        if(iflag.eq.0) then
            ucc=1.0d0/alpha+1.0d0/beta
        else
            ucc=1.0d0
        endif
        dmc=fnn-S4
        dmcp=ZERO
        233 vvd=alpha*beta/2.0d0*(ucc-1.0d0)**2-1.0d0
        ucc=ucc-dmc/vvd
        dmc=fnn-HALF*(2.0d0*ucc-1.0d0/alpha-1.0d0/beta)+SIXTH*alpha*beta*(ucc-1.0d0)**3
        if(abs(dmc).lt.1.0d-15) then
            uc=ucc
        else if(dmcp.ne.ZERO.and.abs(dmc).ge.abs(dmcp)) then
            continue
        else
            dmcp=dmc
            goto 233
        endif
        sc=uc*cxn/axn
        tc=uc*cxn/bxn
        if(iposi.eq.1) then
            p1(:)=oc(:)+sc*(1.0d0-1.0d0/uc)*oa(:) !CEŠÔ
            p2(:)=oc(:)+tc*(1.0d0-1.0d0/uc)*ob(:) !CFŠÔ
            p3(:)=oa(:)+uc*(1.0d0-1.0d0/sc)*oc(:) !AEŠÔ
            p4(:)=ob(:)+uc*(1.0d0-1.0d0/tc)*oc(:) !BFŠÔ
            p5(:)=oa(:)+ob(:)+uc*(1.0d0-1.0d0/sc-1.0d0/tc)*oc(:) !DGŠÔ
            spp(:)=0.2d0*(p1(:)+p2(:)+p3(:)+p4(:)+p5(:))+XL(0,:)
        endif
    elseif (fnn.le.S6.and.fnn.lt.1.0d0) then
        uci=HALF/alpha+sqrt(2.0d0*(1.0d0-fnn)/beta - TWFTH/alpha**2)
        uc=-(uci-1.0d0-1.0d0/alpha-1.0d0/beta)
        sc=uc*alpha
        tc=uc*beta
        if(iposi.eq.1) then
            p1(:)=oa(:)+oc(:)+tc*(1.0d0-1.0d0/uc-1.0d0/sc)*ob(:) !EGŠÔ 
            p2(:)=oc(:)+tc*(1.0d0-1.0d0/uc)*ob(:) !CFŠÔ
            p3(:)=ob(:)+uc*(1.0d0-1.0d0/tc)*oc(:) !BFŠÔ
            p4(:)=oa(:)+ob(:)+uc*(1.0d0-1.0d0/sc-1.0d0/tc)*oc(:) !DGŠÔ
            spp(:)=0.25d0*(p1(:)+p2(:)+p3(:)+p4(:))+XL(0,:)
        endif
    else
        uci=(6.0d0*(1.0d0-fnn)/alpha/beta)**(THIRD)
        uc=-(uci-1.0d0-1.0d0/alpha-1.0d0/beta)
        sc=uc*alpha
        tc=uc*beta
        if(iposi.eq.1) then
            p1(:)=ob(:)+oc(:)+sc*(1.0d0-1.0d0/tc-1.0d0/uc)*oa(:) !FGŠÔ
            p2(:)=oa(:)+oc(:)+tc*(1.0d0-1.0d0/sc-1.0d0/uc)*ob(:) !EGŠÔ
            p3(:)=oa(:)+ob(:)+uc*(1.0d0-1.0d0/sc-1.0d0/tc)*oc(:) !DGŠÔ
            spp(:)=THIRD*(p1(:)+p2(:)+p3(:))+XL(0,:)
        endif
    endif
endif !Full 3D calculation

contains
    
doubleprecision function taiseki2(s,u)

        DoublePrecision,intent(in)::s,u
        !DoublePrecision::ti,si,ui
        if(u.lt.0.0) then !if1
            taiseki2=ZERO
        else if(u.le.1.0d0) then !if1
           if(s.lt.1.0d0) then
             taiseki2=HALF*s*u
           else if(u*(1.0d0-1.0d0/s).lt.1.0d0) then
             taiseki2=u-HALF*u/s
           else 
            continue  
           endif
        else
           if(s.lt.1.0d0) then
             taiseki2=s-HALF*s/u
           else if(u*(1.0d0-1.0d0/s).lt.1.0d0) then
        taiseki2=HALF*u*s*(1.0d0-(1.0d0-1.0d0/s)**2-(1.0d0-1.0d0/u)**2)
           else 
            taiseki2=1.0d0  
           endif
        endif
end function taiseki2    

function nrm_nrm(nrm,pmt1)
doubleprecision,dimension(0:2)::nrm,nrm_nrm
doubleprecision,intent(in)::pmt1
nrm_nrm=nrm
if(abs(nrm(0)).lt.pmt1) nrm_nrm(0)=ZERO
if(abs(nrm(1)).lt.pmt1) nrm_nrm(1)=ZERO
if(abs(nrm(2)).lt.pmt1) nrm_nrm(2)=ZERO
nrm_nrm(0:2)=nrm_nrm(0:2)/sqrt(sum(nrm_nrm(0:2)**2))
continue
end function

subroutine make_XL_ob(iflag)
    integer,intent(in) :: iflag
    integer :: d1, d2, i
    XL(1,0:2) = x0(0:2) !global origin
if (iflag.eq.1) then
    d1 = 0; d2 = 2
    XL(2,0)=dx(0);  XL(2,1)=ZERO;  XL(2,2)=ZERO !A
    XL(3,0)=dx(0);  XL(3,1)=ZERO;  XL(3,2)=dx(2) !E
    XL(4,0)=ZERO;  XL(4,1)=ZERO;  XL(4,2)=dx(2) !C
    XL(5,0)=ZERO;  XL(5,1)=dx(1);  XL(5,2)=ZERO !B       
    XL(6,0)=dx(0);  XL(6,1)=dx(1);  XL(6,2)=ZERO !D
    XL(7,0)=dx(0);  XL(7,1)=dx(1);  XL(7,2)=dx(2) !G
    XL(8,0)=ZERO;  XL(8,1)=dx(1);  XL(8,2)=dx(2) !F
elseif (iflag.eq.0) then
    d1 = 1; d2 = 2
    XL(2,0)=ZERO;  XL(2,1)=dx(1);  XL(2,2)=ZERO !B
    XL(3,0)=ZERO;  XL(3,1)=dx(1);  XL(3,2)=dx(2) !F
    XL(4,0)=ZERO;  XL(4,1)=ZERO;  XL(4,2)=dx(2) !C     
    XL(5,0)=dx(0);  XL(5,1)=ZERO;  XL(5,2)=ZERO !A
    XL(6,0)=dx(0);  XL(6,1)=dx(1);  XL(6,2)=ZERO !D           
    XL(7,0)=dx(0);  XL(7,1)=dx(1);  XL(7,2)=dx(2) !G
    XL(8,0)=dx(0);  XL(8,1)=ZERO;  XL(8,2)=dx(2) !E       
elseif (iflag.eq.2) then
    d1 = 0; d2 = 1
    XL(2,0)=dx(0);  XL(2,1)=ZERO;  XL(2,2)=ZERO !A
    XL(3,0)=dx(0);  XL(3,1)=dx(1);  XL(3,2)=ZERO !D    
    XL(4,0)=ZERO;  XL(4,1)=dx(1);  XL(4,2)=ZERO !B        
    XL(5,0)=ZERO;  XL(5,1)=ZERO;  XL(5,2)=dx(2) !C        
    XL(6,0)=dx(0);  XL(6,1)=ZERO;  XL(6,2)=dx(2) !E
    XL(7,0)=dx(0);  XL(7,1)=dx(1);  XL(7,2)=dx(2) !G
    XL(8,0)=ZERO;  XL(8,1)=dx(1);  XL(8,2)=dx(2) !F   
elseif (iflag.eq.-1) then
    !3D calc
    XL(2,0)=dx(0);  XL(2,1)=ZERO;  XL(2,2)=ZERO
    XL(3,0)=dx(0);  XL(3,1)=dx(1);  XL(3,2)=ZERO
    XL(4,0)=ZERO;  XL(4,1)=dx(1);  XL(4,2)=ZERO
    XL(5,0)=ZERO;  XL(5,1)=ZERO;  XL(5,2)=dx(2)
    XL(6,0)=dx(0);  XL(6,1)=ZERO;  XL(6,2)=dx(2)
    XL(7,0)=dx(0);  XL(7,1)=dx(1);  XL(7,2)=dx(2)
    XL(8,0)=ZERO;  XL(8,1)=dx(1);  XL(8,2)=dx(2)
    !Add back on the global origin..
    do i=2,8
        XL(i,0:2) = XL(i,0:2) + XL(1,0:2)
    enddo
    if(nrm(0).ge.ZERO) then
        if(nrm(1).ge.ZERO.and.nrm(2).ge.ZERO) then
            XL(0,:)=XL(1,:)
        else if(nrm(1).ge.ZERO.and.nrm(2).lt.ZERO) then
            XL(0,:)=XL(5,:)
        else if(nrm(1).lt.ZERO.and.nrm(2).ge.ZERO) then
            XL(0,:)=XL(4,:)
        else if(nrm(1).lt.ZERO.and.nrm(2).lt.ZERO) then
            XL(0,:)=XL(8,:)
        endif
    else
        if(nrm(1).ge.ZERO.and.nrm(2).ge.ZERO) then
            XL(0,:)=XL(2,:)
        else if(nrm(1).ge.ZERO.and.nrm(2).lt.ZERO) then
            XL(0,:)=XL(6,:)
        else if(nrm(1).lt.ZERO.and.nrm(2).ge.ZERO) then
            XL(0,:)=XL(3,:)
        else if(nrm(1).lt.ZERO.and.nrm(2).lt.ZERO) then
            XL(0,:)=XL(7,:)
        endif
    endif
    return
endif
!Add back on the global origin..
do i=2,8
    XL(i,0:2) = XL(i,0:2) + XL(1,0:2)
enddo
if (nrm(d1).gt.ZERO.and.nrm(d2).gt.ZERO) then
    XL(0,:)=XL(1,:)
    ob(:)=XL(5,:)-XL(0,:)
elseif (nrm(d1).gt.ZERO.and.nrm(d2).lt.ZERO) then
    XL(0,:)=XL(4,:)   !C
    ob(:)=XL(8,:)-XL(0,:)
elseif (nrm(d1).lt.ZERO.and.nrm(d2).gt.ZERO) then
    XL(0,:)=XL(2,:)   !A   
    ob(:)=XL(6,:)-XL(0,:)        
else
    XL(0,:)=XL(3,:)   !E 
    ob(:)=XL(7,:)-XL(0,:)  
endif
end subroutine make_XL_ob
!
end subroutine suimen_teimen
    
    

