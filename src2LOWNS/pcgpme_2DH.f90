SUBROUTINE PCGPME_WILL(AL,NL,NA,LL,D,N,N2,B,EPS,PARM,ITER)
!************************************************************************
!*  PRECONDITIONED MULTIPLY CONJUGATED GRADIENT POLYNOMIAL MINIMIZED    *
!*  METHOD FOR FINITE ELEMENT METHOD.                                   *
!*                                                                      *
!*  PARAMETERS:                                                         *
!*   ON ENTRY:                                                          *
!*     AL     NON-ZERO ELEMENTS OF EACH ROW OF THE MATRIX A EXCLUDING   *
!*            DIAGONAL ELEMENTS.                                        *
!*     NL     MAXIMUNM NUMBER OF NON-ZERO ELEMENTS IN EACH ROW OF A.    *
!*     NA     THE LEADING DIMENSION OF THE ARRAY AL.                    *
!*     LL     COLUMN INDEX OF NON-ZERO ELEMENTS OF EACH ROW OF THE      *
!*            MATRIX A EXCLUDING DISGONAL ELEMENTS.                     *
!*     D      DIAGONAL ELEMENTS OF THE MATRIX A.                        *
!*     N      THE ORDER OF THE MATRIX A.                                *
!*     N2     N+N2 IS THE LEADING DIMENSION OF X, WK, DD. USUALLY N2=62 *
!*     B      THE RIGHT HAND SIDE OF THE EQUATIONS.                     *
!*     EPS    THE TOLERLANCE FOR CONVERGENCE.                           *
!*     PARM   PARAMETER FOR CONVERGENCE. ZERO IS PERMITTED.             *
!*     ITER   THE MAXIMUM NUMBER OF ITERATIONS.                         *
!*  ON RETURN:                                                          *
!*     B      THE SOLUTION VECTOR.                                      *
!*  OTHERS:  WORKING PARAMETERS.                                        *
!*                                                                      *
!*  COPYRIGHT:      USHIRO YASUNORI       NOV. 1 1991        VER. 1     *
!************************************************************************
!$    use omp_lib
      use variables, only: ONE, ZERO, TWO
      implicit none 
      integer,intent(in)                  :: N, N2, NA, NL
      real*8,intent(in)                   :: PARM, EPS
      integer,intent(in)                  :: ITER
      real*8,dimension(NA,NL),intent(in)  :: AL
      real*8,dimension(NA,NL)             :: ALTEMP
      real*8,dimension(N),intent(in)      :: D
      real*8,dimension(N),intent(inout)   :: B
      real*8,dimension(N)                 :: DTEMP
      integer,dimension(NA,NL),intent(in) :: LL
      integer,dimension(NA,NL)            :: LLTEMP
      real*8,dimension(N+N2)              :: X, W
      real*8,dimension(N+N2,5)            :: WK
      real*8 :: BN, C1, ALP, AMU, RN, BN0, RQ, BETA, ERR, SS, ERRP
      real*8 :: P1, P2, TEMP, TRUNC
      integer:: K, I, NTHR, IER
      integer:: N1, KK, ITEMP
      real*8,parameter :: L_Limit = 1d-14
      
      IER = 0
      ERRP = ZERO
      N1 = N + 1
      P1 = PARM
      if (PARM < ONE) P1 = TWO * TWO
      P2 = P1 + ONE
!$omp parallel
!$omp single
!$ NTHR = omp_get_num_threads()
!$omp end single
!$omp do schedule(static,(N+N2-1)/NTHR+1)
      do I=1,N+N2
        X(i)    = ZERO
        W(I)    = ZERO
        WK(I,1) = ZERO
        WK(I,2) = ZERO 
        WK(I,3) = ZERO
        WK(I,4) = ZERO
        WK(I,5) = ZERO
        if (i.gt.n) cycle
        DTEMP(i)= D(i)
      enddo
!$omp end do

!$omp do schedule(static,(NA-1)/NTHR+1)
      do i = 1,NA
          ALTEMP(i,:) = AL(i,:)
          LLTEMP(i,:) = LL(i,:)
      enddo
!$omp enddo

! ========== DECOMPOSITION ==========
!*    AL=AL+AU.

!$omp do private(K,TEMP) schedule(static,(N-1)/NTHR+1)
      do I=1,N
        do K=1,NL
          if(LLTEMP(I,K)==0) LLTEMP(I,K)=N1+mod(I,64)
          W(I)=W(I)+dabs(ALTEMP(I,K))
        enddo
        TEMP=DTEMP(I)
        W(I)=P2/(max(W(I),TEMP)+TEMP*P1)
!*    AL=AL*W
        do K=1,NL
          ALTEMP(I,K)=ALTEMP(I,K)*W(I)
        enddo
!*    D=D*W
        TEMP=W(I)
        DTEMP(I)=DTEMP(I)*TEMP
        B(I)=B(I)*TEMP
      enddo
!$omp enddo
      
!*    SORT.
!$omp single
      do K=1,NL-1
        do KK=K+1,NL
          do I=1,N
            if( dabs(ALTEMP(I,K)) < dabs(ALTEMP(I,KK)) ) then
!              call swap(ALTEMP(I,K),ALTEMP(I,KK))
!              call swap(LLTEMP(I,K),LLTEMP(I,KK))
              TEMP=ALTEMP(I,K)
              ITEMP=LLTEMP(I,K)
              ALTEMP(I,K)=ALTEMP(I,KK)
              LLTEMP(I,K)=LLTEMP(I,KK)
              ALTEMP(I,KK)=TEMP
              LLTEMP(I,KK)=ITEMP
            endif
          enddo
        enddo
      enddo
!$omp end single
! =========================================

!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        WK(I,1)=B(I)
      enddo
!$omp end do
!$omp end parallel

!*    B=PM*B.
      CALL LDUSUB(WK,WK)

!*    Q=PM*A*X.
      CALL AXSUB(X,WK(1,3))
      CALL LDUSUB(WK(1,3),WK(1,3))

!*    BN=(B,B) , P=R0=R=B-Q , C1=(R,R).
      BN=0.0d0
      C1=0.0d0
!$omp parallel
!$omp do private(TEMP) reduction(+:BN,C1) schedule(static,(N-1)/NTHR+1)
      do I=1,N
        TEMP=WK(I,1)
        BN=BN+TEMP*TEMP
        TEMP=TEMP-WK(I,3)
        B(I)=TEMP
        C1=C1+TEMP*TEMP
        WK(I,1)=TEMP
        WK(I,2)=TEMP
      enddo
!$omp enddo
!$omp end parallel
      BN0=BN
      
!**   PMCGPM ITERATION.
      do K=1,ITER
!*    Q=PM*A*P.
        CALL AXSUB(WK(1,1),WK(1,3))
        CALL LDUSUB(WK(1,3),WK(1,3))
        
!*    ALP=C1/(R0,Q).

        RQ=0.0d0
!$omp parallel
!$omp do reduction(+:RQ) schedule(static,(N-1)/NTHR+1)
        do I=1,N
          RQ=RQ+WK(I,2)*WK(I,3)
        enddo
!$omp enddo
!$omp single
        if(isnan(RQ)) then
          write(6,*) 'RQ=',RQ
        endif

        ALP=C1/RQ
!$omp end single
        
!*    E=R-ALP*Q.
!$omp do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          WK(I,4)=B(I)-ALP*WK(I,3)
        enddo
!$omp enddo
!$omp end parallel
        
!*    V=PM*A*E.
        CALL AXSUB(WK(1,4),WK(1,5))
        CALL LDUSUB(WK(1,5),WK(1,5))
        
!*    AMU=(E,V)/(V,V).
        AMU = ZERO
        SS  = ZERO
!$omp parallel private(TEMP)
!$omp do reduction(+:AMU,SS) schedule(static,(N-1)/NTHR+1)
        do I=1,N
          TEMP=WK(I,5)
          AMU=AMU+WK(I,4)*TEMP
          SS=SS+TEMP*TEMP
        enddo
!$omp enddo
!$omp single
        if(isnan(AMU)) then            
          write(6,*) 'amu0=',amu,ss,'err=',err,k
        endif
        AMU=AMU/SS

!*    X=X+ALP*P+AMU*E , R=E-AMU*V , RN=(R,R) , C1=(R0,R).
        C1 = ZERO
        RN = ZERO
!$omp end single
!$omp do reduction(+:RN,C1) schedule(static,(N-1)/NTHR+1)
        do I=1,N
          TEMP=WK(I,4)
          X(I)=X(I)+ALP*WK(I,1)+AMU*TEMP
          TEMP=TEMP-AMU*WK(I,5)
          B(I)=TEMP
          RN=RN+TEMP*TEMP
          C1=C1+TEMP*WK(I,2)
        enddo
!$omp enddo
!$omp end parallel

        ERR=DSQRT(RN/BN)
        
#ifdef BUG
       if(abs(bn).lt.1.0d-10) write(6,*) 'bn=',bn,'rn=',rn,err
     
        if(abs(bn).lt.1.0d-20.or.abs(rn).gt.1.0d+20) then
           write(6,*) 'k=',k,'bn,rn,err=',bn,rn,err
           stop
        endif
#endif
#ifdef BUGWIN
        if(abs(bn).lt.1.0d-20.or.abs(rn).gt.1.0d+20) then
           write(6,*) 'k=',k,'bn,rn,err=',bn,rn,err
        endif

!
!---write error values
        if (mod(k,100).eq.0) then
          write(6,6000) k ,err,rn,bn
          call flush(6)
        endif
 6000   format(1h ,'itrno.=',i6,1x,'err=',e17.10,1x,e17.10,1x,e17.10)
#endif
!
        IF (K.ge.10.and.abs(ERR-ERRP).LE.ERR*1.0d-6) then
!                write(6,*) k,err,errp,err-errp
          GO TO 110
        endif
        IF (ERR.LE.EPS) GO TO 110
        ERRP = ERR
        IF (amu.eq.ZERO) GO TO 110  !‚±‚êˆÈãŒvŽZ‚µ‚Ä‚à–³—
        
        if (C1.eq.ZERO) goto 110  !‚±‚êˆÈãŒvŽZ‚µ‚Ä‚à–³—
        BETA = C1 / ( AMU * RQ ) 
!$omp parallel do schedule(static,(N-1)/NTHR+1)
        do I = 1,N
          WK(I,1) = B(I) + BETA * ( WK(I,1) - AMU * WK(I,3) )
        enddo
!$omp end parallel do
      enddo
      write(6,*) 'bn0=',bn0,'err',err
      IER = 3000
      
110   continue
      ! Write out the solution vector, suppress small values,
      ! & search surrounding nodes for small numerical errors
      do i = 1,n
         if (abs(X(i)).le.L_Limit) X(i) = ZERO
         do k = 1,NL
            if (LL(i,k).eq.i) cycle
            if (LL(i,k).eq.0) exit
            if (X(i).eq.X(LL(i,k))) cycle
            if (abs(X(i)-X(LL(i,k))).le.L_Limit) then
                X(i) = X(LL(i,k))
            endif
            B(i) = X(i)
         enddo
      enddo
      RETURN
      
      contains
!
      subroutine LDUSUB(P,Q)
!*    Q=PM*P
      implicit none
      real*8,dimension(N+N2),intent(in)  :: P
      real*8,dimension(N+N2),intent(out) :: Q
      real*8,dimension(N+N2)             :: W
      integer :: I

!$omp parallel
!$omp do schedule(static,(N+N2-1)/NTHR+1)
      do I=1,N+N2
        W(I) = ZERO
      enddo
!$omp enddo

!*    J=1
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        W(I)=P(I)-ALTEMP(I,1)*P(LLTEMP(I,1))
      enddo
!$omp enddo

!*    J=2
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        Q(I)=W(I)-ALTEMP(I,2)*W(LLTEMP(I,2))
      enddo
!$omp enddo

!*    J=3
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        W(I)=Q(I)-ALTEMP(I,3)*Q(LLTEMP(I,3))
      enddo
!$omp enddo

!*    J=4
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        Q(I)=W(I)-ALTEMP(I,4)*W(LLTEMP(I,4))
      enddo
!$omp enddo

!*    J=5,6
      if(NL>4) then
!$omp do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          W(I)=Q(I)-ALTEMP(I,5)*Q(LLTEMP(I,5))-ALTEMP(I,6)*Q(LLTEMP(I,6))
        enddo
!$omp enddo
!$omp do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          Q(I)=W(I)
        enddo
!$omp enddo
      endif
!$omp end parallel
      end subroutine


      subroutine AXSUB(X,Y)
!*    Y=A*X.
      implicit none
      real*8,dimension(N+N2),intent(in)  :: X
      real*8,dimension(N+N2),intent(out) :: Y
      integer :: I

      if (NL==6) then   ! Full 3D
!$omp parallel do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          Y(I)=X(I)*DTEMP(I)+ALTEMP(I,1)*X(LLTEMP(I,1))+ALTEMP(I,2)*X(LLTEMP(I,2)) &
              +ALTEMP(I,3)*X(LLTEMP(I,3))+ALTEMP(I,4)*X(LLTEMP(I,4))               &
              +ALTEMP(I,5)*X(LLTEMP(I,5))+ALTEMP(I,6)*X(LLTEMP(I,6))
        enddo
!$omp end parallel do
      elseif (NL==4) then   ! Vertical 2D
!$omp parallel do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          Y(I)=X(I)*DTEMP(I)+ALTEMP(I,1)*X(LLTEMP(I,1))+ALTEMP(I,2)*X(LLTEMP(I,2)) &
              +ALTEMP(I,3)*X(LLTEMP(I,3))+ALTEMP(I,4)*X(LLTEMP(I,4))
        enddo
!$omp end parallel do
      else
        stop 'AXSUB in PCGMPME: NL is not equal to 4 or 6'
      endif

      return

      end subroutine
!
      end subroutine PCGPME_WILL