subroutine PCGPME(AL,NL,NA,LL,D,N,N2,B,X,EPS,PARM,ITER,IER)
!
!************************************************************************
!*  PRECONDITIONED MULTIPLY CONJUGATED GRADIENT POLINOMIAL MINIMIZED    *
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
!*     X      THE INITIAL VALUES OF THE SOLUTION. ZEROS ARE PERMITTED.  *
!*     EPS    THE TOLERLANCE FOR CONVERGENCE.                           *
!*     PARM   PARAMETER FOR CONVERGENCE. ZERO IS PERMITTED.             *
!*     ITR    THE MAXIMUM NUMBER OF ITERATIONS.                         *
!*  ON RETURN:                                                          *
!*     X      THE SOLUTION VECTOR.                                      *
!*     EPS    THE ERROR OF THE SOLUTION ON RETURN.                      *
!*     ITR    THE NUMBER OF ITERATION ON RETURN.                        *
!*  OTHERS:  WORKING PARAMETERS.                                        *
!*                                                                      *
!*  COPYRIGHT:      USHIRO YASUNORI       NOV. 1 1991        VER. 1     *
!************************************************************************

      use variables, only: NTHR, ZERO, TINY, SMALL, ONE, TWO 
      use omp_lib
      implicit none
      integer,intent(in) :: NL,NA,N,N2
      real*8,dimension(NA,NL),intent(inout) :: AL
      real*8,dimension(N),intent(inout) :: D,B
      real*8,dimension(N+N2),intent(inout) :: X
      real*8,dimension(2),intent(inout) :: EPS
      integer,dimension(NA,NL),intent(inout) :: LL
      integer,dimension(NA,NL) :: LLTEMP
      real*8,intent(in) :: PARM
      integer,intent(out) :: IER
      integer,intent(inout) :: ITER
      real*8,dimension(N+N2,5) :: WK
      real*8,dimension(N+N2) :: W
      real*8 :: BN,C1,ALP,AMU,RN,BN0,RQ,BETA,ERR,SS,ERRP
      integer :: K,I,N1,KK,ITEMP
      real*8 :: P1,P2,TEMP

      if (EPS(1) <= TINY) EPS(1) = SMALL
      if (ITER <= 0) ITER = 10 * (dsqrt(dfloat(N))+1)
      IER = 0
      ERRP = ZERO
      N1 = N + 1
      P1 = PARM
      if (PARM < ONE) P1 = TWO * TWO
      P2 = P1 + ONE
!$omp parallel
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N+N2
        W (I)  = ZERO
        WK(I,1)= ZERO
        WK(I,2)= ZERO
        WK(I,3)= ZERO
        WK(I,4)= ZERO
        WK(I,5)= ZERO
      enddo
!$omp end do
! Store original LL for later
!$omp do schedule(static,(NA-1)/NTHR+1)
      do I = 1,NA
          LLTEMP(I,:) = LL(I,:)
      enddo
!$omp enddo
      
!$omp do private(K,TEMP) schedule(static,(N-1)/NTHR+1)
      do I=1,N
        do K=1,NL
          if (LL(I,K)==0) LL(I,K)=N1+mod(I,64)
          W(I) = W(I) + dabs(AL(I,K))
        enddo
        TEMP = D(I)
        W(I) = P2 / ( max(W(I),TEMP) + TEMP*P1 )
!*    AL=AL*W
        do K=1,NL
          AL(I,K) = AL(I,K) * W(I)
        enddo
!*    D=D*W
        TEMP = W(I)
        D(I) = D(I) * TEMP
        B(I) = B(I) * TEMP
      enddo
!$omp enddo
      
!*    SORT.
!$omp do private(K,KK,TEMP,ITEMP) schedule(static,(N-1)/NTHR+1)
      do I=1,N
        do K=1,NL-1
          do KK=K+1,NL
            if ( dabs(AL(I,K)) < dabs(AL(I,KK)) ) then
              TEMP=AL(I,K)
              ITEMP=LL(I,K)
              AL(I,K)=AL(I,KK)
              LL(I,K)=LL(I,KK)
              AL(I,KK)=TEMP
              LL(I,KK)=ITEMP
            endif
          enddo
        enddo
      enddo
!$omp end do
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
      BN=ZERO
      C1=ZERO
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

        RQ=ZERO
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

        ALP = C1/RQ
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
        AMU=ZERO
        SS=ZERO
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
        C1=ZERO
        RN=ZERO
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
        
        ERR = dsqrt(RN/BN)
#ifdef BUGWIN
        if(abs(bn).lt.1.0d-20.or.abs(rn).gt.1.0d+20) then
           write(6,*) 'k=',k,'bn,rn,err=',bn,rn,err
        endif
!
!---write error values
        if (mod(k,100).eq.0) then
          write(6,6000) k ,err,rn,bn
        endif
 6000   format(1h ,'itrno.=',i6,1x,'err=',e17.10,1x,e17.10,1x,e17.10)
#endif
!
        if (K >= 10 .and. abs(ERR-ERRP) <= ERR * 1.0d-6) then
          !write(6,*) k,err,errp,err-errp
          goto 110
        endif
        if (ERR <= EPS(1)) goto 110
        ERRP = ERR
        if (amu == ZERO) then
          write(6,*) 'return due to amu=',amu
          goto 110  !‚±‚êˆÈãŒvŽZ‚µ‚Ä‚à–³—
        endif

!*      BETA=C1/(AMU*RQ).
        if (C1.eq.ZERO) goto 110  !‚±‚êˆÈãŒvŽZ‚µ‚Ä‚à–³—
        BETA = C1/AMU/RQ
!$omp parallel do schedule(static,(N-1)/NTHR+1)
        do I=1,N
            WK(I,1) = B(I) + BETA * (WK(I,1)-AMU*WK(I,3))
        enddo
!$omp end parallel do
      enddo
!**   END OF PMCGPM ITERATION.

      write(6,*) 'bn0=',BN0
      IER=3000

  110 continue
      ITER=K
      EPS(1)=ERR
      
      return


      contains


      subroutine LDUSUB(P,Q)
!*    Q=PM*P
      implicit none
      real*8,dimension(N+N2),intent(in)::P
      real*8,dimension(N+N2),intent(out)::Q
      real*8,dimension(N+N2)::W
      integer::I

!$omp parallel
!$omp do schedule(static,(N+N2-1)/NTHR+1)
      do I=1,N+N2
        W(I) = ZERO
      enddo
!$omp enddo

!*    J=1
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        W(I) = P(I)-AL(I,1)*P(LL(I,1))
      enddo
!$omp enddo

!*    J=2
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        Q(I) = W(I)-AL(I,2)*W(LL(I,2))
      enddo
!$omp enddo

!*    J=3
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        W(I) = Q(I)-AL(I,3)*Q(LL(I,3))
      enddo
!$omp enddo

!*    J=4
!$omp do schedule(static,(N-1)/NTHR+1)
      do I=1,N
        Q(I) = W(I)-AL(I,4)*W(LL(I,4))
      enddo
!$omp enddo

!*    J=5,6
      if (NL>4) then
!$omp do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          W(I) = Q(I)-AL(I,5)*Q(LL(I,5))-AL(I,6)*Q(LL(I,6))
        enddo
!$omp enddo
!$omp do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          Q(I) = W(I)
        enddo
!$omp enddo
      endif
!$omp end parallel
      end subroutine


      subroutine AXSUB(X,Y)
!*    Y=A*X.

      implicit none
      real*8,dimension(N+N2),intent(in) :: X
      real*8,dimension(N+N2),intent(out) :: Y
      integer :: I

      if (NL==6) then   ! Full 3D
!$omp parallel do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          Y(I) = X(I)*D(I)+AL(I,1)*X(LL(I,1))+AL(I,2)*X(LL(I,2)) &
                        +AL(I,3)*X(LL(I,3))+AL(I,4)*X(LL(I,4))   &
                        +AL(I,5)*X(LL(I,5))+AL(I,6)*X(LL(I,6))
        enddo
!$omp end parallel do
      elseif (NL==4) then   ! Vertical 2D
!$omp parallel do schedule(static,(N-1)/NTHR+1)
        do I=1,N
          Y(I) = X(I)*D(I)+AL(I,1)*X(LL(I,1))+AL(I,2)*X(LL(I,2)) &
                        +AL(I,3)*X(LL(I,3))+AL(I,4)*X(LL(I,4))
        enddo
!$omp end parallel do
      else
        stop 'AXSUB in PCGMPME: NL is not equal to 4 or 6'
      endif

      return

      end subroutine

end subroutine PCGPME