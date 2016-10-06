!%%%%%%%%%%%%%%%%%%%% FILE: 2Ddata_io.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads input data and control file
subroutine read_control_file2D
    use variables, only: LNUM, tide_level, input2D, inputxyz2D, inputf2D, IM
    use arrays, only: L
    implicit none
    integer :: LN
    !
    if (IM.eq.3) return ! No 2DH calc.
    read(5,*) tide_level ; print *,'tide level [m] = ', tide_level
    !Read number of Layers, allocate and read each layer information
    read(5,*) LNUM ; allocate(L(LNUM)) ; print *,'LNUM = ', LNUM
    do LN = 1,LNUM
        read(5,*) L(LN)%dt,L(LN)%cl_set
        read(5,*) L(LN)%SWE,L(LN)%DISP,L(LN)%GEOC,L(LN)%BC
        read(5,*) L(LN)%Dmin,L(LN)%CD,L(LN)%IC
        read(5,*) L(LN)%Mandef
        if (IM.eq.1.or.LN.lt.LNUM) then
            read(5,*) L(LN)%coupling_dim, L(LN)%write_out
        else
            L(LN)%coupling_dim = 1 !(No further inner layer)
            L(LN)%write_out = 1    !(output by default)
        endif
    enddo
    !Read the 2D input files
    read(5,'(A255)') input2D    ; print *,'input2D=',trim(input2D)
    read(5,'(A255)') inputxyz2D ; print *,'input_xyz2D=',trim(inputxyz2D)
    read(5,'(A255)') inputf2D   ; print *,'inputf2D=',trim(inputf2D)
end subroutine read_control_file2D
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                          SUBROUTINE: READ_DATA2D                        &&!
!&&    THIS SUBROUTINE READS THE DATA FROM THE INPUT FILES AND              &&!
!&&    CONVERTS THEM INTO THE PARAMETERS USED IN THE CALCULATION            &&!
!&&                                                        - PRINGLE (2015) &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! 
subroutine read_data2D
    use variables, only: LNUM, t, fg_xyz, inputxyz2D, inputf2D, input2D,      &
                         ZERO, ONE
    use arrays, only: L, alloc_XY, alloc_F
    implicit none
    integer :: LN, nn, inn, i, j, k, na1, na2, icc, icc2
    integer,allocatable :: intemp(:,:)
    real*8 :: dummy, gomi
    character*1 :: c1
    character*2 :: c2
    logical :: exist
!
!----------- Read the XY data -------------------------------------------------
na1 = len_trim(inputxyz2D) ; na2 = len_trim(inputf2D)
write(6,*) '[HOST]read STdata',inputxyz2D(1:na1)
if (inputxyz2D(na1:na1).eq.'n') then
    fg_xyz = 1 !New way
else 
    fg_xyz = 0 !Old way
endif
!
!Loop over all the layers 
DO LN = 1,LNUM
    if (LNUM.eq.1) then
        !Just open the normal file name
        OPEN(8,FILE=inputxyz2D(1:na1),status='old',form='unformatted')
    else
        if (LN.lt.10) then
            write(c1,'(I1)') LN
            inquire(file = inputxyz2D(1:na1-4-fg_xyz)//'_'//c1//              &
                           inputxyz2D(na1-3-fg_xyz:na1), exist = exist )       
        elseif (LN.ge.10) then
            write(c2,'(I2)') LN
            inquire(file = inputxyz2D(1:na1-4-fg_xyz)//'_'//c2//              &
                           inputxyz2D(na1-3-fg_xyz:na1), exist = exist )
        endif
        if (exist) then
            if (LN.lt.10) OPEN(8,FILE = inputxyz2D(1:na1-4-fg_xyz)//'_'//c1// &
                                        inputxyz2D(na1-3-fg_xyz:na1),         &
                                        status = 'old',form = 'unformatted' )
            if (LN.ge.10) OPEN(8,FILE = inputxyz2D(1:na1-4-fg_xyz)//'_'//c2// &
                                        inputxyz2D(na1-3-fg_xyz:na1),         &
                                        status = 'old',form = 'unformatted' )
        else  !Read from F data name    
            if (LN.lt.10) OPEN(8,FILE = inputf2D(1:na2-6)//'_'//c1//          &
                                        inputxyz2D(na1-3-fg_xyz:na1),         &
                                        status = 'old',form = 'unformatted' ) 
            if (LN.ge.10) OPEN(8,FILE = inputf2D(1:na2-6)//'_'//c2//          &
                                        inputxyz2D(na1-3-fg_xyz:na1),         &
                                        status = 'old',form = 'unformatted' )
        endif
    endif
    !Read bounds of the layer
    read(8) L(LN)%is,L(LN)%js,L(LN)%ks,L(LN)%ie,L(LN)%je,L(LN)%ke
    write(6,*) 'Layer No.', LN
    write(6,*) L(LN)%is,L(LN)%js,L(LN)%ks,L(LN)%ie,L(LN)%je,L(LN)%ke
    if (L(LN)%js.eq.L(LN)%je-1) then
        L(LN)%DIM = 1
    else
        L(LN)%DIM = 2
    endif
    write(6,*) 'Dimension of Leapfrog calculation:',L(LN)%DIM
    !Allocate the arrays in the layer
    call alloc_XY(LN,L(LN)%ie+1,L(LN)%je+1,L(LN)%ks,L(LN)%ke)
    !Read the x values
    read(8) (L(LN)%X(i),i=L(LN)%is-2,L(LN)%ie+1)
    !Read the y values
    read(8) (L(LN)%Y(j),j=L(LN)%js-2,L(LN)%je+1)
    !Read the z values
    read(8) (L(LN)%Z(k),k=L(LN)%ks-2,L(LN)%ke+1)
    !Read the extra data relating to where we have or don't have fluid
    if (fg_xyz.eq.1) then
        read(8) icc
        !Initialise nf and nfb
        L(LN)%NF = 0 ; L(LN)%NFB = 0
        do nn = 1,icc
            read(8) i,k,j,L(LN)%NF(i,k,j),L(LN)%NFB(i,k,j)
        enddo
    endif
    CLOSE(8)
ENDDO
!
!----------- Read the MAIN data -----------------------------------------------
na1 = len_trim(input2D)
write(6,*) '[HOST]read VPFdata',input2D(1:na1)
na2 = len_trim(inputf2D)
!
!Loop over all the layers 
DO LN = 1,LNUM
    if (LNUM.eq.1) then
        !Just open the normal file name
        OPEN(8,FILE=input2D(1:na1),status='old',form='unformatted')
    else
        !Open the file name with layer number attached
        if (LN.lt.10) then
            write(c1,'(I1)') LN
            inquire(file=input2D(1:na1)//'_'//c1,exist=exist)       
        elseif (LN.ge.10) then
            write(c2,'(I2)') LN
            inquire(file=input2D(1:na1)//'_'//c2,exist=exist)
        endif
        if (exist) then
            if (LN.lt.10) OPEN(8,FILE = input2D(1:na1)//'_'//c1,              &
                                        status = 'old',form = 'unformatted' )
            if (LN.ge.10) OPEN(8,FILE = input2D(1:na1)//'_'//c2,              &
                                        status = 'old',form = 'unformatted' )
        else  !Read from F data name    
            if (LN.lt.10) OPEN(8,FILE = inputf2D(1:na2-6)//'_'//c1,           &
                                        status = 'old',form = 'unformatted' )
            if (LN.ge.10) OPEN(8,FILE = inputf2D(1:na2-6)//'_'//c2,           &
                                        status = 'old',form = 'unformatted' )
        endif
    endif
    read(8) t
    if (fg_xyz.eq.0) then
        !Read NF and NFB data
        read(8) (((L(LN)%NF(I,k,j),K=L(LN)%ks-1,L(LN)%KE),                    &
                                  J=L(LN)%js-1,L(LN)%je),I=L(LN)%is-1,L(LN)%IE)
        read(8) (((L(LN)%NFB(I,k,j),K=L(LN)%ks-1,L(LN)%KE),                   &
                                  J=L(LN)%js-1,L(LN)%je),I=L(LN)%is-1,L(LN)%IE)
         ! Output xyzn data if required
         if (t.eq.ZERO) then
          if (LN.lt.10) then
              inquire(file=input2D(1:na1)//'_'//c1//'.xyzn',exist=exist)      
          else
              inquire(file=input2D(1:na1)//'_'//c2//'.xyzn',exist=exist)      
          endif
          if (.not.exist) then
            if (LN.lt.10) OPEN(9,FILE = input2D(1:na1)//'_'//c1//'.xyzn',     &
                                      status = 'unknown',form = 'unformatted' )
            if (LN.ge.10) OPEN(9,FILE = input2D(1:na1)//'_'//c2//'.xyzn',     &
                                      status = 'unknown',form = 'unformatted' )
                write(9) L(LN)%is,L(LN)%js,L(LN)%ks,L(LN)%ie,L(LN)%je,L(LN)%ke
                write(9) (L(LN)%X(i),i=L(LN)%is-2,L(LN)%ie+1)
                write(9) (L(LN)%Y(j),j=L(LN)%js-2,L(LN)%je+1)
                write(9) (L(LN)%Z(k),k=L(LN)%ks-2,L(LN)%ke+1)
                icc = 0
                do i=L(LN)%is-1,L(LN)%ie
                    do j=L(LN)%js-1,L(LN)%je
                        do k=L(LN)%ks-1,L(LN)%ke
                            if (L(LN)%nf(i,k,j).eq.-1) icc = icc + 1
                        enddo
                    enddo
                enddo 
                write(9) icc    
                write(6,*) 'num of nf=',icc
                do i=L(LN)%is-1,L(LN)%ie
                    do j=L(LN)%js-1,L(LN)%je
                        do k=L(LN)%ks-1,L(LN)%ke
                            if (L(LN)%nf(i,k,j).eq.-1) then
                                write(9) i,k,j,L(LN)%nf(i,k,j),L(LN)%nfb(i,k,j)
                            endif
                        enddo
                    enddo
                enddo
            close(9)
          endif
        endif
    endif
    !
    allocate(intemp(2,L(LN)%ie*L(LN)%je))
    L(LN)%inne=0;L(LN)%mn=0
    !Find the values of arrays to convert 1D array to 2D and vice versa
    do j = L(LN)%js-1,L(LN)%je
        do i = L(LN)%is-1,L(LN)%ie
           k = L(LN)%ks
           !
           if (L(LN)%nf(i,k,j).ge.0.or.(L(LN)%nf(i,k,j).eq.-1                 &
               .and.L(LN)%nfb(i,k,j).ge.1.and.L(LN)%nfb(i,k,j).lt.200)) then
                L(LN)%inne=L(LN)%inne+1 
                L(LN)%mn(i,j)=L(LN)%inne
                intemp(:,L(LN)%inne)= [ i, j ]
           endif  
           if (L(LN)%nf(i,k,j).eq.-1.and.L(LN)%nfb(i,k,j).ge.200) then
                L(LN)%nfb(i,k,j) = 0
           endif  
        enddo
    enddo
    !
    write(6,*) 'Layer No.', LN
    write(6,*) 'Size i,j=',maxval(intemp(1,:)),maxval(intemp(2,:))
    write(6,*) 'Size inne=',L(LN)%inne
    !
    call alloc_F(LN,L(LN)%DIM,L(LN)%inne)
    allocate(L(LN)%IN(2,0:L(LN)%inne)) ; L(LN)%IN(:,0) = 0
    L(LN)%in(:,1:L(LN)%inne) = intemp(:,1:L(LN)%inne)
    deallocate(intemp)

    if (fg_xyz.eq.0) then
        !Read in old data format
        do nn = 1,L(LN)%inne
            !Read the VPF data
            read(8) L(LN)%U(1,nn),dummy,gomi,L(LN)%F(nn),L(LN)%ETAn(nn)
            if (L(LN)%DIM.eq.2) L(LN)%U(2,nn) = dummy
        enddo
    else
        L(LN)%F = ZERO; L(LN)%ETAn = ZERO; L(LN)%U = ZERO
        !Read in new data format
        read(8) icc,icc2
        do inn = 1,icc
            read(8) nn,L(LN)%U(1,nn),dummy,gomi,L(LN)%F(nn),L(LN)%ETAn(nn),   &
                    L(LN)%NF(L(LN)%in(1,nn),L(LN)%ks,L(LN)%in(2,nn)),         &
                    L(LN)%NFB(L(LN)%in(1,nn),L(LN)%ks,L(LN)%in(2,nn))
            if (L(LN)%DIM.eq.2) L(LN)%U(2,nn) = dummy
        enddo
    endif
    !
    if (t.gt.ZERO.and.L(LN)%dt.eq.ZERO) then
        read(8) L(LN)%dt
        write(6,*) 'read dt=',L(LN)%dt   
    endif
    CLOSE(8)
ENDDO 
write(6,*) '[HOST]read end VPFdata',input2D(1:na1)
!
!----------- Read the F data --------------------------------------------------
write(6,*) '[HOST]read fdata',inputf2D(1:na2)
!
!Loop over all the layers 
DO LN = 1,LNUM
    L(LN)%H = ZERO
    if (LNUM.eq.1) then
        !Just open the normal file name
        OPEN(8,FILE=inputf2D(1:na2),status='old',form='unformatted')  
    else
        !Open the file name with layer number attached
        if (LN.lt.10) then
            write(c1,'(I1)') LN
            OPEN(8,FILE=inputf2D(1:na2-6)//'_'//c1//'.fdata',                 &
                 status='old',form='unformatted')  
        elseif (LN.ge.10) then
            write(c2,'(I2)') LN
            OPEN(8,FILE=inputf2D(1:na2-6)//'_'//c2//'.fdata',                 &
                 status='old',form='unformatted') 
        endif
    endif
    if (fg_xyz.eq.1) then
        ! Adopting New format for fdata
        L(LN)%H = ONE ; L(LN)%ZK = ONE ; L(LN)%MAN = ONE !Default vals.
        L(LN)%ZK(0) = L(LN)%Z(L(LN)%ke) ; L(LN)%H(:,0) = -L(LN)%Z(L(LN)%ke)
        read(8) icc,icc2
        do inn = 1,icc
            read(8) nn,L(LN)%H(1,nn),dummy,L(LN)%MAN(nn),L(LN)%ZK(nn) 
            if (L(LN)%DIM.eq.2) L(LN)%H(2,nn) = dummy
        enddo
    else
        ! Old format
        do nn = 1,L(LN)%inne
            read(8) L(LN)%H(1,nn),dummy,L(LN)%MAN(nn),L(LN)%ZK(nn) 
            if (L(LN)%DIM.eq.2) L(LN)%H(2,nn) = dummy
        enddo
    endif
    CLOSE(8)
ENDDO
write(6,*) '[HOST]read end  fdata',inputf2D(1:na2)
!
!Initialising max data
do LN = 1,LNUM
    L(LN)%ETAmax = -1d4 ; L(LN)%Umax = -1d4
enddo
!------------------- Reading MAX data if available ----------------------------
if (t.gt.ZERO) then
    !Max data
    !Loop over all the layers 
    DO LN = 1,LNUM  
        if (LNUM.eq.1) then
            !Just open the normal file name
            do i = 6,9
                inquire(file=input2D(1:na1-i)//'tempMAX.DAT',exist=exist)
                if (exist) then
                    OPEN(8,FILE=input2D(1:na1-i)//'tempMAX.DAT',status='old')
                    exit
                endif
            enddo
        else
            !Open the file name with layer number attached
            do i = 6,9
                if (LN.lt.10) then
                    write(c1,'(I1)') LN
                    inquire(file=input2D(1:na1-i)//'tempMAX_'//c1//'.dat',    &
                            exist=exist)
                elseif (LN.ge.10) then
                    write(c2,'(I2)') LN
                    inquire(file=input2D(1:na1-i)//'tempMAX_'//c2//'.dat',    &
                            exist=exist)
                endif
                if (exist) then
                    if (LN.lt.10) OPEN(8,FILE=input2D(1:na1-i)//'tempMAX_'//  &
                                      c1//'.dat',status='old') 
                    if (LN.ge.10) OPEN(8,FILE=input2D(1:na1-i)//'tempMAX_'//  &
                                      c2//'.dat',status='old')
                    exit
                endif
            enddo
        endif
        if (exist) then
            ! For DAT format
            read(8,*) ! Skip reading the 
            read(8,*) ! x and y vectors
            do while (.true.) ! Read the remaining information
                read(8,*,END=123) i,j,L(LN)%ETAmax(L(LN)%mn(i,j)),gomi,       &
                                  L(LN)%Umax(1,L(LN)%mn(i,j)),dummy
                if (L(LN)%DIM.eq.2) L(LN)%Umax(2,L(LN)%mn(i,j)) = dummy
            enddo
123         CLOSE(8)
        else
            write(6,*) 'Warning: Could not find max data, skipping. LN=',LN
        endif
    ENDDO
    if (exist) write(6,*) '[HOST]read end maxdata (DAT)',                     &
                          input2D(1:na1-i)//'tempMAX.DAT'
endif
!
!============ END OF READING DATA =============================================
endsubroutine read_data2D
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&                         SUBROUTINE: OUTPUT2D                            &&!
!&&    THIS SUBROUTINE OUTPUTS THE CURRENT,MAX OR MCF DATA INTO             &&!
!&&    FILES FOR VISUALIZATION                             - PRINGLE (2015) &&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
subroutine output2D(cflag)
    use arrays, only: L, locout
    use variables, only: LNUM, fg_xyz, ZERO, t, g, nu, locnum, title2D
    implicit none
    character*3,intent(in) :: cflag
    integer :: na1, nn, i, k , j, LN, u, icc, icc2
    real*8 :: dummy
    character*1 :: c1
    character*2 :: c2
    na1 = len_trim(title2D)
!Creating the necessary file based on cflag
!$omp parallel do private(u,nn,i,k,j,c1,c2,icc,dummy)
DO LN = 1,LNUM
    if (L(LN)%WRITE_OUT == 0) cycle !Do not output
    if (cflag.eq.'NOW') then
        if (LNUM.eq.1) then
            OPEN(newunit=u,FILE = title2D(1:na1),                             &
                 form='unformatted', status= 'unknown')
        else
            if (LN.lt.10) then
                write(c1,'(I1)') LN
                OPEN(newunit=u,FILE=title2D(1:na1)//'_'//c1,                  &
                     form='unformatted', status= 'unknown')
            elseif (LN.ge.10) then
                write(c2,'(I2)') LN
                OPEN(newunit=u,FILE=title2D(1:na1)//'_'//c2,                  &
                     form='unformatted', status= 'unknown')
            endif
        endif
    elseif (cflag.eq.'MAX') then
        if (LNUM.eq.1) then
            OPEN(newunit=u,FILE=title2D(1:na1)//'MAX.DAT', status= 'unknown')
        else
            if (LN.lt.10) then
                write(c1,'(I1)') LN
                OPEN(newunit=u,FILE=title2D(1:na1)//'MAX_'//c1//'.dat',       &
                     status= 'unknown')
            elseif (LN.ge.10) then
                write(c2,'(I2)') LN
                OPEN(newunit=u,FILE=title2D(1:na1)//'MAX_'//c2//'.dat',       &
                     status= 'unknown')
            endif
        endif
    endif
    if (cflag.eq.'NOW') then
        write(u) T,g,nu
        if (fg_xyz.eq.0) then
            !Old data format
            write(u) (((L(LN)%NF(I,k,j),K=L(LN)%ks-1,L(LN)%KE),               &
                        J=L(LN)%js-1,L(LN)%je),I=L(LN)%is-1,L(LN)%IE)
            write(u) (((L(LN)%NFB(I,k,j),K=L(LN)%ks-1,L(LN)%KE),              &
                      J=L(LN)%js-1,L(LN)%je),I=L(LN)%is-1,L(LN)%IE)
        else
            !New data format
            icc = 0
            do nn = 1,L(LN)%inne
                if (L(LN)%NF(L(LN)%in(1,nn),L(LN)%ks,L(LN)%in(2,nn)).eq.0.and.&
                   sum(L(LN)%U(:,nn)**2).eq.ZERO.and.L(LN)%F(nn).eq.ZERO) cycle
                icc =  icc + 1
            enddo
            write(u) icc,0 
        endif
        if (fg_xyz.eq.0) then
            !Old data format
            do nn = 1,L(LN)%inne     
                if (L(LN)%DIM.eq.1) dummy = ZERO
                if (L(LN)%DIM.eq.2) dummy = L(LN)%U(2,nn)
                write(u) L(LN)%U(1,nn),dummy,ZERO,L(LN)%F(nn),L(LN)%ETAn(nn)
            enddo
        else
            !New data format
            do nn = 1,L(LN)%inne
                if (L(LN)%NF(L(LN)%in(1,nn),L(LN)%ks,L(LN)%in(2,nn)).eq.0.and.&
                   sum(L(LN)%U(:,nn)**2).eq.ZERO.and.L(LN)%F(nn).eq.ZERO) cycle
                if (L(LN)%DIM.eq.1) dummy = ZERO
                if (L(LN)%DIM.eq.2) dummy = L(LN)%U(2,nn)
               write(u) nn,L(LN)%U(1,nn),dummy,ZERO,L(LN)%F(nn),L(LN)%ETAn(nn)&   
                        ,L(LN)%NF(L(LN)%IN(1,nn),L(LN)%ks,L(LN)%IN(2,nn)),    & 
                         L(LN)%NFB(L(LN)%IN(1,nn),L(LN)%ks,L(LN)%IN(2,nn))
            enddo
        endif
        write(u) L(LN)%DT
    elseif (cflag.eq.'MAX') then
        ! For DAT format
        write(u,'(9999(F13.3))') (L(LN)%X(i),i = L(LN)%is,L(LN)%ie)
        write(u,'(9999(F13.3))') (L(LN)%Y(j),j = L(LN)%js,L(LN)%je)
        do nn = 1,L(LN)%inne
            i = L(LN)%in(1,nn) ; j = L(LN)%in(2,nn)  
            if (L(LN)%NF(i,L(LN)%ks,j).eq.-1) cycle
            if (L(LN)%ETAmax(nn)-L(LN)%ZK(nn).le.L(LN)%Dmin) cycle
            if (L(LN)%DIM.eq.1) dummy = ZERO
            if (L(LN)%DIM.eq.2) dummy = L(LN)%Umax(2,nn)
            write(u,78) i,j,L(LN)%ETAmax(nn),L(LN)%ETAmax(nn)-L(LN)%ZK(nn),   &
                        L(LN)%Umax(1,nn),dummy
        enddo
    endif
    CLOSE(u)
ENDDO
!$omp end parallel do
if (cflag.eq.'NOW') write(6,*) 'write end =',title2D(1:na1)
if (cflag.eq.'MAX') write(6,*) 'write end =',title2D(1:na1)//'MAX.dat'
78 FORMAT (2I5,4E15.7)
end subroutine output2D   

subroutine outputF2D(name)
    use arrays, only: L 
    use variables, only: LNUM, inputf2D, ZERO, ONE, fg_xyz
    implicit none
    character*255,intent(in) :: name
    integer :: na, nn, LN, i, k, j, u, icc
    character*1 :: c1
    character*2 :: c2
    real*8 :: temp1, temp2, temp3
    logical :: exist
    !
    na = len_trim(name)
!$omp parallel do private(u,nn,i,k,j,c1,c2,icc,temp1,temp2,temp3)
    DO LN = 1,LNUM
        if (L(LN)%WRITE_OUT == 0) cycle !Do not output
        if (LNUM.eq.1) then
            inquire(FILE=name(1:na)//'.fdata',exist = exist)
            if (exist) cycle
            OPEN(newunit=u,FILE=name(1:na)//'.fdata',                         &
                 form='unformatted', status= 'unknown')
        else
            if (LN.lt.10) then
                write(c1,'(I1)') LN
                inquire(FILE=name(1:na)//'_'//c1//'.fdata',exist = exist)
                if (exist) cycle
                OPEN(newunit=u,FILE=name(1:na)//'_'//c1//'.fdata',            &
                     form='unformatted', status= 'unknown')
            elseif (LN.ge.10) then
                write(c2,'(I2)') LN
                inquire(FILE=name(1:na)//'_'//c2//'.fdata',exist = exist)
                if (exist) cycle
                OPEN(newunit=u,FILE=name(1:na)//'_'//c2//'.fdata',            &
                     form='unformatted', status= 'unknown')
            endif
        endif
        !Write out the new topography data
        icc = 0
        do nn = 1,L(LN)%inne
            i = L(LN)%in(1,nn) ; j = L(LN)%in(2,nn)
            temp1 = ONE ; temp2 = ONE
            if (i.eq.L(LN)%is-1) then
                 temp1 = ZERO
            elseif (-L(LN)%H(1,nn).gt.L(LN)%ZK(nn).and.                       &
                -L(LN)%H(1,nn).gt.L(LN)%ZK(L(LN)%mn(i-1,j))) then
                temp1 = temp1 + (L(LN)%Z(L(LN)%KS) + L(LN)%H(1,nn)) / L(LN)%DZ
            endif
            if (L(LN)%DIM.eq.1.or.j.eq.L(LN)%js-1) then
                temp2 = ZERO
            elseif (L(LN)%DIM.eq.2.and.-L(LN)%H(2,nn).gt.L(LN)%ZK(nn).and.    &
                    -L(LN)%H(2,nn).gt.L(LN)%ZK(L(LN)%mn(i,j-1))) then
                temp2 = temp2 + (L(LN)%Z(L(LN)%KS) + L(LN)%H(2,nn)) / L(LN)%DZ
            endif
            temp3 = ONE + (L(LN)%Z(L(LN)%KS) - L(LN)%ZK(nn)) / L(LN)%DZ
            if (temp1.eq.ONE.and.temp2.eq.ONE.and.                            &
                L(LN)%MAN(nn).eq.ONE.and.temp3.eq.ONE)  cycle
            if (fg_xyz.eq.1) then
                !New format
                icc = icc+1
            else
                !Old format
                write(u) temp1 , temp2, L(LN)%MAN(nn), temp3    
            endif
        enddo
        if (fg_xyz.eq.1) then
            !New format
            write(u) icc,0
            do nn = 1,L(LN)%inne
                i = L(LN)%in(1,nn) ; j = L(LN)%in(2,nn)
                temp1 = ONE ; temp2 = ONE
                if (i.eq.L(LN)%is-1) then
                     temp1 = ZERO
                elseif (-L(LN)%H(1,nn).gt.L(LN)%ZK(nn).and.                   &
                        -L(LN)%H(1,nn).gt.L(LN)%ZK(L(LN)%mn(i-1,j))) then
                    temp1 = temp1 + (L(LN)%Z(L(LN)%KS) + L(LN)%H(1,nn))       &
                                   / L(LN)%DZ
                endif
                if (L(LN)%DIM.eq.1.or.j.eq.L(LN)%js-1) then
                    temp2 = ZERO
                elseif (L(LN)%DIM.eq.2.and.-L(LN)%H(2,nn).gt.L(LN)%ZK(nn).and.&
                        -L(LN)%H(2,nn).gt.L(LN)%ZK(L(LN)%mn(i,j-1))) then
                    temp2 = temp2 + (L(LN)%Z(L(LN)%KS) + L(LN)%H(2,nn))       &
                                   / L(LN)%DZ
                endif
                temp3 = ONE + (L(LN)%Z(L(LN)%KS) - L(LN)%ZK(nn)) / L(LN)%DZ
                if (temp1.eq.ONE.and.temp2.eq.ONE.and.                        &
                    L(LN)%MAN(nn).eq.ONE.and.temp3.eq.ONE)  cycle
                write(u) nn, temp1 , temp2, L(LN)%MAN(nn), temp3  
            enddo
        else
            !Old format
            do j = L(LN)%js-1,L(LN)%je
                do k = L(LN)%ks-1,L(LN)%ke
                    do i = L(LN)%is-1,L(LN)%ie
                        if (L(LN)%nf(i,k,j).eq.-1.and.                        &
                            L(LN)%nfb(i,k,j).ge.200) then
                            write(u) ZERO,ZERO,ZERO
                        endif
                    enddo
                enddo
            enddo
        endif
        CLOSE(u)
    ENDDO
!$omp end parallel do
    if (.not.exist) write(6,*) 'write end =',name(1:na)//'.fdata'
    !
endsubroutine outputF2D
!
subroutine set_title2D
    use variables, only: nkai, t, tw, ko_t, ko_tmp, ko_out, testa, title2D,   &
                         icong, n, t_end, ZERO
    implicit none
    integer       :: n1, n2, nslash, ns, ne
    character*255 :: num
    
    ko_t = ko_t + 1
    call owari(testa,n2) ; call hajime(testa,n1) ; call sura(testa,nslash)
    if (t.lt.ZERO.or.icong.gt.100) then
        title2D = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'2D.error'
        ko_t  = -1
    elseif (ko_t.eq.ko_tmp.or.t.ge.t_end.or.n.eq.nkai) then
        call ttonum(t,num,ns,ne)
        title2D = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'2D.t'//num(ns:ne)
        ko_t = 0
    elseif (ko_out.eq.1) then ! sfdata Çç◊Ç©Ç≠èoóÕÇ∑ÇÈèÍçá
        call ttonum(t,num,ns,ne)
        title2D = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'2D.t'//num(ns:ne)
    else
        title2D = testa(n1:n2)//'/VIS/'//testa(nslash+1:n2)//'2D.temp'
    endif
end subroutine set_title2D     