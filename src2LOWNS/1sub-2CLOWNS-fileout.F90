!%%%%%%%%%%%%%%%%%%%%%%% FILE: 1sub-2CLOWNS-fileout.F90 %%%%%%%%%%%%%%%%%%%%%%%
subroutine read_max_data_MAT
    use arrays, only: mn2d, FSLmax, U_3Dmax
    use variables, only: js, je, is, ie, input
    implicit none
    integer :: fnum, i, j, nx
    real*8  :: dummy, gomi
    logical :: exist
    call owari(INPUT,nx)
    do i = 6,9
        inquire(file = INPUT(1:nx-i)//'MAX.DAT',exist = exist)
        if (exist) exit
    enddo
    if (.not.exist) then
        write(6,*) 'Warning: Could not find 3D MAX data'
        return
    endif
    open(newunit = fnum,file = INPUT(1:nx-i)//'MAX.DAT',status='old')
    read(fnum,*) ! Skip reading the 
    read(fnum,*) ! x and y vectors
    do while (.true.) ! Read the remaining information
        read(fnum,*,END=123) i, j, FSLmax(mn2d(i,j)), gomi,                   &
                             U_3Dmax(0,mn2d(i,j)), dummy
        if (ubound(U_3Dmax,1).eq.1) U_3Dmax(1,mn2d(i,j)) = dummy 
    enddo
123 close(fnum)
    write(6,*) '[HOST]read end maxdata (DAT)',INPUT(1:nx-i)//'MAX.DAT'
endsubroutine read_max_data_MAT
!    
subroutine write_max_data_MAT
    use arrays, only: mn2d, FSLmax, U_3Dmax, IWD, x, in2d
    use variables, only: js, je, is, ie, title, inn2d, ZERO
    implicit none
    integer :: fnum, nx, nn, i, j
    real*8  :: dummy, gomi
    !
    call owari(title,nx)
    open(newunit = fnum,file = title(1:nx)//'MAX.DAT', status = 'unknown')
        write(fnum,'(9999(F13.3))') (X(0,i),i = is,ie)
        write(fnum,'(9999(F13.3))') (X(1,j),j = js,je)
        do nn = 1, inn2d
            if (IWD(nn)+FSLmax(nn).le.ZERO) cycle
            i = in2d(0,nn) ; j = in2d(1,nn)      
            if (ubound(U_3Dmax,1).eq.0) dummy = ZERO
            if (ubound(U_3Dmax,1).eq.1) dummy = U_3Dmax(1,nn)
            write(fnum,77) i, j, FSLmax(nn), IWD(nn) + FSLmax(nn),            &
                           U_3Dmax(0,nn), dummy
        enddo
    close(fnum)
77  FORMAT (2I5,4E15.7)
    write(6,*) '[HOST]read end maxdata (DAT)',title(1:nx)//'MAX.DAT'
endsubroutine write_max_data_MAT  
! Writes out local outputs for wave phenomenon    
subroutine local_data_output(iflag)
    use variables, only: testa, output_file, locnum
    use arrays, only: locout, filenums
    use ifport, only: makedirqq
    implicit none
    integer, intent(in) :: iflag
    integer :: unit, n1
    logical :: exist, result
    character*255 :: outfolder
    !
    !=========================================================================!
    !       This subroutine outputs time series of data                       ! 
    !  (water levels, velocities, momentum fluxes, turbulent quantities)      ! 
    !  at whatever locations are specified in the data-file : output_file     !
    !=========================================================================!
    !
    IFLAG_IF: if (iflag.eq.0) then
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !----------- Reading output data file and initialising ---------------!
        !------------------------ iflag = 0 ----------------------------------!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        ! 
        n1 = len_trim(testa)
        !Create the output directory
        inquire(file = testa(1:n1)//'/', exist = exist)
        if (.not.exist) result = makedirqq(testa(1:n1))
        !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !  Creat the VIS subdirectory and                                     !
        !  Write out S2D,M2D and window.init files for visualization          !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        outfolder = testa(1:n1)//'/VIS/'
        n1  =len_trim(outfolder)
        inquire(file = outfolder(1:n1),exist = exist)
        if (.not.exist) result = makedirqq(outfolder(1:n1))
        !Call subroutine that creates S2D,M2D window.init files
        if (.not.exist) call visfiles(outfolder)
        !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !     Making the DAT subdirectory                                     !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        outfolder = outfolder(1:n1-5)//'/DAT/'
        n1 = len_trim(outfolder)
        inquire(file = outfolder(1:n1),exist = exist)
        if (.not.exist) result = makedirqq(outfolder(1:n1))
        !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !     Read the output_file data if required                           !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        ! Return if no outputs to file are required
        if (output_file(1:3).eq.'NON') then
            write(6,*) 'No output_file; skipping make DAT files'
            return
        endif
        ! Otherwise...
        ! Define locations and what to output 
        ! at the various locations
        open(9,file = output_file,status='old')
        !
            read(9,*) locnum
            allocate(locout(1:locnum),filenums(1:locnum))
            read(9,*) !Don't read header line
            read(9,*) locout
        !
        close(9)
        !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !     Convert x and y coordinates to and i and j etc                  !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !
        call FINDXY
        !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !     Open files in the DAT subdirectory                              !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !
        call OPEN_DAT(outfolder)
!
    elseif (iflag.eq.1) then
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !----------- Writing current information to the output files ---------!
        !------------------------ iflag = 1 ----------------------------------!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !Except where no outputs are required           
        if (output_file(1:3).eq.'NON') return
        call WRITE_OUT_DAT
    endif IFLAG_IF
!
end subroutine local_data_output  
!
subroutine visfiles(outfolder)
    use variables, only: LNUM, fg_xyz, input, inputxyz, inputf0, inputxyz2D,  &
                         inputf2D, is, ie, js, je, ks, ZERO, IM, testa, inflow
    use arrays, only: L, x
    use ifport, only: makedirqq
    implicit none
    character*255,intent(in) :: outfolder
    real*8 :: han_x, han_y
    integer :: n1, na1, na2, na3, ns, ichop, ixyz, i0, k0, j0, LN, unit
    character*255 :: movie, filedummy
    character*10 :: jpegdir
    character*1  :: c1
    character*2  :: c2
    logical :: exist, result
    !-------------------------------------------------------------------------!
    ! This subroutine automatically outputs the S2D, M2D and                  !
    !     window.init files into the VIS outfolder for visualization          !
    !-------------------------------------------------------------------------!
    !
    n1  = len_trim(outfolder)
    call sura(testa,ns) ; na3 = len_trim(testa)
    jpegdir = './JPEG'
    ichop = 0 !Indicates opening of the window at water level
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !              Print out the 2DH visualization handles                    ! 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    if (IM.eq.1.or.IM.eq.2) then
        na1 = len_trim(inputxyz2D) ; na2 = len_trim(inputf2D)
        !S2D_2D
        filedummy = outfolder(1:n1)//'S2D_2D'
        inquire(file = filedummy, exist = exist)
        if (.not.exist) then
            open(newunit=unit,file = filedummy,status='new',action='write')
#include "S2D_2D.inc"        
            close(unit)
        endif
        !M2D_2D
        filedummy=outfolder(1:n1)//'M2D_2D'
        inquire(file = filedummy, exist = exist)
        if (.not.exist) then
            open(newunit=unit,file=filedummy,status='new',action='write')
#include "M2D_2D.inc"        
            close(unit)
        endif 
        !Window_2D.init
        do LN = 1,LNUM
            if (LN.lt.10) then
                write(c1,'(I1)') LN
                filedummy = outfolder(1:n1)//'window_'//c1//'.init'
            elseif (LN.ge.10) then
                write(c2,'(I2)') LN
                filedummy = outfolder(1:n1)//'window_'//c2//'.init'
            endif
            inquire(file = filedummy, exist = exist)
            if (.not.exist) then    
                !Preparing variables for window.init outputs
                ixyz = L(LN)%DIM - 1 !X-Z for 1D, X-Y plane for 2D
                i0 = 0.5 * (L(LN)%is+L(LN)%ie) 
                j0 = 0.5 * (L(LN)%js+L(LN)%je) ; k0 = L(LN)%ks !i,j,k number 
                han_x = 0.52d0 * (L(LN)%X(L(LN)%is)+L(LN)%X(L(LN)%ie))
                han_y = ZERO
                open(newunit=unit,file=filedummy,status='new',action='write')
#include "windowinit.inc"
                close(unit)  
            endif
        enddo
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !              Print out the 3D visualization handles                     ! 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    if (IM.eq.1.or.IM.eq.3) then
        na1 = len_trim(inputxyz) ; na2 = len_trim(inputf0)
        na3 = len_trim(input)
        !S2D(obj)
        filedummy = outfolder(1:n1)//'S2Dobj'
        inquire(file = filedummy, exist = exist)
        if (.not.exist) then
            open(newunit=unit,file = filedummy,status='new',action='write')
#include "S2Dobj.inc"
            close(unit)  
        endif
        !S3D
        filedummy = outfolder(1:n1)//'S3D'
        inquire(file = filedummy,exist = exist)
        if (.not.exist) then
            open(newunit=unit,file=filedummy,status='new',action='write')
#include "S3D.inc"
            close(unit)
        endif
        !M2D(obj)
        filedummy = outfolder(1:n1)//'M2Dobj'
        inquire(file = filedummy,exist = exist)
        if (.not.exist) then
            open(newunit=unit,file=filedummy,status='new',action='write')
#include "M2Dobj.inc"
            close(unit)
        endif
        !M3D
        filedummy = outfolder(1:n1)//'M3D'
        inquire(file = filedummy,exist = exist)
        if (.not.exist) then
            open(newunit=unit,file=filedummy,status='new',action='write')
#include "M3D.inc"
            close(unit)
        endif
        !
        !Window.init
        filedummy = outfolder(1:n1)//'window.init'
        inquire(file = filedummy, exist = exist)
        if (.not.exist) then    
            !Preparing variables for window.init outputs
            ixyz = 0  !X-Z plane
            i0 = 0.5 * (is+ie) ; j0 = 0.5 * (js+je) ; k0 = ks !i,j,k number
            han_x = 0.52d0 * (x(0,is)+x(0,ie)) ; han_y = ZERO
            open(newunit=unit,file=filedummy,status='new',action='write')
#include "windowinit.inc"
            close(unit)  
        endif
    endif
    !Create movie directory
    movie = outfolder(1:n1)//'/JPEG'
    na1 = len_trim(movie)
    inquire(file = movie,exist = exist)
    if (.not.exist) result = makedirqq(movie(1:na1))
!   
endsubroutine visfiles
!
subroutine FINDXY
    use variables, only: locnum, is, js, ie, je, ks, ke, LNUM
    use arrays, only: locout, x, nff, mn, L
    implicit none
    integer :: nn, i, j, k, LN, iss, iee, jss, jee
    real*8  :: xx_m, yy_m, xx_p, yy_p, x_m, y_m
    logical :: found
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !             Converts x and y coordinates to and i and j,                !
    !                   and finds the correct layer etc                       !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    FINDXY_LOOP: do nn = 1,locnum
        if (locout(nn)%locname(1:5).eq.'RUNUP') cycle
        if (locout(nn)%locname(1:4).eq.'CTAN') cycle
        if (locout(nn)%locname(1:6).eq.'ENERGY') cycle
        ! For AREA input
        if (locout(nn)%locname(1:4).eq.'AREA') then
            LN = locout(nn)%layer 
            if (LN.eq.0) stop 'Need to define layer number for AREA' 
            if (locout(nn)%x_int.ne.0) then
                if (locout(nn)%ks.eq.0) stop 'Define ks for AREA'
                if (L(LN)%DIM.eq.2) then
                    if (locout(nn)%y_int.eq.0) stop 'Define y_int for AREA'
                    if (locout(nn)%ke.eq.0) stop 'Define ke for AREA'
                endif
            endif
            cycle
        endif
        ! For XSEC input
        if (locout(nn)%locname(1:4).eq.'XSEC') then
            LN = locout(nn)%layer 
            if (LN.eq.0) stop 'Need to define layer number for XSEC'
            if (locout(nn)%ks.eq.0) then
                if (LN.le.LNUM) then
                     !2D Layer
                    if (locout(nn)%x_int.eq.0) locout(nn)%ks = L(LN)%is
                    if (locout(nn)%y_int.eq.0) locout(nn)%ks = L(LN)%js
                else !3D Layer
                    if (locout(nn)%x_int.eq.0) locout(nn)%ks = is
                    if (locout(nn)%y_int.eq.0) locout(nn)%ks = js                
                endif
            endif
            if (locout(nn)%ke.eq.0) then
                if (LN.le.LNUM) then
                     !2D Layer
                    if (locout(nn)%x_int.eq.0) locout(nn)%ke = L(LN)%ie-1
                    if (locout(nn)%y_int.eq.0) locout(nn)%ke = L(LN)%je-1
                else !3D Layer
                    if (locout(nn)%x_int.eq.0) locout(nn)%ke = ie-1
                    if (locout(nn)%y_int.eq.0) locout(nn)%ke = je-1                
                endif
            endif 
            cycle
        endif
        !If specified as x and y coordinates, find the
        !i and j integers of the output locations for faster printing
        if (locout(nn)%x_int.eq.0.and.locout(nn)%y_int.eq.0) then
            !Use the specified coordinates in m
            x_m = locout(nn)%x_m ; y_m = locout(nn)%y_m
            found = .false.
            LAY_DO: do LN = LNUM + 1, 1, -1 
                !Loop over the layers from inside to out so we get
                !the information from the most detailed layer
                if (LN.gt.LNUM) then
                    !3D Layer
                    iss = is; iee = ie
                    jss = js; jee = je 
                else
                    !2D Layer
                    iss = L(LN)%is; iee = L(LN)%ie
                    jss = L(LN)%js; jee = L(LN)%je
                endif
                do i = iss + 1, iee
                    if (LN.gt.LNUM) then
                        xx_p = X(0,i)     ; xx_m = X(0,i-1)
                    elseif (LN.le.LNUM) then
                        xx_p = L(LN)%X(i) ; xx_m = L(LN)%X(i-1)
                    endif
                    if (xx_p.ge.x_m.and.xx_m.lt.x_m) then
                        locout(nn)%x_int = i-1
                        do j = jss + 1, jee
                            if (LN.gt.LNUM) then
                                yy_p = X(1,j)     ; yy_m = X(1,j-1)
                            elseif (LN.le.LNUM) then
                                yy_p = L(LN)%Y(j) ; yy_m = L(LN)%Y(j-1)
                            endif
                            if (yy_p.ge.y_m.and.yy_m.lt.y_m) then
                                locout(nn)%y_int = j-1
                                found = .true.; exit
                            endif
                        enddo
                        exit
                    endif
                enddo
                if (found) then
                    locout(nn)%layer = LN; exit
                endif
            enddo LAY_DO  
            if (.not.found) then
                write(6,*) nn
                stop 'Error: the desired location is outside the bounds of the calculation domain'
            endif
        elseif (locout(nn)%x_int.eq.0.or.locout(nn)%y_int.eq.0) then
            stop 'Only one of the integer values are defined for a location output'
        endif
        if (locout(nn)%layer.le.LNUM) then
            if (locout(nn)%locname(1:2).eq.'VP'.or.                           &
                locout(nn)%locname(1:2).eq.'PS') stop 'PS or VP in 2D domain'
            cycle !(In 2D domain)
        endif
        !Get the vertical bounds of 
        !Lower bound
        i = locout(nn)%x_int; j = locout(nn)%y_int
        if (locout(nn)%ks.eq.0) then
            do k = ks,ke-1
                if (nff(mn(i,k,j))%f.eq.-1) cycle
                locout(nn)%ks = k
                exit
            enddo
        endif
        !Upper bound
        if (locout(nn)%ke.eq.0) then
            do k = ke-1,locout(nn)%ks,-1
                if (nff(mn(i,k,j))%f.eq.-1) cycle
                locout(nn)%ke = k
                exit
            enddo
        endif   
    enddo FINDXY_LOOP
!        
end subroutine FINDXY
        !
subroutine OPEN_DAT(outfolder)
    use variables, only: locnum, testa, t, HALF, ONE, LNUM
    use arrays, only: locout, dim3d, filenums, x, L
    use type_list, only: header
    implicit none
    character*255,intent(in) :: outfolder
    character*255 :: filedummy
    real*8,allocatable :: dummy(:,:)
    real*8 :: angle
    integer :: n1, n2, n3, nn, i, j, k, nslash, rows, jmax, filenum, LN, unit
    logical :: exist
    type(header) :: txt
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !          Open files in the DAT subdirectory to write the                !  
    !           data to. May be an existing file where we keep                !
    !              the information up to the current time                     !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Initialise header text
    txt%t     = '    t (s)   '
    txt%x     = '    x (m)   '
    txt%eta   = '  eta (m)   '
    txt%depth = '  depth (m) '
    txt%u     = '  u (m/s)   '
    txt%v     = '  v (m/s)   '
    txt%w     = '  w (m/s)   '
    txt%p     = '  p (kPa)   '
    txt%mom   = '  M (m2/s)  '
    txt%tke   = ' TKE (m2/s2)'
    txt%e     = '  e (m2/s3) '
    txt%nu_t  = ' nu_t (m2/s)' 
    txt%tang  = ' FS slope ()'
    txt%area  = '  Area (m2) '
    txt%Mass= ' Mass (m3)  '
    txt%Ep  = ' Ep (m4)    '
    txt%Ek  = ' Ek (m5s-2) '
    !
    n1  = len_trim(outfolder)
    n2 = len_trim(testa) ; call sura(testa,nslash)
    OUTPUT_LOOP: do nn = 1,locnum
        n3 = len_trim(locout(nn)%locname)
        LN = locout(nn)%layer
        !Set the file names
        filedummy = outfolder(1:n1)//locout(nn)%locname(1:n3)//'.dat'
        !Check if the file already exists or not
        inquire(file = filedummy, exist = exist)
        !Open the file
        IF_EXIST: if (exist) then
            !IF the file already exists than we need to 
            !Read the current data into a dummy array
            rows = -2 ; jmax = 5 !Default
            if (locout(nn)%locname(1:5).eq.'RUNUP') jmax = 6
            if (locout(nn)%locname(1:4).eq.'TURB') jmax = 6
            if (.not.dim3d) jmax = jmax - 1
            if (locout(nn)%locname(1:4).eq.'AREA') jmax = 2
            if (locout(nn)%locname(1:4).eq.'CTAN') jmax = 3
            if (locout(nn)%locname(1:6).eq.'ENERGY') jmax = 4
            if (locout(nn)%locname(1:4).eq.'XSEC') then
                if (locout(nn)%x_int.eq.0.or.locout(nn)%y_int.eq.0) then
                    jmax = locout(nn)%ke + 2 - locout(nn)%ks
                else
                    jmax = 1 + (max(locout(nn)%ke - locout(nn)%y_int,         &
                                    locout(nn)%ks - locout(nn)%x_int) + 1 ) * 2
                endif
            endif
            if (locout(nn)%locname(1:2).eq.'VP') then
                jmax =  1 + 3 * (locout(nn)%ke - locout(nn)%ks + 1)
            elseif (locout(nn)%locname(1:2).eq.'PS') then
                jmax =  2 + locout(nn)%ke - locout(nn)%ks
            endif
            !Open the old file
            open(newunit = unit,file = filedummy,status = 'old',action ='read')
                do while (.true.)
                    read(unit,*,END=123) 
                    rows = rows+1
                enddo           
                !
123             allocate(dummy(rows,jmax))
                rewind(unit)
                read(unit,*) !Skip header
                read(unit,*) !Skip header
                do i = 1,rows
                    read(unit,*) (dummy(i,j),j=1,jmax)
                enddo
            close(unit)
            !Open the existing file
            open(newunit = filenum,file = filedummy,                          &
                 status = 'old', action = 'write')
        else
            !IF the file does not exist then lets open a new one
            open(newunit = filenum,file = filedummy,                          &
                 status = 'new', action = 'write')
        endif IF_EXIST
        filenums(nn) = filenum
        !
        !Write banner to file 
        if (locout(nn)%locname(1:5).eq.'RUNUP') then
            ! Finds the location of the shoreline and give its details
            if (locout(nn)%x_int.eq.0) & !For runup along the x-axis
                write(filenum,19) 'Calc. Case: ',testa(nslash+1:n2),          &
                                  '  Loc.: Layer =',LN,                       &
                                  'I= Shoreline  J=',locout(nn)%y_int 
            if (locout(nn)%y_int.eq.0) & !For runup along the y-axis
                write(filenum,19) 'Calc. Case: ',testa(nslash+1:n2),          &
                                  '  Loc.: Layer =',LN,                       & 
                                  'J= Shoreline  I=',locout(nn)%x_int
            write(filenum,*) txt%t,txt%x,txt%eta,txt%depth,txt%u,txt%v 
        elseif (locout(nn)%locname(1:4).eq.'XSEC') then
            ! Cross-section of free-surface
            if (locout(nn)%x_int.eq.0) then
                !Cross-section along the x-axis
                write(filenum,19) 'Calc. Case: ',testa(nslash+1:n2),          &
                                  '  Loc.: Layer =',LN,                       &
                                  'I= All cells  J=',locout(nn)%y_int 
                if (locout(nn)%locname(5:5).eq.'M') then
                    write(filenum,*) txt%t,txt%mom,                           &
                                     ' X values are printed first:'
                else
                    write(filenum,*) txt%t,txt%eta,                           &
                                     ' X values are printed first:'
                endif
                if (.not.exist) then
                    if (LN.le.LNUM) then
                        write(filenum,'(9999(F12.4))') t,                     &
                                      ( HALF * ( L(LN)%X(i) + L(LN)%X(i+1) ), &
                                        i = locout(nn)%ks,locout(nn)%ke )
                    else
                        write(filenum,'(9999(F12.4))') t,                     &
                                      ( HALF * ( X(0,i) + X(0,i+1) ),         &
                                        i = locout(nn)%ks,locout(nn)%ke )   
                    endif
                endif  
            elseif (locout(nn)%y_int.eq.0) then
                !Cross-section along the y-axis
                write(filenum,19) 'Calc. Case: ',testa(nslash+1:n2),          &
                                  '  Loc.: Layer =',LN,                       &
                                  'J= All cells  I=',locout(nn)%x_int 
                if (locout(nn)%locname(5:5).eq.'M') then
                    write(filenum,*) txt%t,txt%mom,                           &
                                     ' Y values are printed first:'
                else
                    write(filenum,*) txt%t,txt%eta,                           &
                                     ' Y values are printed first:'
                endif
                if (.not.exist) then
                    if (LN.le.LNUM) then
                        write(filenum,'(9999(F12.4))') t,                     &
                                      ( HALF * ( L(LN)%Y(j) + L(LN)%Y(j+1) ), &
                                        j = locout(nn)%ks,locout(nn)%ke )
                    else
                        write(filenum,'(9999(F12.4))') t,                     &
                                      ( HALF * ( X(1,j) + X(1,j+1) ),         &
                                        j = locout(nn)%ks,locout(nn)%ke )    
                    endif
                endif
            else
                !Cross-section along angle
                write(filenum,18) 'Calc. Case: ',testa(nslash+1:n2),          &
                                  '  Loc.: Layer =',LN,                       &
                                  'I= ',locout(nn)%x_int,'J= ',locout(nn)%y_int
                if (locout(nn)%locname(5:5).eq.'M') then
                    write(filenum,*) txt%t,txt%mom,                           &
                                     ' X,Y values are printed first:'
                else
                    write(filenum,*) txt%t,txt%eta,                           &
                                     ' X,Y values are printed first:'
                endif
                angle = real(locout(nn)%ke - locout(nn)%y_int)                &
                      / real(locout(nn)%ks - locout(nn)%x_int)
                if (.not.exist) then
                    if (LN.le.LNUM) then
                        if (abs(angle).ge.ONE) then
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( ( HALF * ( L(LN)%X(i) + L(LN)%X(i+1) ),         &
                            i = locout(nn)%x_int +                            &
                            nint(real(j-locout(nn)%y_int) / angle),           &
                            locout(nn)%x_int +                                &
                            nint(real(j-locout(nn)%y_int) / angle) ),         &
                            j = locout(nn)%y_int,locout(nn)%ke),              &
                            ( HALF * ( L(LN)%Y(j) + L(LN)%Y(j+1) ),           &           
                            j = locout(nn)%y_int, locout(nn)%ke )
                        elseif (abs(angle).lt.ONE) then 
                            write(filenum,'(9999(F12.4))') t,                 &
                            (HALF * ( L(LN)%X(i) + L(LN)%X(i+1) ),            &
                            i = locout(nn)%x_int, locout(nn)%ks),             &
                            ( (HALF * ( L(LN)%Y(j) + L(LN)%Y(j+1) ),          &
                            j = locout(nn)%y_int +                            &
                            nint(real(i-locout(nn)%x_int) * angle),           &
                            locout(nn)%y_int + nint(real(i-locout(nn)%x_int) *&
                            angle)),i = locout(nn)%x_int,locout(nn)%ks )   
                        endif
                    else
                        if (abs(angle).ge.ONE) then
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( (HALF * (X(0,i) + X(0,i+1) ),                   &
                            i = locout(nn)%x_int +                            &
                            nint(real(j-locout(nn)%y_int) / angle),           &
                            locout(nn)%x_int +                                &
                            nint(real(j-locout(nn)%y_int) / angle) ),         &
                            j = locout(nn)%y_int, locout(nn)%ke),             &
                            ( HALF * ( X(1,j) + X(1,j+1) ),                   &
                            j = locout(nn)%y_int ,locout(nn)%ke )
                        elseif (abs(angle).lt.ONE) then 
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( HALF * ( X(0,i) + X(0,i+1) ),                   &
                            i = locout(nn)%x_int,locout(nn)%ks ),             &
                            ( (HALF*(X(1,j) + X(1,j+1) ),j = locout(nn)%y_int+&
                            nint(real(i-locout(nn)%x_int) * angle),           & 
                            locout(nn)%y_int +                                &
                            nint(real(i-locout(nn)%x_int) * angle) ),         &
                            i = locout(nn)%x_int,locout(nn)%ks )   
                        endif
                    endif
                endif            
            endif 
        elseif (locout(nn)%locname(1:4).eq.'CTAN'.or.                         &
                locout(nn)%locname(1:6).eq.'ENERGY') then
            write(filenum,20) 'Calc. Case: ',testa(nslash+1:n2),              &
                              ' Loc.: Layer =', LN,                           &
                              'J= All cells  I= All Cells '
            if (locout(nn)%locname(1:4).eq.'CTAN') then
                write(filenum,*) txt%t,txt%x,txt%tang
            elseif (locout(nn)%locname(1:6).eq.'ENERGY') then
                write(filenum,*) txt%t,txt%Mass,txt%Ep,txt%Ek
            endif
        elseif (locout(nn)%locname(1:4).eq.'AREA') then
            write(filenum,20) ' Calc. Case: ',testa(nslash+1:n2),             &
                              ' Loc.: Layer =', LN,                           &
                              'I,J= All cells: Inund. Area'
            write(filenum,*) txt%t,txt%area
        elseif (locout(nn)%locname(1:2).eq.'VP'.or.                           &
                locout(nn)%locname(1:2).eq.'PS') then
            ! May be used to get velocity profile in vertical direction
            ! over a cell column or even at just a single cell
            write(filenum,18) 'Calc. Case: ',testa(nslash+1:n2),              &
                              ' Loc.: Layer =', LN,                           &
                              'I= ',locout(nn)%x_int,'J= ',locout(nn)%y_int
            if (locout(nn)%locname(1:2).eq.'VP') then
                write(filenum,*) txt%t,txt%u,txt%v,txt%w,                     & 
                                 ' Z values are printed first:'
            elseif (locout(nn)%locname(1:2).eq.'PS') then
                write(filenum,*) txt%t,txt%p,                                 &
                                 ' Z values are printed first:'
            endif
            !
            if (.not.exist) then
                write(filenum,'(999(F12.4))') t,(( HALF * (X(2,k) + X(2,k+1)),&
                                                   k = locout(nn)%ks,         &
                                                   locout(nn)%ke ), i = 1,3 )  
            endif 
        else ! Will give various data at a specified location
                !(can be thought of as a wave gauge)
            write(filenum,18) 'Calc. Case: ',testa(nslash+1:n2),              &
                              ' Loc.: Layer =', LN,                           &
                              'I= ',locout(nn)%x_int,'J= ',locout(nn)%y_int
            if (locout(nn)%locname(1:4).eq.'TURB') then
                write(filenum,*) txt%t,txt%u,txt%w,txt%tke,txt%e,txt%nu_t
            else
                write(filenum,*) txt%t,txt%eta,txt%depth,txt%u,txt%v
            endif
        endif
        !
        if (exist) then 
        !IF we have an existing file write the 
        !existing data up until the current time
            do i = 1,rows
                if (dummy(i,1).gt.t) exit
                if (locout(nn)%locname(1:4).eq.'XSEC') then
                    write(filenum,'(9999(F12.4))') (dummy(i,j),j=1,jmax)
                elseif (locout(nn)%locname(1:5).eq.'RUNUP') then
                    write(filenum,17) (dummy(i,j),j=1,jmax)
                elseif (locout(nn)%locname(1:2).eq.'VP'.or.                   &
                        locout(nn)%locname(1:2).eq.'PS') then
                    write(filenum,'(999(F12.4))') (dummy(i,j),j = 1,jmax )
                elseif (locout(nn)%locname(1:4).eq.'AREA') then
                    write(filenum,21) (dummy(i,j),j=1,jmax)
                elseif (locout(nn)%locname(1:4).eq.'CTAN') then
                    write(filenum,'(3(F12.4))') (dummy(i,j),j=1,jmax)
                elseif (locout(nn)%locname(1:6).eq.'ENERGY') then
                    write(filenum,'(F12.4,3(E14.5))') (dummy(i,j),j=1,jmax)
                else
                    write(filenum,16) (dummy(i,j),j=1,jmax)
                endif
            enddo
            deallocate(dummy)
        endif
    enddo OUTPUT_LOOP  
    16 FORMAT (5f12.4)
    17 FORMAT (f12.4,f12.3,4f12.4)
    18 FORMAT (A12,A,A15,I3,2x,A3,I5,2x,A3,I5)  
    19 FORMAT (A12,A,A15,I3,2x,A16,I5)  
    20 FORMAT (A12,A,A15,I3,2x,A27)  
    21 FORMAT (f12.4,f14.1)
    22 FORMAT (f12.4,5e12.3) 
end subroutine OPEN_DAT
!
subroutine WRITE_OUT_DAT
    use variables, only: locnum, is, ie, je, js, ke, ks, t, isfn2, LNUM,      &
                         ONE, HALF, ZERO, LITTLE, inns, inne
    use arrays, only: filenums, locout, x, f, nff, mn, dx, Flux_3D, mn2d, FSL,&
                      u, p, r, d, in, nor, IWD, sfn, U_3D, L, dim3d, kn, jn
    implicit none
    integer :: nn2d, nn, i, j, k, ii, jj, iss, iee, filenum, LN, nnn, nip, njp
    real*8  :: xx_p, yy_p, xx_m, yy_m, x_m, y_m
    real*8  :: angle, Uvel, Vvel, area, ETA, Depth
    real*8  :: Mass, massnow, Ep, Ek, Uav, Vav, Wav
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !      Writes out the water levels, velocities and fluxes                 !  
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !
    WRITE_DAT_LOOP: do nn = 1, locnum
    !
        filenum = filenums(nn)
    !
        i  = locout(nn)%x_int; j = locout(nn)%y_int
        LN = locout(nn)%layer
        if (locout(nn)%locname(1:5).eq.'RUNUP') then
            if (i.eq.0) then
                ! Runup-up along x-axis
                jj = 0
                if (LN.le.LNUM) then
                    iee = L(LN)%ie; iss = L(LN)%is
                else
                    iee = ie; iss = is
                endif
            elseif (j.eq.0) then
                ! Runup-up along y-axis
                jj = 1
                if (LN.le.LNUM) then
                    iee = L(LN)%je; iss = L(LN)%js
                else
                    iee = je; iss = js
                endif
            endif
            do ii = iee-1,iss,-1
                if (jj.eq.0) i = ii
                if (jj.eq.1) j = ii
                if (LN.le.LNUM) then
                    if (jj.eq.0) xx_p = L(LN)%X(ii+1)
                    if (jj.eq.1) xx_p = L(LN)%Y(ii+1)
                    if (L(LN)%F(L(LN)%mn(i,j)).ne.ZERO) exit
                elseif (LN.gt.LNUM) then
                    xx_p = x(jj,ii+1)
                    do k = ke-1,ks,-1
                        if (nff(mn(i,k,j))%f.eq.-1) exit 
                            !Setting min depth to 1mm
                        if (f(mn(i,k,j))*dx(2,k).lt.LITTLE.or.                &
                            nff(mn(i,k,j))%f.eq.0) cycle
                        goto 7
                    enddo
                endif
            enddo 
        elseif (locout(nn)%locname(1:4).eq.'XSEC') then
            if (i.ne.0.and.j.ne.0) then
                angle = real(locout(nn)%ke - locout(nn)%y_int)                &
                      / real(locout(nn)%ks - locout(nn)%x_int)
            endif
            if (locout(nn)%locname(5:5).eq.'M') then
                ! Cross-section of momentum fluxes            
                if (LN.le.LNUM) then
                        !2D domain                          
                    if (i.eq.0) then               ! Cross-section along x-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                             ( L(LN)%Q(L(LN)%DIM,L(LN)%mn(i,j)),              &
                               i = locout(nn)%ks , locout(nn)%ke )
                    elseif (j.eq.0) then           ! Cross-section along y-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                             ( L(LN)%Q(1,L(LN)%mn(i,j)),                      &
                               j = locout(nn)%ks , locout(nn)%ke )    
                    else                                    
                        if (abs(angle).ge.ONE) then ! Cross-section along angle
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( (L(LN)%Q(1,L(LN)%mn(i,j)),i = locout(nn)%x_int +&
                            nint(real(j-locout(nn)%y_int) / angle),           &
                            locout(nn)%x_int +                                &
                            nint(real(j-locout(nn)%y_int) / angle)),          &
                            j = locout(nn)%y_int , locout(nn)%ke ),           &
                            ( (L(LN)%Q(2,L(LN)%mn(i,j)),i = locout(nn)%x_int +&
                            nint(real(j-locout(nn)%y_int) / angle),           &
                            locout(nn)%x_int +                                &  
                            nint(real(j-locout(nn)%y_int) / angle) ),         &
                            j = locout(nn)%y_int , locout(nn)%ke )
                        elseif (abs(angle).lt.ONE) then 
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( (L(LN)%Q(1,L(LN)%mn(i,j)),j = locout(nn)%y_int +&
                            nint(real(i-locout(nn)%x_int) * angle),           &
                            locout(nn)%y_int +                                &
                            nint(real(i-locout(nn)%x_int) * angle)),          &
                            i = locout(nn)%x_int , locout(nn)%ks),            &
                            ( (L(LN)%Q(2,L(LN)%mn(i,j)),j = locout(nn)%y_int +&
                            nint(real(i-locout(nn)%x_int) * angle),           &
                            locout(nn)%y_int +                                &
                            nint(real(i-locout(nn)%x_int) * angle) ),         &
                            i = locout(nn)%x_int , locout(nn)%ks )  
                        endif
                    endif
                else
                    if (i.eq.0) then               ! Cross-section along x-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                                       ( Flux_3D(ubound(Flux_3D,1),mn2d(i,j)),&
                                         i = locout(nn)%ks , locout(nn)%ke )    
                    elseif (j.eq.0) then           ! Cross-section along y-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                                       ( Flux_3D(0,mn2d(i,j)),                &
                                         j = locout(nn)%ks , locout(nn)%ke )    
                    else      
                        if (abs(angle).ge.ONE) then ! Cross-section along angle
                            write(filenum,'(9999(F12.4))') t,                 &
                            ((Flux_3D(0,mn2d(i,j)),i = locout(nn)%x_int +     &
                            nint(real(j-locout(nn)%y_int)/angle),             &
                            locout(nn)%x_int +                                &
                            nint(real(j-locout(nn)%y_int) / angle)),          &
                            j = locout(nn)%y_int,locout(nn)%ke),              &
                            ( (Flux_3D(1,mn2d(i,j)),i = locout(nn)%x_int +    &
                            nint(real(j-locout(nn)%y_int)/angle),             & 
                            locout(nn)%x_int +                                &
                            nint(real(j-locout(nn)%y_int) / angle)),          &
                            j = locout(nn)%y_int , locout(nn)%ke )
                        elseif (abs(angle).lt.ONE) then 
                            write(filenum,'(9999(F12.4))') t,                 &
                            ((Flux_3D(0,mn2d(i,j)),j = locout(nn)%y_int +     &
                            nint(real(i-locout(nn)%x_int)*angle),             &
                            locout(nn)%y_int +                                &
                            nint(real(i-locout(nn)%x_int) * angle)),          &
                            i = locout(nn)%x_int , locout(nn)%ks),            &
                            ( (Flux_3D(1,mn2d(i,j)),j = locout(nn)%y_int +    &
                            nint(real(i-locout(nn)%x_int)*angle),             &
                            locout(nn)%y_int +                                & 
                            nint(real(i-locout(nn)%x_int) * angle)),          &
                            i = locout(nn)%x_int , locout(nn)%ks )  
                        endif
                    endif
                endif
            else
                ! Cross-section of free-surface 
                if (LN.le.LNUM) then  !2D domain
                    if (i.eq.0) then               ! Cross-section along x-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                             ( L(LN)%ETAn(L(LN)%mn(i,j)),                     &
                               i = locout(nn)%ks , locout(nn)%ke )
                    elseif (j.eq.0) then           ! Cross-section along y-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                             ( L(LN)%ETAn(L(LN)%mn(i,j)),                     &
                               j = locout(nn)%ks , locout(nn)%ke )  
                    else           
                        if (abs(angle).ge.ONE) then ! Cross-section along angle
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( (L(LN)%ETAn(L(LN)%mn(i,j)),i = locout(nn)%x_int+&
                            nint(real(j-locout(nn)%y_int) * angle),           &
                            locout(nn)%x_int +                                &
                            nint(real(j-locout(nn)%y_int) * angle) ),         &
                            j = locout(nn)%y_int , locout(nn)%ke )
                        elseif (abs(angle).lt.ONE) then 
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( (L(LN)%ETAn(L(LN)%mn(i,j)),j = locout(nn)%y_int+&
                            nint(real(i-locout(nn)%x_int) * angle),           &
                            locout(nn)%y_int +                                &   
                            nint(real(i-locout(nn)%x_int) * angle) ),         &
                            i = locout(nn)%x_int , locout(nn)%ks )                     
                        endif
                    endif
                else ! 3D domain
                    if (i.eq.0) then               ! Cross-section along x-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                             ( FSL(mn2d(i,j)),i = locout(nn)%ks, locout(nn)%ke)
                    elseif (j.eq.0) then          ! Cross-section along y-axis
                        write(filenum,'(9999(F12.4))') t,                     &
                             ( FSL(mn2d(i,j)),j = locout(nn)%ks, locout(nn)%ke)  
                    else                                  
                        if (abs(angle).ge.ONE) then ! Cross-section along angle
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( (FSL(mn2d(i,j)),i = locout(nn)%x_int +          &
                            nint(real(j-locout(nn)%y_int)  *angle),           &
                            locout(nn)%x_int +                                & 
                            nint(real(j-locout(nn)%y_int) * angle) ),         &
                            j = locout(nn)%y_int, locout(nn)%ke )
                        elseif (abs(angle).lt.ONE) then 
                            write(filenum,'(9999(F12.4))') t,                 &
                            ( (FSL(mn2d(i,j)),j = locout(nn)%y_int +          &
                            nint(real(i-locout(nn)%x_int) * angle),           &
                            locout(nn)%y_int +                                &
                            nint(real(i-locout(nn)%x_int) * angle) ),         &
                            i = locout(nn)%x_int, locout(nn)%ks )                 
                        endif
                    endif   
                endif
            endif
            cycle
        elseif (locout(nn)%locname(1:4).eq.'AREA') then
            ! Calculates the flooded area
            area = ZERO
            if (LN.le.LNUM) then
                do nnn = 1,L(LN)%inne
                    i = L(LN)%in(1,nnn) ; j = L(LN)%in(2,nnn)
                    if (L(LN)%coupling_dim.eq.2.and.                          &
                        i.ge.L(LN)%is2.and.i.lt.L(LN)%ie2.and.                &
                        j.ge.L(LN)%js2.and.j.lt.L(LN)%je2) then
                        ! Cycle inside of the next layer domain
                        if (L(LN)%mn2(0,i,j).eq.0) cycle   
                    endif
                    ! Cycle inside of specified region (if exists)
                    if (locout(nn)%x_int.ne.0) then
                        if (i.ge.locout(nn)%x_int.and.i.le.locout(nn)%ks.and. &
                           ( L(LN)%DIM.eq.1.or.(j.ge.locout(nn)%y_int.and.    &
                             j.le.locout(nn)%ke) ) ) cycle
                    endif
                    if (L(LN)%F(nnn).gt.ZERO) then
                        area = area + sum(L(LN)%DX(:,nnn)**2)
                    endif
                enddo
            endif
            write(filenum,18) t, area
            cycle
        elseif (locout(nn)%locname(1:2).eq.'VP') then
            !Velocity profile over the cell column
            write(filenum,'(999(F12.4))') t, ( u(:,mn(i,k,j)),                &
                                             k = locout(nn)%ks, locout(nn)%ke )  
            cycle
        elseif (locout(nn)%locname(1:2).eq.'PS') then
            !Pressure profile over the cell column
            write(filenum,'(999(F12.4))') t, ( p(mn(i,k,j)),                  &
                                             k = locout(nn)%ks, locout(nn)%ke )  
            cycle 
        endif
7       if (locout(nn)%locname(1:4).eq.'TURB') then
            ! Get turb characteristics at bottom cell close to bed
            do k = ks,ke-1
                if (nff(mn(i,k,j))%f.eq.-1) cycle     ! Cycle 5 mm from the bed
                if (nff(mn(i,k-1,j))%f.eq.-1.and.dx(2,k).lt.0.01d0) cycle 
                if (nff(mn(i+1,k,j))%f.lt.0) then
                    Uvel = u(0,mn(i,k,j))
                elseif (nff(mn(i+1,k,j))%f.ge.0) then
                    Uvel = HALF * ( u(0,mn(i,k,j)) +  u(0,mn(i+1,k,j)) )
                endif
                write(filenum,19) t, Uvel, u(2,mn(i,k+1,j)),                  &
                                  r(mn(i,k,j)), d(mn(i,k,j))
                exit
            enddo
            cycle
        elseif (locout(nn)%locname(1:4).eq.'CTAN') then
            ! Cross-section of free-surface tangent             
            i = in(sfn(maxloc(abs(nor(sfn(1:isfn2),0)),1)))
            write(filenum,'(4(F12.4))') t, HALF*(x(0,i)+x(0,i+1)),            &
                                        maxval(abs(nor(sfn(1:isfn2),0)))
            !,maxval(abs(nor(sfn(1:isfn2),1))),maxval(abs(nor(sfn(1:isfn2),2)))
            cycle
        elseif (locout(nn)%locname(1:6).eq.'ENERGY') then
            ! Calculate wave energies and mass
            Mass = ZERO ; Ep = ZERO ; Ek = ZERO
            do nnn = inns,inne
                if (nff(nnn)%f.eq.-1) cycle
                i = in(nnn) ; j = jn(nnn) ; k = kn(nnn)
                massnow = F(nnn) * dx(0,i) * dx(1,j) * dx(2,k)
                ! Mass
                Mass = Mass + massnow
                ! Potential energy
                Ep   = Ep   + massnow * HALF * ( x(2,k) + x(2,k+1) )
                ! Kinetic energy
                Uav  = HALF * ( u(0,nnn) + u(0,mn(i+1,k,j)) )
                Vav  = HALF * ( u(1,nnn) + u(1,mn(i,k,j+1)) )
                Wav  = HALF * ( u(2,nnn) + u(2,mn(i,k+1,j)) )
                Ek   = Ek   + massnow * HALF * ( Uav * Uav + Vav * Vav + Wav * Wav )
            enddo
            write(filenum,'(F12.4,3(E14.5))') t, Mass, Ep, Ek
            cycle
        endif
        ! OUTPUTS the free surface, depth and depth-averaged velocities
        if (LN.le.LNUM) then !2D domain
            nn2d  = L(LN)%mn(i,j) ; nip = L(LN)%mn(i+1,j)
            ETA   = L(LN)%ETAn(nn2d)
            Depth = L(LN)%ETAn(nn2d) - L(LN)%ZK(nn2d)
            Uvel  = HALF * ( L(LN)%U(1,nn2d) + L(LN)%U(1,nip) )
            if (L(LN)%DIM.eq.2) then
                njp   = L(LN)%mn(i,j+1)
                Vvel  = HALF * ( L(LN)%U(2,nn2d) + L(LN)%U(2,njp) )
            else
                Vvel  = ZERO
            endif
        else                 !3D domain
            nn2d  = mn2d(i,j) ; nip = mn2d(i+1,j)
            ETA   = FSL(nn2d)
            Depth = IWD(nn2d) + FSL(nn2d)
            Uvel  = HALF * ( U_3D(0,nn2D) + U_3D(0,nip) )
            if (dim3d) then
                njp   = mn2d(i,j+1)
                Vvel  = HALF * ( U_3D(1,nn2D) + U_3D(1,njp) )
            else
                Vvel  = ZERO
            endif
        endif
        if (locout(nn)%locname(1:5).eq.'RUNUP') then
            write(filenum,17) t, xx_p, ETA, Depth, Uvel, Vvel
        else
            write(filenum,16) t, ETA, Depth, Uvel, Vvel
        endif
    !
    enddo WRITE_DAT_LOOP
    16 FORMAT (5f12.4)
    17 FORMAT (f12.4,f12.3,4f12.4)
    18 FORMAT (f12.4,f14.1)
    19 FORMAT (f12.4,5e12.3)   
end subroutine WRITE_OUT_DAT