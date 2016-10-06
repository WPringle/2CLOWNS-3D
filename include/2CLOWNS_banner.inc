
  !Top Banner
  write(6,*)      '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  write(6,*)      '                             2CLOWNS - 3D (2016)                              '
  write(6,*)      '          ==================================================================  '
  write(6,*)      '                                                                              '
  write(6,*)      '           -INCLUDES A 3D RANS BASED MODEL USING THE VOF METHOD FOR           '
  write(6,*)      '            FREE SURFACE AND K-EPS MODEL FOR TURBULENCE MODELLING             '
  write(6,*)      '                                                 ORIGINAL BY YONEYAMA (1998)  '
  write(6,*)      '                                                                              ' 
  write(6,*)      '           -CREATION OF 2DH MODEL AND TWO-WAY COUPLING WITH 3D MODEL          '
  write(6,*)      '            FOR WAVE-LIKE BEHAVIOUR AND TSUNAMI SOURCES, PLUS                 '
  write(6,*)      '            VARIOUS IMPROVEMENTS AND CHANGES IN THEORY AND CODING LAYOUT      '
  write(6,*)      '                                                 BY WILLIAM PRINGLE (2016)    '
  write(6,*)      '                                                                              ' 
  write(6,*)      '           -FOR USE OF INPUT/OUTPUT FILES AND VISUALIZATION PROGRAMS          '
  write(6,*)      '            DEVELOPED BY NOZOMU YONEYAMA                                      '
  write(6,*)      '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  write(6,*)      '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  !
  !input Info
  !====================================================================================================================
  !Basic Parameters
  !----------------
  if (IM.eq.1) then
	write(6,*) 'IM = 1: Performing Hybrid 2DH-3D Calculation'
  elseif (IM.eq.2) then
	write(6,*) 'IM = 2: Performing 1DH or 2DH Calculation'
  elseif (IM.eq.3) then
	write(6,*) 'IM = 3: Performing 2DV or 3D Calculation'
  endif
  write(6,*) 'Max. No. of time steps=',nkai
  write(6,*) 'Info/Local Data/Temp output cycle=',MWR, DWR, ko_tmp
  write(6,*) 'output sfdata (yes =1)', ko_out
  write(6,*) 'trange=',ZERO,'->',t_end,'[sec]'
  write(6,*) 'Time of first output=',tw0,'[sec]'
  write(6,*) 'Data Output interval=',dtw,'[sec]'
  write(6,*) 'Gravity=',g,'[ms^-2]'
  if (IM.eq.1.or.IM.eq.3) then
  write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%% 3D Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	  write(6,*) 'Allowable limit for error in continuity eqn',D0,'[ms^-2]'
	  write(6,*) 'Allowable error in Poisson Pressure matrix solver=',eps0
	  write(6,*) 'Viscosity=',nu,'[m^2s^-1]'
	  write(6,*) 'init turbulent intensity =',TI_b,'nu_t ratio=',nut_b
	  if (dt.eq.ZERO) then
		write(6,*) 'Courant set limit=',cl_set,'Max. timestep limit=',dtmax,'[sec]'
	  else
		write(6,*) 'Init. Timestep=',dt,'[sec], Max. timestep limit=',dtmax,'[sec]'
	  endif
  endif
  write(6,*) '%%%%%%%%%%%%%%%%%%%%%% Input/Output files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  !-----------
       dummyc='                                                           }' 
	   if (IM.eq.1.or.IM.eq.3) then
		   call owari(input,na1)
		   write(6,*) 'input   = ',input(1:na1),dummyc(10+na1:60) 
		   call owari(inputf0,na1) 
		   write(6,*) 'inputF3D= ',inputf0(1:na1),dummyc(10+na1:60),'3D Input files'     
		   call owari(inputxyz,na1) 
		   write(6,*) 'inputXYZ= ',inputxyz(1:na1),dummyc(10+na1:60)   
		   call owari(DRIINIT,na1) 
		   write(6,*) 'inputDRI= ',DRIINIT(1:na1),dummyc(10+na1:60)  
	   endif
	   if (IM.eq.1.or.IM.eq.2) then
		   call owari(input2D,na1)
		   write(6,*) 'input2D = ',input2D(1:na1),dummyc(10+na1:60) 
		   call owari(inputf2D,na1) 
		   write(6,*) 'inputF2D= ',inputf2D(1:na1),dummyc(10+na1:60),'2D Input files'     
		   call owari(inputxyz2D,na1) 
		   write(6,*) 'inputX2D= ',inputxyz2D(1:na1),dummyc(10+na1:60)  
	   endif  
       call owari(testa,na1)
       write(6,*) 'output name= ',testa(1:na1),dummyc(13+na1:59),'Output folder'
	   call owari(output_file,na1) 
       write(6,*) 'Output loc. file= ',output_file(1:na1),dummyc(13+na1:59)
  !
  write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  !-------------------
  if (IM.eq.1.or.IM.eq.2) then
	write(6,*) '%%%% 2DH Outer Layer %%%%'
	write(6,*) 'North = ',BC(2)%North,'  South = ',BC(2)%South,'  East= ',BC(2)%East,'  West= ',BC(2)%West
  endif
  if (IM.eq.1.or.IM.eq.3) then
	write(6,*) '%%%% 3D  Layer %%%%%%%%%%'
	write(6,*) 'North = ',BC(3)%North,'  South = ',BC(3)%South,'  East= ',BC(3)%East,'  West= ',BC(3)%West
  endif
  !
  write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%% Inflow Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  !-------------------
  write(6,*) 'Inflow wave type=',inflow%wavetype
  write(6,*) 'Target wave height=',inflow%h_wave,'[m]'
  write(6,*) 'Mean water depth=',inflow%h_shoki,'[m]'
  write(6,*) 'Angle of input wave relative to positive x-axis=',inflow%angle,'[deg]'
  write(6,*) 'Wave period=',inflow%period,'[s]'
