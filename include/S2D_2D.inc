
if (LNUM.eq.1) then
	!When using the built fdata
	write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1)//' ../../'//inputf2D(1:na2)//' none $1 $2'
else
	if (fg_xyz.eq.0) write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1-4)//'_$2.xyz ../../'//inputf2D(1:na2-6)//'_$2.fdata none $1_$2 window_$2.init'
	if (fg_xyz.eq.1) then
	if (inflow%wavetype.eq.'GC') then
		write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1-5)//'_$2.xyzn '//testa(ns+1:na3)//'2D.t0.000_$2.fdata none $1_$2 window_$2.init'
	elseif (inflow%wavetype.eq.'GT') then
		write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1-5)//'_$2.xyzn $1_$2.fdata none $1_$2 window_$2.init'
	else
		write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1-5)//'_$2.xyzn ../../'//inputf2D(1:na2-6)//'_$2.fdata none $1_$2 window_$2.init'
	endif
	endif
endif