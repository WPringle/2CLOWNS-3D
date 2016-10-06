if (LNUM.eq.1) then
	write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1)//' ../../'//inputf2D(1:na2)//' none $2 $1.t?.??? $1.t??.?? $1.t???.? $1.t????.?'
else
	if (fg_xyz.eq.0) write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1-4)//'_$2.xyz ../../'//inputf2D(1:na2-6)//'_$2.fdata non window_$2.init $1.t?.???_$2 $1.t??.??_$2 $1.t???.?_$2 $1.t????.?_$2'
	if (fg_xyz.eq.1) then
	if (inflow%wavetype(1:1).eq.'G') then
		write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1-5)//'_$2.xyzn '//testa(ns+1:na3)//'2D.t0.000_$2.fdata none window_$2.init $1.t?.???_$2 $1.t??.??_$2 $1.t???.?_$2 $1.t????.?_$2'
	else
		write(unit,'(A$)') '../../../../Vis/execute/GEN64.exe ../../'//inputxyz2D(1:na1-5)//'_$2.xyzn ../../'//inputf2D(1:na2-6)//'_$2.fdata none window_$2.init $1.t?.???_$2 $1.t??.??_$2 $1.t???.?_$2 $1.t????.?_$2'
	endif
	endif
endif
