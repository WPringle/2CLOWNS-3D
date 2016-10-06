       subroutine ttonum2(t,num,ns,ne)
       implicit none
       doubleprecision,intent(in)::t
       integer::nt
       integer,intent(out)::ns,ne
       character*255:: num
       character*255:: num1

       ns=1
       ne=5
       nt=nint(t*10000.)

       if(nt.lt.10) then
!CC -1 keta
       write(num1,'(I1)') nt
       ne=6
       num='0.000'//num1(1:1)
!       
       else if(nt.lt.100) then
!CC -1 keta
       write(num1,'(I2)') nt
       ne=6
       num='0.00'//num1(1:2)
!
       else if(nt.lt.1000) then
!CC -1 keta
       ne=6
       write(num1,'(I3)') nt
       num='0.0'//num1(1:3)
       
       else if(nt.lt.1.0d4) then
!CC 1 keta
       write(num1,'(I4)') nt
       num='0.'//num1(1:4)
       ne=6       
       else if(nt.lt.1.0d5) then
!CC 2 keta
       write(num,'(f6.4)') t
       ns=1
       ne=6
       else if(nt.lt.1.0d6) then
!CC 2 keta
       write(num,'(f7.4)') t
       ne=7       
       else if(t.lt.999.9) then
!CC 3 keta
       write(num,101) t

!c 101   format(f6.2)
!c       ne=6
 101   format(f5.1)
       else if(t.lt.9999.9) then
!CC 3 keta
       write(num,104) t
 104   format(f5.0)

       else if(t.lt.99999.9) then
!CC 3 keta
       write(num,105) t
 105   format(f6.0)
       ne=6
     
       else if(t.lt.999999.9) then
!CC 3 keta
       write(num,106) t
 106   format(f7.0)
       ne=7
     
       else 
       write(6,*) 'error of file name'
       stop

       end if

       return
       end
