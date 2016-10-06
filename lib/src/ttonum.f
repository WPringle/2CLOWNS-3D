       subroutine ttonum(t,num,ns,ne)
       real*8 t
       character*255 num
       character*255 num1

       ns=1
       ne=5

       if(t.lt.0.009) then
CC -1 keta
       write(num1,108) nint(t*1000.)
 108   format(I1)
       num='0.00'//num1(1:1)
       else if(t.lt.0.099) then
CC -1 keta
       write(num1,107) nint(t*1000.)
 107   format(I2)
       num='0.0'//num1(1:2)
       else if(t.lt.0.999) then
CC -1 keta
       write(num1,103) nint(t*1000.)
 103   format(I3)
       num='0.'//num1(1:3)
       else if(t.lt.9.99) then
CC 1 keta
       write(num1,109) t

       num=num1(2:6)
 109   format(f6.3)

       else if(t.lt.99.9) then
CC 2 keta
       write(num,102) t
 102   format(f5.2)
       else if(t.lt.999.9) then
CC 3 keta
       write(num,101) t
 101   format(f5.1)
       else if(t.lt.9999.9) then
CC 3 keta
       write(num,104) t
 104   format(f6.1)
       ne=6
       else if(t.lt.99999.9) then
CC 3 keta
       write(num,105) t
 105   format(f7.1)
       ne=7
     
       else if(t.lt.999999.9) then
CC 3 keta
       write(num,106) t
 106   format(f8.1)
       ne=8
     
       else 
       write(6,*) 'error of file name'
       stop

       end if

       return
       end
