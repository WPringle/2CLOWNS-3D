function u3ss(a11,a21,a31,y0,y1,y2,y3)
       implicit real*8 (a-h,o-z)
       a12=a11**2
       a22=a21**2
       a32=a31**2
       a13=a11**3
       a23=a21**3
       a33=a31**3
       detem=a11*a22*a33&
            +a21*a32*a13&
            +a31*a12*a23&
            -a11*a32*a23&
            -a21*a12*a33&
            -a31*a22*a13
      if(detem.ne.0.0d0) then
  
      b11=a22*a33-a32*a23
      b12=a32*a13-a12*a33
      b13=a12*a23-a22*a13
      u3ss1=(b11*(y0-y1)+b12*(y2-y1)+b13*(y3-y1))/detem
      as1=(y1-y2)/(0.0d0-a21)
      as2=(y1-y3)/(0.0d0-a31)
      if(as1*as2.ge.0.0d0) then
        if(as1*u3ss1.le.0.0d0) then
          u3ss=as1
        else
          u3ss=u3ss1
        end if
      else 
          u3ss=u3ss1
      end if
      else
      u3ss=-1.d10       
      end if
end function 
