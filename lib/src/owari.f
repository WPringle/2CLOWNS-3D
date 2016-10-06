       subroutine owari(text,nx)

       character*255 text

       do 100 n=255,1,-1
       
c      write(6,*) text(n:n),'n=',n
        if(text(n:n).ne.' ') then
          nx=n
          return
        end if
        
 100   continue
       
       write(6,*) 'owaranai!!!'
        return
       end

       subroutine hajime(text,nx)

       character*255 text

       do 100 n=1,255
       
c      write(6,*) text(n:n),'n=',n
        if(text(n:n).ne.' ') then
          nx=n
          return
        end if
        
 100   continue
       
       write(6,*) 'hajimaranai!!!'
        return
        end

       subroutine sura(text,nx)
       character*255 text

       do 110 n=1,255
       
        if(text(n:n).eq.'/') then
       !write(6,*) text(n:n),'n=',n
          nx=n
        end if
        
 110   continue
       
       !write(6,*) 'suranai!!!'
        return
       end
  
       subroutine numplu(number,text,text1)
       character*255 suuji,text,text1
       if(number.lt.1000) then
          write(6,*) 'number lt 1000'
          return
       else if(number.ge.10000) then
          write(6,*) 'number ge 10000'
          return
       end if
       write(suuji,405) number
 405   format(i4)
        call owari(text,mm)
        text1=text(1:mm)//suuji
        return
       end

       function nyMOD(N,MWR)
       if(mwr.eq.0) then 
           nymod=1 
       else
           nymod=MOD(N,MWR)
       end if
       return
       end

       subroutine pritime(idate)
       implicit real*8 (a-h,o-z) 
       CHARACTER*8 idate
#ifdef SCS
       call DATE(idate) 
       call TIME(itime) 
       itime=int(itime/1000)
       ihour=int(itime/3600)
       imint=int((itime-ihour*3600)/60)
       isec= itime-ihour*3600-imint*60
       write(6,392) idate,ihour,imint,isec
 392   format(1h ,A8,1x,i2.2,':',i2.2,':',i2.2)
#else
      CHARACTER*10 itime
      CALL DATE_AND_TIME(idate,itime)
      write(6,*) idate(1:4)//'-'//idate(5:6)//'-'//idate(7:8)//' , '
     1          //itime(1:2)//':'//itime(3:4)//':'//itime(5:6)
#endif
!       call flush(6)
       return
       end
