
      function PSHU(PP,PS,DXS,DXP,FS)
      implicit none
      doubleprecision::PSHU
      doubleprecision,intent(in)::PP,PS,DXS,DXP,FS
  !    doubleprecision::DL,DS
 !     REAL*8 dl,dxs,dxp,ds,fs,pshu,ps,pp
   !     DL= (DXS+DXP)*0.5d0   
   !     DS=DXP*0.5d0+DXS*FS    
        pshu=( PS-PP )/(DXP*0.5d0+DXS*FS )*((DXS+DXP)*0.5d0 )+PP       
      END FUNCTION
