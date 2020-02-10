      subroutine simpson(XY,npontos,H,S)
      
      Integer Npontos,I
      Real*8 XY(npontos),H,S
      
      S=0.0
      do I=1,Npontos-2,2
        S=S+XY(I)+4.0*XY(i+1)+XY(I+2)
      enddo
      S=S*H/3.0
      Return
      end subroutine simpson
