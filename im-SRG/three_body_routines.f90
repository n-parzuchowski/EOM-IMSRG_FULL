module three_body_routines 
  use basic_IMSRG
  
  
  contains


real(8) function overlap_3b(a1,b1,c1,Jab1,Tab1,a2,b2,c2,Jab2,Tab2,jtot,ttot,jbas) 
  ! INCLUDES ANTI-SYMMETRY FACTOR
  implicit none 
  
  type(spd) :: jbas
  integer :: a1,b1,c1,a2,b2,c2
  integer :: Jab1,Jab2,Tab1,Tab2,jtot,ttot
  integer :: ja,jb,jc
  
  if ( a1 == a2 ) then 
     
     if ( b1 == b2) then 
        
        
        
        if (c1 == c2) then 
           
           if ( Jab1 == Jab2 ) then 
              if (Tab1 == Tab2) then 
                 overlap_3b = 1.d0 
              else 
                 overlap_3b = 0.d0 
              end if
           else
              overlap_3b=0.d0 
           end if 
           
        else
           
           overlap_3b = 0.d0
        
        end if 
  
     else if (b1 == c2 )  then 
        
        if (c1 == b2 ) then 
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = (-1) **(( jc + jb + Jab1 + Jab2 +Tab1 + Tab2)/2 ) &
                * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(jb,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2) 
           
        else 
           
           overlap_3b = 0.d0 
        end if 
     else 
        overlap_3b = 0.d0 
     end if 
  else if (a1 == b2) then 
     
     if (b1 == a2 ) then 
        
        if ( c1 == c2)  then 
           
           if (Jab1 == Jab2) then
              if (Tab1 == Tab2) then 
                 
                 ja = jbas%jj(a1) 
                 jb = jbas%jj(b1) 
           
                 overlap_3b =  (-1)**((ja+jb-Jab1-Tab1)/2) 
              else 
                 overlap_3b = 0.d0 
              end if 
           else 
              overlap_3b = 0.d0 
           end if 
        else
           overlap_3b = 0.d0 
        end if 
        
     else if (b1 == c2) then 
        if (c1 == a2 ) then 
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = (-1) **(( jc + jb + Jab1 + Jab2 +Tab1 + Tab2)/2 ) &
                * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(jb,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2) &
                * (-1) ** (( ja + jc  - Jab2 -Tab2 )/2) 
           
        else 
           overlap_3b = 0.d0 
        end if 
     else 
        overlap_3b = 0.d0 
     end if 
     
  else if (a1 == c2)  then 
     
     if (b1 == b2) then 
        
        if ( c1 == a2) then
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = -1*sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(ja,jb,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)
        else
           overlap_3b = 0.d0
        end if 
        
     else if (b1 == a2) then 
        if ( c1 == b2) then 
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = -1*sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(ja,jb,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2) &
                * (-1)**((jb +jc-Jab2-Tab2)/2)
           
        else 
           overlap_3b = 0.d0 
        end if 
     end if 
  else 
     overlap_3b = 0.d0
  end if 
  
end function

end module
        
           
           
     
     
        
        
  
  
  
  
  



