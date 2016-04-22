 Module M_MPI_WaveGen_3D

use M_General, only: nx,ny,nz,pi,x,y,z,phi,landa,gz,dt,tpstar,lx,ly,Lz,Froude


implicit none


real(8),dimension ( 0:1,1:ny-1,1:nz-1) :: uWa
real(8),dimension (-1:0,1:ny-1,0:nz+1) :: wWa
real(8),dimension(:),Allocatable       :: xgage,ygage,zgage
real(8),dimension(1:ny-1)              :: ZGageGen
real(8)                                :: hait,wavegen,UCurrent,HH


contains


subroutine Wave_Gen_Constant_Ini()
implicit none 

include "Par_Constant_Wave.txt"

return

end subroutine 


subroutine Regularwave ()

implicit none 

real(8)  :: HH1,period,kwave,wwave,LRamp
integer  :: k,j





kwave=2*pi/(landa)
period=1.0/(  1/(2*pi)* sqrt(kwave*abs(gz)*tanh(kwave*hait)) )
wwave=2*pi/(period)


if ( ((tpstar-1)*dt).lt.(0.1*period) ) then 
HH1=((tpstar-1)*dt)/(0.1*period)*HH
else 
HH1=HH
end if

print*,"wave height=",HH1,"tpstar=",tpstar

!do k=1, nz-1 
!print*, z(k)
!end do 


    ZGageGen(:)=0 
    do j=1,ny-1 
      do k=1,nz-1
       if ((0.5*( phi(1,j,k+1)+phi(0,j,k+1) ).ge.0 )) then 
          !! may be more precsiness can be done above !!      
             Zgagegen(j)=z(k)   +  &
          &   ( z(k+1)-z(k) )/ ( 0.5*( phi(1,j,k+1)+phi(0,j,k+1) ) &
          &     -0.5*( phi(1,j,k)+phi(0,j,k) ) ) *( -0.5*( phi(1,j,k)+phi(0,j,k) ) ) 
        
          exit 
        end if 
      end do 
     end do 
    

UWa=0
do k=1,nz-1
do j=1,ny-1 
       

 uWa(1,j,k)=UCurrent+&
 & HH1/2* wwave *cos(-wwave*(tpstar*dt))*cosH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait)!+ &
    !& 3/4*((pi*HH1)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))* &
    !& cosH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )
 
 uWa(0,j,k)=UCurrent+&
 & HH1/2* wwave*cos(kwave*(0.5*(x(0)+x(-1)))-wwave*(tpstar*dt))*cosH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait) !+ &
    !& 3/4*((pi*HH1)**2)/(period*landa)*cos(2*kwave*(0.5*(x(0)+x(-1)))-2*wwave*(tpstar*dt))* &
    !& cosH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )
                        
      end do
end do 




WWa=0
 do k=0,nz+1
 do j=1,ny-1   

 wWa(0,j,k)=HH1/2*wwave*  &
& sin(kwave*x(0)-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait) !+ &
     !& 3/4*((pi*HH1)**2)/(period*landa)*sin(2*kwave*x(0)-2*wwave*(tpstar*dt))* &
     !& sinH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )

 wWa(-1,j,k)=HH1/2*wwave* &
& sin(kwave*x(-1)-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait) !+ &
     !& 3/4*((pi*HH1)**2)/(period*landa)*sin(2*kwave*x(-1)-2*wwave*(tpstar*dt))* &
     !& sinH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )
 
 
 end do   
 end do 



!! ramp for air velocity to reach zero within half of the distance from top
do j=1,ny-1
 
 LRamp=0.5*(Lz-ZGageGen(j))

  do k=1,nz-1
       
   if (z(k).gt.ZGageGen(j)) then 

    uWa(1,j,k)=uWa(1,j,k)* 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )  
    uWa(0,j,k)=uWa(0,j,k)* 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )  
 
   end if 
  end do 
end do 



 do j=1,ny-1 
  
  LRamp=0.5*(Lz-ZGageGen(j))
  do k=0,nz+1

   if (z(k).gt.ZGageGen(j)) then 

     wWa(0,j,k) =wWa(0,j,k) * 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )
     wWa(-1,j,k)=wWa(-1,j,k)* 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )
   end if
 
  end do 

 end do 




return

end subroutine 
 








subroutine Wavegage(s,e,comm1d)

   IMPLICIT NONE
   integer,intent(in)                                             :: s, e
   integer,intent(in)                                             :: comm1d
   integer igage,jgage,i,j,k,kk




   do kk=1, size(xgage)  


  
     igage=1000
     zgage_cpu(kk)=0
     if (Xgage(kk).ge.x(s).AND.Xgage(kk).le.x(e)) then

 
       do  i=s,e
         if (x(i).le.xgage(kk).AND.x(i+1).ge.xgage(kk)) then 
           igage=i
           exit 
         end if 
       end do
       do  j=1,ny-1
         if (y(j).le.ygage(kk).AND.y(j+1).ge.ygage(kk)) then 
           jgage=j
           exit 
         end if 
       end do
 
       zgage(kk)=0
       do k=1,nz-1
         if (0.5*( phi(igage,jGage,k)+phi(igage-1,jGage,k) ).le.0.AND.0.5*( phi(igage,jGage,k+1)+phi(igage-1,jGage,k+1) ).ge.0 ) then 
         !! may be more precsiness can be done above !!      
           zgage_cpu(kk)=z(k)+ &
           &   ( z(k+1)-z(k) )/ ( 0.5*( phi(igage,jGage,k+1)+phi(igage-1,jGage,k+1) ) &
           &     -0.5*( phi(igage,jGage,k)+phi(igage-1,jGage,k) ) ) *( -0.5*( phi(igage,jGage,k)+phi(igage-1,jGage,k) ) ) 
           exit 
         end if 
       end do 

       if (igage.eq.1000) then 
         print *, "error in wave gage"
       end if 
   

    end if
  



  end do 

  if (Parallel) then
      call MPI_Allreduce( ZGage_CPU ,Zgage,size(xgage),MPI_Double_Precision, MPI_Sum, comm1d, ierr )
  end if

  write(75,2020) tp*dt,zgage(1),zgage(2),zgage(3),zgage(4),zgage(5),zgage(6),zgage(7), &
  &                      zgage(8),zgage(9),zgage(10),zgage(11),zgage(12),zgage(13),zgage(14)
  2020  format (15(1x,e15.7))
  call flush (75)


             
return        
  end subroutine  




end module 





