module M_MPI_Mesh_3D
  use M_MPI_General_3D, only:x,y,z,hx,hy,hz,nx,ny,nz,lx,ly,lz,Chv,Bp,parallel
  use M_MPI_Exch_3D

  implicit none 
  integer,private :: meshgen
  real(8),private :: Lx1,Lx2,Lx3,Lx4,Lz1,Lz2,Lz3,Ly1,Ly2,Ly3
  real(8),private :: hxx1,hxx2,hxx3,hxx4,hyy1,hyy2,hyy3,hzz1,hzz2,hzz3
  real(8),private :: Amesh12,Amesh14,Amesh21,Amesh23,Amesh31,Amesh32,Amesh33,betam12,betam14,betam21,betam23
  real(8),private :: betam31,betam32,betam33,Dmesh12,Dmesh14,Dmesh21,Dmesh23,Dmesh31,Dmesh32,Dmesh33
  integer,private :: nx1,nx2,nx3,ny1,ny2,nz1,nz2





  contains 
    subroutine meshgenerator(s,e,MyRank,comm1d,nbr_right,nbr_left)

    integer,intent(in)                                             ::    s,e
    integer,intent(in)                                             ::    MyRank,comm1d,nbr_right,nbr_left

    double precision,dimension(:),allocatable                      ::    XE
    dOUBLE PRECISION,dimension(0:ny-1)                             ::    YE
    dOUBLE PRECISION,dimension(0:nz-1)                             ::    zE
    dOUBLE PRECISION,dimension(1:nx-1)                             ::    XTo

    dOUBLE PRECISION                                               ::    SPoint_1,SPoint_2,SPoint_3
    dOUBLE PRECISION                                               ::    SPoint_1_cpu,SPoint_2_cpu,SPoint_3_cpu
    integer                                                        ::    i,j,k
    integer                                                        ::    ierr


    include 'Par_MPI_Mesh_3D.txt'
    if (meshgen.eq.1) then 





      do i=s-2,e+2
        x(i)=real(i)*Lx/(nx-1)-0.5*lx/(nx-1)
      end do
      hx(:)=Lx/(nx-1)

      do j=-1,ny+1
        y(j)=real(j)*Ly/(ny-1)-0.5*ly/(ny-1)
      end do
      hy(:)=Ly/(ny-1)

      do k=-1,nz+1
        z(k)=real(k)*Lz/(nz-1)-0.5*lz/(nz-1)
      end do
      hz(:)=Lz/(nz-1)


   


   else





    allocate(XE(s-1:e))
    SPoint_1_CPU=0
    SPoint_2_CPU=0
    SPoint_3_CPU=0

    nx2=nx2+nx1
    nx3=nx3+nx2
    nz2=nz2+nz1    
    ny2=ny2+ny1
    
    hxx1=1.0/dble(nx1-1) ; hxx2=1.0/dble(nx2-nx1) ;  hxx3=1.0/dble(nx3-nx2) ; hxx4=1.0/dble(nx-nx3)
    hyy1=1.0/dble(ny1-1) ; hyy2=1.0/dble(ny2-ny1) ;  hyy3=1.0/dble(ny-ny2)  
    hzz1=1.0/dble(nz1-1) ; hzz2=1.0/dble(nz2-nz1) ;  hzz3=1.0/dble(nz-nz2)

    Amesh12=1/(2*betam12)*log(  (  1+( exp(betam12)-1 )*Dmesh12/Lx2 )/(  1+( exp(-betam12)-1 )*Dmesh12 /Lx2   )    ) 
    Amesh14=1/(2*betam14)*log(  (  1+( exp(betam14)-1 )*Dmesh14/Lx4 )/(  1+( exp(-betam14)-1 )*Dmesh14 /Lx4   )    ) 


    Amesh21=1/(2*betam21)*log(  (  1+( exp(betam21)-1 )*Dmesh21/Ly1 )/(  1+( exp(-betam21)-1 )*Dmesh21 /Ly1   )    ) 
    Amesh23=1/(2*betam23)*log(  (  1+( exp(betam23)-1 )*Dmesh23/Ly3 )/(  1+( exp(-betam23)-1 )*Dmesh23 /Ly3   )    )

    Amesh31=1/(2*betam31)*log(  (  1+( exp(betam31)-1 )*Dmesh31/Lz1 )/(  1+( exp(-betam31)-1 )*Dmesh31 /Lz1   )    ) 
    Amesh32=1/(2*betam32)*log(  (  1+( exp(betam32)-1 )*Dmesh32/Lz2 )/(  1+( exp(-betam32)-1 )*Dmesh32 /Lz2   )    )  
    Amesh33=1/(2*betam33)*log(  (  1+( exp(betam33)-1 )*Dmesh33/Lz3 )/(  1+( exp(-betam33)-1 )*Dmesh33 /Lz3   )    )  

    
 
     if ( (s.ge.0.AND.s.le.(nx1-1)).OR.(e.ge.0.AND.e.le.(nx1-1) ).OR.(s.le.0.AND.e.ge.(nx1-1)) )then
       do i=max(s,0),min(nx1-1,e)
         xE(i)=hxx1*Lx1*dble(i)
       end do
     end if  
    
     if ((nx1-1).ge.s.AND.(nx1-1).le.e) then
       SPoint_1_CPU=xE(nx1-1)
     end if
     if (parallel) then
       call MPI_Allreduce( SPoint_1_CPU ,SPoint_1  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
     endif


     if ( (s.ge.nx1.AND.s.le.(nx2-1)).OR.(e.ge.nx1.AND.e.le.(nx2-1)).OR.(s.le.nx1.AND.e.ge.(nx2-1)) )then
       do i=max(nx1,s),min(nx2-1,e)
         xE(i)=SPoint_1+Dmesh12* (  1+ (  sinh ( betam12*(dble(i-(nx1-1))*hxx2-Amesh12) )  )/sinh(betam12*Amesh12)  )
       end do
     end if 

     if ((nx2-1).ge.s.AND.(nx2-1).le.e) then
       SPoint_2_CPU=xE(nx2-1)
     end if
     if (parallel) then
       call MPI_Allreduce( SPoint_2_CPU ,SPoint_2  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
     end if


     if ( (s.ge.nx2.AND.s.le.(nx3-1)).OR.(e.ge.nx2.AND.e.le.(nx3-1)).OR.(s.le.nx2.AND.e.ge.(nx3-1)) )then
       do i=max(nx2,s),min(nx3-1,e) 
         xE(i)=SPoint_2+hxx3*Lx3*dble( i-(nx2-1) )
       end do 
     end if
 
     if ((nx3-1).ge.s.AND.(nx3-1).le.e) then
       SPoint_3_CPU=xE(nx3-1)
     end if
     if (parallel) then
       call MPI_Allreduce( SPoint_3_CPU ,SPoint_3  ,1 ,MPI_Double_Precision, MPI_Sum, comm1d, ierr )
     end if


     if ( (s.ge.nx3.AND.s.le.(nx-1)).OR.(e.ge.nx3.AND.e.le.(nx-1)).OR.(s.le.nx3.AND.e.ge.(nx-1)) )then
       do i=max(nx3,s),min(nx-1,e) 
         xE(i)=SPoint_3+Dmesh14* (  1+ ( sinh ( betam14*(dble(i-(nx3-1)) *hxx4-Amesh14) )  )/sinh(betam14*Amesh14)  )
       end do
     end if  
     
     
      
     yE(0)=0
     do j=1,ny1-1 
       yE(j)=          Dmesh21* (  1+ (  sinh ( betam21*(dble(j)        *hyy1-Amesh21) )  )/sinh(betam21*Amesh21)  )
     end do 
     do j=ny1,ny2-1
       yE(j)=yE(ny1-1)+hyy2*Ly2*dble(j-(ny1-1))
     end do 
     do j=ny2,ny-1
       yE(j)=yE(ny2-1)+Dmesh23* (  1+ (  sinh ( betam23*(dble(j-(ny2-1))*hyy3-Amesh23) )  )/sinh(betam23*Amesh23)  )
     end do 
     
     

     ZE(0)=0
     do k=1,nz1-1
       zE(k)=         Dmesh31* (  1+ (   sinh ( betam31*(dble(k)        *hzz1-Amesh31) )  )/sinh(betam31*Amesh31)  )
     end do 
     do k=nz1,nz2-1
       zE(k)=zE(nz1-1)+Dmesh32* (  1+ (   sinh ( betam32*(dble(k-(nz1-1))*hzz2-Amesh32) )  )/sinh(betam32*Amesh32)  )
     end do 
     do k=nz2,nz-1
       zE(k)=zE(nz2-1)+Dmesh33* (  1+ (   sinh ( betam33*(dble(k-(nz2-1))*hzz3-Amesh33) )  )/sinh(betam33*Amesh33)  )
     end do 

   

     if (Parallel) then
       call XE_Exch(XE, s, e, comm1d, nbr_left, nbr_right)
     end if

     do i=s,e
       x(i)=XE(i-1)+0.5*(XE(i)-XE(i-1))
       hx(i)=XE(i)-XE(i-1)
     end do
     if (s==1) then 
       x(-1)=-3*x(1)
       x(0)=-x(1)
       hx(0)=hx(1)
       hx(-1)=hx(0)
     end if
 
     if (e==(nx-1)) then 
       x(nx)=x(nx-1)+ ( x(nx-1)-x(nx-2) )
       x(nx+1)=x(nx)+ ( x(nx)-x(nx-1) )
       hx(nx)=hx(nx-1)
       hx(nx+1)=hx(nx)
     end if 



    do j=1,ny-1
      y(j)=YE(j-1)+0.5*(YE(j)-YE(j-1))
      hy(j)=YE(j)-YE(j-1)
    end do   
    y(-1)=-3*y(1)
    y(0)=-y(1)
    y(ny)=y(ny-1)+ ( y(ny-1)-y(ny-2) )
    y(ny+1)=y(ny)+ ( y(ny)-y(ny-1) )

    hy(0)=hy(1)
    hy(-1)=hy(0)
    hy(ny)=hy(ny-1)
    hy(ny+1)=hy(ny)



    do k=1,nz-1
      z(k)=zE(k-1)+0.5*(zE(k)-zE(k-1))
      hz(k)=zE(k)-zE(k-1)
    end do   
    z(-1)=-3*z(1)
    z(0)=-z(1)
    z(nz)=z(nz-1)+ ( z(nz-1)-z(nz-2) )
    z(nz+1)=z(nz)+ ( z(nz)-z(nz-1) )

    hz(0)=hz(1)
    hz(-1)=hz(0)
    hz(nz)=hz(nz-1)
    hz(nz+1)=hz(nz)


    if (Parallel) then 
      call Xhx_Exch(X,  s, e, comm1d, nbr_left, nbr_right)
      call Xhx_Exch(hx, s, e, comm1d, nbr_left, nbr_right)
    end if

 



 end if 




 call MPI_GATHER (X(s),e-s+1,MPI_dOUBLE_PRECISION ,XTo(1),e-s+1,MPI_dOUBLE_PRECISION,0,comm1d,ierr )  
 if (MyRank==0) then 


   OPEN(205,file='meshx.CSV')
   OPEN(215,file='meshy.CSV')
   OPEN(225,file='meshz.CSV')
   OPEN(235,file='grideval2DXZ.plt')

   do i=1,nx-1
     write(205,*) i,',',xTo(i)
   end do 
   call flush (205)

   do j=1,ny-1
     write(215,*) j,',',y(j)
   end do
   call flush (215)

   do k=1,nz-1
     write(225,*) k,',',z(k)
   end do 
   call flush (225)

   write(235,*) 'zone i=',nx-1,' k=',nz-1
   do k=1,nz-1 
     Do i=1,nx-1
      write(235,350) xTo(i),z(k)
     end do 
   end do
   350 format (3(1x,e15.7))
   call flush (235)

 end if 



 return 
  end subroutine 
  
  


     
end module 






