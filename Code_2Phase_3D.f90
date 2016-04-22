
program Code_2Phase_3D
  use mpi
  use M_MPI_General_3D
  use M_MPI_Exch_3D
  use M_MPI_Mesh_3D
  use M_MPI_WaveGen_3D
  use M_MPI_Solver_3D
  use M_MPI_FrSrf_3D
  implicit none


  Double precision,dimension (:,:,:),allocatable       ::   advectu,advectv,advectw
  Double precision,dimension (:,:,:),allocatable       ::   Tx,Ty,Tz   
  Double precision,dimension (:,:,:),allocatable       ::   divv
  Double precision,dimension (:,:,:),allocatable       ::   ro,p,miuv

  Double precision,dimension (:,:,:),allocatable       ::   dpux,dmux,dpuy,dmuy,dpuz,dmuz
  Double precision,dimension (:,:,:),allocatable       ::   dpvx,dmvx,dpvy,dmvy,dpvz,dmvz
  Double precision,dimension (:,:,:),allocatable       ::   dpwx,dmwx,dpwy,dmwy,dpwz,dmwz

  Double precision,dimension (:,:,:),allocatable       ::   Apx,Amx,Apy,Amy,Apz,Amz,AP,Q
  Double precision                                     ::   div,divmax,maxp
  Double precision                                     ::   StartTime,EndTime
  integer                                     ::   s,e
  integer                                     ::   i,j,k,tp,Count_Plot
  
  integer                                     ::   comm1d,nbr_left,nbr_right,myid,MyRank,numprocs,ierr
  integer,dimension(1)                        ::   DIMS
  logical,dimension(1)                        ::   PERIODS



  if (Parallel) then 

    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

    print *, "my id is",myid,"of",numprocs
    PERIODS=.FALSE.
    DIMS=numprocs

    call MPI_CART_CREATE(MPI_COMM_WORLD,1,DIMS, PERIODS,.true., comm1d, ierr)
    call MPI_COMM_RANK( comm1d, MyRank, ierr )
    call MPI_Cart_shift(comm1d,0,1,nbr_left,nbr_right,ierr)



   
  else

    numprocs=1
    MyRank=0
    myid=0

  endif

  if (myid.eq.0) then
    call CPU_TIME(StartTime)
    print*,"Start time is", StartTime
  end if



  s=1+MyRank*((nx-1)/numprocs)
  e=s+((nx-1)/numprocs)-1
  print*,s,e,"MyRank",MyRank,"Rank on my left=",nbr_left,"rank on my right",nbr_right
  call flush (6)

   

   count_Plot=10000
   call General_Constant_Ini(s,e)
   call meshgenerator(s,e,MyRank,comm1d,nbr_right,nbr_left)
   call Wave_Gen_Constant_Ini()  
      !if (Parallel) then
   !    call Xhx_Exch(X,  s, e, comm1d, nbr_left, nbr_right)
   !    call Yhx_Exch(hx, s, e, comm1d, nbr_left, nbr_right)
   !endif


   allocate(  divv(s:e,1:ny-1,1:nz-1) )
   allocate(  p(s-1:e+1,0:ny,0:nz),ro(s-1:e+1,0:ny,0:nz),miuv(s-1:e+1,0:ny,0:nz) )
   allocate(  Tx(s:e,1:ny-1,1:nz-1),Ty(s:e,1:ny-1,1:nz-1),Tz(s:e,1:ny-1,1:nz-1) )
   allocate(  dpux(s-1:e+1,0:ny,0:nz),dmux(s-1:e+1,0:ny,0:nz),dpuy(s-1:e+1,0:ny,0:nz),dmuy(s-1:e+1,0:ny,0:nz),dpuz(s-1:e+1,0:ny,0:nz),dmuz(s-1:e+1,0:ny,0:nz)  )
   allocate(  dpvx(s-1:e+1,0:ny,0:nz),dmvx(s-1:e+1,0:ny,0:nz),dpvy(s-1:e+1,0:ny,0:nz),dmvy(s-1:e+1,0:ny,0:nz),dpvz(s-1:e+1,0:ny,0:nz),dmvz(s-1:e+1,0:ny,0:nz)  )
   allocate(  dpwx(s-1:e+1,0:ny,0:nz),dmwx(s-1:e+1,0:ny,0:nz),dpwy(s-1:e+1,0:ny,0:nz),dmwy(s-1:e+1,0:ny,0:nz),dpwz(s-1:e+1,0:ny,0:nz),dmwz(s-1:e+1,0:ny,0:nz)  )
   allocate(  advectu(s:e,1:ny-1,1:nz-1),advectv(s:e,1:ny-1,1:nz-1),advectw(s:e,1:ny-1,1:nz-1)  )
   allocate(  Amx(s:e,1:ny-1,1:nz-1),Apx(s:e,1:ny-1,1:nz-1),Amy(s:e,1:ny-1,1:nz-1),Apy(s:e,1:ny-1,1:nz-1),Amz(s:e,1:ny-1,1:nz-1),Apz(s:e,1:ny-1,1:nz-1),  &
          &  Ap (s:e,1:ny-1,1:nz-1),Q  (s:e,1:ny-1,1:nz-1)  )

   call Fluid_Dynamic_Ini(p,s,e)
   call Boundarycond(s,e)   



   call Levelset_Ini(MyRank,comm1d,s,e)

   if (MyRank==0) then 
     OPEN(35,file='Result.plt')
     OPEN(45,file='CavityU.plt')
     OPEN(55,file='CavityV.plt')
     OPEN(65,file='FreeSurface.CSV')
     OPEN(75,file='WaveGage.plt')
     write(35,*) 'variables="x","y","z","u","v","w","p","ro","phi"'
     write(75,'(A200)') 'variables="t","Gage1","Gage2","Gage3","Gage4","Gage5","Gage6","Gage7",&
                &        "Gage1b","Gage2b","Gage3b","Gage4b","Gage5b","Gage6b","Gage7b"'

     write(65,'(A50)') "time,GaugeZ"
     !write(35,*) 'variables="x","y","u","v","w","p","ro"'
   endif 

   do tp=1,tstepE !time step      
      
     print*,"time step=",tp
     call flush (6)

           call Proprety(ro,rodrop,roAir,s,e)
           call Proprety(miuv,miudrop,miuAir,s,e) 

           call Advection_Ini(dmux,dpux,dmuy,dpuy,dmuz,dpuz,    &
                            & dmvx,dpvx,dmvy,dpvy,dmvz,dpvz,    &
                            & dmwx,dpwx,dmwy,dpwy,dmwz,dpwz,s,e )
                      
          
           if (Parallel) then
            
             call Dy_A_Exch(Dmux,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
             call Dy_A_Exch(Dpux,nx,ny,nz, s, e, comm1d, nbr_left, nbr_right)
             call Dy_A_Exch(Dmvx,nx,ny,nz, s, e, comm1d, nbr_left, nbr_right)
             call Dy_A_Exch(Dpvx,nx,ny,nz, s, e, comm1d, nbr_left, nbr_right)
             call Dy_A_Exch(Dmwx,nx,ny,nz, s, e, comm1d, nbr_left, nbr_right)
             call Dy_A_Exch(Dpwx,nx,ny,nz, s, e, comm1d, nbr_left, nbr_right)
           end if
          
                      call Advection(dmux,dpux,dmuy,dpuy,dmuz,dpuz, &
                        & dmvx,dpvx,dmvy,dpvy,dmvz,dpvz, &
                        & dmwx,dpwx,dmwy,dpwy,dmwz,dpwz, &
                        & advectu,advectv,advectw,s,e    )

           call Viscoset(miuv,ro,Tx,Ty,Tz,s,e)
           
           Do i=s,e 
             Do j=1, ny-1 
               do k=1, nz-1
                 U(i,j,k)=0 !U(i,j,k)+( advectu(i,j,k)+Tx(i,j,k) )*dt
               end do 
             end do 
           end do  

           do j=1,ny-1 
             do i=s,e 
               do k=1, nz-1
                 V(i,j,k)=V(i,j,k)+( advectv(i,j,k)+Ty(i,j,k) )*dt
               end do 
             end do 
           end do  

           do k=1, nz-1
             do j=1,ny-1 
               do i=s,e 
                 W(i,j,k)=W(i,j,k)+( advectW(i,j,k)+Tz(i,j,k)+gz )*dt  
               end do 
             end do 
           end do  

          
           call BoundCond_UVW_Poisson(s,e) !! The velocities are defined and are not based on the 
          !if (s==1) then
          !  w(3,-1,0)=1.0
          !   w(3,3,0)=10.0
          !   w(3,6,5)=20.0
          ! end if
           
          ! if (e==(nx-1)) then
          !   w(4,0,-1)=2.5
          !   w(4,3,0)=3.5
          !   w(4,6,7)=4.5 
          ! end if

           if (Parallel) then
             call UVPhi_Exch(U,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
             call UVPhi_Exch(V,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
             call UVPhi_Exch(W,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
           end if

          ! if (e==(nx-1)) then
           !if (s==1) then
          !   print*,w(3,-1,0),s,e
          !   print*,w(3,3,0),s,e
          !   print*,w(3,6,5),s,e
          ! end if
           
          ! if (s==1) then
          !   print*,w(4,0,-1),s,e
          !   print*,w(4,3,0),s,e
          !   print*,w(4,6,7),s,e
          ! end if

           

           call Poisson_COF(ro,dt,Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,s,e) 

           if ( Parallel ) then
             call Poisson_Solver_Parallel_3D(Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,pdif,beta,P,s,e,comm1d,nbr_left,nbr_right)
           else
             call Poisson (Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,beta,pdif,p,s,e)
           end if !!!
                    

           do i=s,e
             do j=1,ny-1 
               do k=1, nz-1
                 U(i,j,k)=0 !U(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )   )*( p(i,j,k)-p(i-1,j,k) )*dt/( x(i)-x(i-1)  )
               end do 
             end do 
           end do 

           do i=s,e 
             do j=1,ny-1 
               do k=1, nz-1
                 V(i,j,k)=V(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )   )*( p(i,j,k)-p(i,j-1,k) )*dt/( y(j)-y(j-1)  )
               end do 
             end do 
           end do 

           do i=s,e 
             do j=1,ny-1 
               do k=1, nz-1
                 W(i,j,k)=W(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )   )*( p(i,j,k)-p(i,j,k-1) )*dt/( z(k)-z(k-1)  )
               end do 
             end do 
           end do
 
           call Boundarycond(s,e)    !! it is only for divergence in the next line other wise it could be merged with after solid if the solid is not on the boundaries

           if (Parallel) then
             call UVPhi_Exch(U,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
             call UVPhi_Exch(V,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
             call UVPhi_Exch(W,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
           end if

           divmax=0
           do i=s,e 
             do j=1,ny-1 
               do k=1,nz-1

                 div= (U(i+1,j,k)-U(i,j,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
                 &(v(i,j+1,k)-v(i,j,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+&
                 &(W(i,j,k+1)-W(i,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) 
    
                 if(abs(div).gt.divmax)then
                   divmax=abs(div)
                 end if 

               end do 
             end do 
           end do 
           print*,'DIVERGENCE before correction=',DIVMAX

           call levelset (dt,s,e,comm1d, nbr_left,  nbr_right)
           call reinitialize(s,e,comm1d, nbr_left,  nbr_right)

           if (Parallel) then
             call PrintData(p,ro,plot,tp,TStepE,dt,MyRank,numprocs,s,e,comm1d)
             call Wavegage(s,e,comm1d)
           endif 
           !count_Plot=count_Plot+1
           !if (count_Plot.ge.plot)then

           !  count_Plot=0
           !  print*, "Write data"
           !  write(35,*) 'ZONE T="',tp,'" SOLUTIONTIME=',tp*dt,'i=',nx-1,' j=',ny-1,' k=',nz-1 
           !  do k=1,nz-1 
           !    do j=1,ny-1 
           !      do i=s,e
           !        write(35,122) x(i),y(j),z(k),0.5*(u(i,j,k)+u(i+1,j,k)),0.5*(v(i,j,k)+v(i,j+1,k)), &
           !        &            0.5*(w(i,j,k)+w(i,j,k+1)),p(i,j,k)
           !      end do 
           !    end do
           !  end do
           !  122 format (7(1x,e15.7))

           !endif 




   end do !time step


   if (myid.eq.0) then
     write(*,*)'end'

     call CPU_TIME(EndTime)
     print*,"Simulation time is",EndTime-StartTime,"seconds"
   end if


   if (Parallel) then
     call MPI_FINALIZE(ierr)
   end if



end program 




!!!!!!!!!!!!!!!!!!!!!!! subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine BoundaryCond(s,e)
    use M_MPI_General_3D, only: u,v,w,nx,ny,nz
    use M_MPI_WaveGen_3D

    implicit none 
    integer,intent(in)                  ::   s,e

    integer i,j,k


    if (wavegen) then
      call regularwave()
    end if
 
    if (s==1) then
      do j=1, ny-1 
        do k=1, nz-1  
          U(1,j,k)=0       !U(2,j,k)      
          U(0,j,k)=-U(2,j,k)
          U(-1,j,k)=-U(3,j,k)
        end do 
      end do
     end if 
   
     if (e==(nx-1)) then
       do j=1 ,ny-1 
         do k=1,nz-1 
           u(nx,j,k)=0 
           u(nx+1,j,k)=-u(nx-1,j,k)
         end do 
       end do
     endif  

     
     do i=s,e
       do k=1,nz-1
         U(i,0,k)=U(i,1,k)  !cavity          
         U(i,-1,k)=U(i,2,k)        
         U(i,ny,k)=U(i,ny-1,k)
         U(i,ny+1,k)=U(i,ny-2,k)
       end do
     end do 

     if (s==1) then
       do i=-1,0
         do k=1,nz-1
           U(i,0,k)=U(i,1,k)       !cavity   
           U(i,-1,k)=U(i,2,k)        
           U(i,ny,k)=U(i,ny-1,k)
           U(i,ny+1,k)=U(i,ny-2,k)
         end do 
       end do 
     end if 
     
     if (e==(nx-1)) then
       do i=nx,nx+1
         do k=1 , nz-1
           U(i,0,k)=U(i,1,k)      !cavity      
           U(i,-1,k)=U(i,2,k)        
           U(i,ny,k)=U(i,ny-1,k)
           U(i,ny+1,k)=U(i,ny-2,k)
         end do 
       end do 
     end if   
   
 
     do i=s,e
       do j=-1,ny+1                   !! Edit April 2015
         U(i,j,-1)=-U(i,j,2)
         U(i,j,0)=-U(i,j,1)
         U(i,j,nz)=-U(i,j,nz-1)
         U(i,j,nz+1)=-U(i,j,nz-2)
       end do 
     end do 

     if (s==1) then
       do i=-1,0
         do j=-1,ny+1                  !! Edit April 2015
           U(i,j,-1)=-U(i,j,2)
           U(i,j,0)=-U(i,j,1)
           U(i,j,nz)=-U(i,j,nz-1)
           U(i,j,nz+1)=-U(i,j,nz-2)
         end do 
       end do 
     end if 

     if (e==(nx-1)) then
       do i=nx,nx+1
         do j=-1,ny+1                   !! Edit April 2015 
           U(i,j,-1)=-U(i,j,2)
           U(i,j,0)=-U(i,j,1)
           U(i,j,nz)=-U(i,j,nz-1)
           U(i,j,nz+1)=-U(i,j,nz-2)
         end do 
       end do
     end if       



      do i=s,e
        do k=1, nz-1
          V(i,1,k)=0
          V(i,0,k)=-V(i,2,k)
          v(i,-1,k)=-v(i,3,k)
          V(i,ny,k)=0
          V(i,ny+1,k)=-V(i,ny-1,k)
        end do 
      end do 


      if (s==1) then
        do j=-1,ny+1 
          do k=1 , nz-1
            V(0,j,k)=V(1,j,k)
            V(-1,j,k)=V(2,j,k)
           end do 
         end do
       end if 
 
       if (e==(nx-1)) then 
         do j=-1,ny+1
           do k=1,nz-1
             V(nx,j,k)=V(nx-1,j,k)
             V(nx+1,j,k)=V(nx-2,j,k)
           end do 
         end do 
       endif 


       do i=s,e
         do j=-1,ny+1       
           V(i,j,-1)=-V(i,j,2)
           V(i,j,0) =-V(i,j,1)
           V(i,j,nz)=-V(i,j,nz-1)
           V(i,j,nz+1)=-V(i,j,nz-2)
         end do
       end do

       if (s==1) then 
         do i=-1,0
           do j=-1,ny+1       
             V(i,j,-1)=-V(i,j,2)
             V(i,j,0) =-V(i,j,1)
             V(i,j,nz)=-V(i,j,nz-1)
             V(i,j,nz+1)=-V(i,j,nz-2)
           end do
         end do
       end if
 
       if (e==(nx-1)) then 
         do i=nx,nx+1
           do j=-1,ny+1       
             V(i,j,-1)=-V(i,j,2)
             V(i,j,0) =-V(i,j,1)
             V(i,j,nz)=-V(i,j,nz-1)
             V(i,j,nz+1)=-V(i,j,nz-2)
           end do
         end do
       end if 


 
       do j=1 ,ny-1 
         do i=s,e
           W(i,j,1)=0
           w(i,j,0)=-w(i,j,2)
           w(i,j,-1)=-w(i,j,3)
           W(i,j,nz)=0
           W(i,j,nz+1)=-W(i,j,nz-1)
         end do 
       end do
 
 
       if (s==1) then           
          do k=-1,nz+1 
            do j=1,ny-1 
              w(0,j,k)=w(1,j,k)
              w(-1,j,k)=w(2,j,k)
            end do 
          end do
        end if 

        if (e==(nx-1)) then  
          do k=-1,nz+1 
            do j=1,ny-1 
              w(nx,j,k)=w(nx-1,j,k)
              w(nx+1,j,k)=w(nx-2,j,k)
            end do 
          end do
        end if 


        do i=s,e 
          do k=-1,nz+1
            w(i,0,k)=w(i,1,k)         !cavity 
            w(i,-1,k)=w(i,2,k)
            w(i,ny,k)=w(i,ny-1,k)
            w(i,ny+1,k)=w(i,ny-2,k)
          end do 
        end do

        if (s==1) then 
          do i=-1,0
            do k=-1,nz+1
              w(i,0,k)=w(i,1,k)     !cavity 
              w(i,-1,k)=w(i,2,k)
              w(i,ny,k)=w(i,ny-1,k)
              w(i,ny+1,k)=w(i,ny-2,k)
            end do 
          end do
        end if 

        if (e==(nx-1)) then
          do i=nx,nx+1 
            do k=-1,nz+1
              w(i,0,k)=w(i,1,k)   !cavity 
              w(i,-1,k)=w(i,2,k)
              w(i,ny,k)=w(i,ny-1,k)
              w(i,ny+1,k)=w(i,ny-2,k)
            end do 
          end do
        endif 
 



          return 
        end subroutine 





Subroutine BoundCond_UVW_Poisson(s,e)

  use M_MPI_General_3D, only: u,v,w,nx,ny,nz

  implicit none 

  integer,intent(in)                                       ::   s,e
  integer i,j,k


   if (s==1) then 
     do j=1 ,ny-1 
       do k=1,nz-1
         u(1,j,k)=0
       end do 
     end do
   end if 
 
   if (e==(nx-1)) then 
     do j=1 ,ny-1 
       do k=1,nz-1 
         u(nx,j,k)=0 
       end do 
     end do
   end if  

   do i=s,e 
     do k=1, nz-1
       V(i,1,k)=0
       V(i,ny,k)=0
     end do 
   end do 

   do j=1 ,ny-1 
     do i=s,e    
       W(i,j,1)=0
       W(i,j,nz)=0
     end do 
   end do 
    

   return 
 end subroutine  



    subroutine Fluid_Dynamic_Ini(p,s,e)
      use M_MPI_General_3D,              only: ny,nz

      implicit none
      integer,intent(in)                                       ::   s,e
      Double precision,dimension (s-1:e+1,0:ny,0:nz)                    ::   p
      
      p(s-1:e+1,0:ny,0:nz)=0.0

      print*, "Fluid Dynamic Initial  condition finished"


      return 
    end subroutine 

    subroutine Advection_Ini(dmux,dpux,dmuy,dpuy,dmuz,dpuz, &
                           & dmvx,dpvx,dmvy,dpvy,dmvz,dpvz,    &
                           & dmwx,dpwx,dmwy,dpwy,dmwz,dpwz,s,e )

      use M_MPI_General_3D, only:u,v,w,hx,hy,hz,nx,ny,nz
      implicit none
   
      integer,intent(in)                                       ::   s,e
      Double precision,dimension (s-1:e+1,0:ny,0:nz),intent(out)        ::   dpux,dmux,dpuy,dmuy,dpuz,dmuz
      Double precision,dimension (s-1:e+1,0:ny,0:nz),intent(out)        ::   dpvx,dmvx,dpvy,dmvy,dpvz,dmvz
      Double precision,dimension (s-1:e+1,0:ny,0:nz),intent(out)        ::   dpwx,dmwx,dpwy,dmwy,dpwz,dmwz
     
      integer i,j,k
     do i=s,e 
       do j=0,ny 
         do k=0,nz

           Dpux(i,j,k)=( u(i+1,j,k)-u(i,j,k)   )/( hx(i)   )
           Dmux(i,j,k)=( u(i,j,k)  -u(i-1,j,k) )/( hx(i-1) )
           Dpuy(i,j,k)=( u(i,j+1,k)-u(i,j,k)   )/( 0.5d0*( hy(j+1)+hy(j)) )
           Dmuy(i,j,k)=( u(i,j,k)  -u(i,j-1,k) )/( 0.5d0*( hy(j)+hy(j-1)) )
           Dpuz(i,j,k)=( u(i,j,k+1)-u(i,j,k)   )/( 0.5d0*( hz(k+1)+hz(k)) )
           Dmuz(i,j,k)=( u(i,j,k)  -u(i,j,k-1) )/( 0.5d0*( hz(k)+hz(k-1)) )

           Dpvx(i,j,k)=( v(i+1,j,k)-v(i,j,k)   )/( 0.5d0*( hx(i+1)+hx(i)) )
           Dmvx(i,j,k)=( v(i,j,k)  -v(i-1,j,k) )/( 0.5d0*( hx(i)+hx(i-1)) )
           Dpvy(i,j,k)=( v(i,j+1,k)-v(i,j,k)   )/( hy(j)   )
           Dmvy(i,j,k)=( v(i,j,k)  -v(i,j-1,k) )/( hy(j-1) )
           Dpvz(i,j,k)=( v(i,j,k+1)-v(i,j,k)   )/( 0.5d0*( hz(k+1)+hz(k)) )
           Dmvz(i,j,k)=( v(i,j,k)  -v(i,j,k-1) )/( 0.5d0*( hz(k)+hz(k-1)) )

           DpWx(i,j,k)=( W(i+1,j,k)-W(i,j,k)   )/( 0.5d0*( hx(i+1)+hx(i)) )
           DmWx(i,j,k)=( W(i,j,k)  -W(i-1,j,k) )/( 0.5d0*( hx(i)+hx(i-1)) )
           DpWy(i,j,k)=( W(i,j+1,k)-W(i,j,k)   )/( 0.5d0*( hy(j+1)+hy(j)) )
           DmWy(i,j,k)=( W(i,j,k)  -W(i,j-1,k) )/( 0.5d0*( hy(j)+hy(j-1)) )
           DpWz(i,j,k)=( W(i,j,k+1)-W(i,j,k)   )/( hz(k)   )
           DmWz(i,j,k)=( W(i,j,k)  -W(i,j,k-1) )/( hz(k-1) )

         end do 
       end do 
     end do 
    
     if (s==1) then
       i=0
       do j=0,ny 
         do k=0,nz
          
           Dpux(i,j,k)=( u(i+1,j,k)-u(i,j,k)   )/( hx(i)   )
           Dmux(i,j,k)=( u(i,j,k)  -u(i-1,j,k) )/( hx(i-1) )
           Dpuy(i,j,k)=( u(i,j+1,k)-u(i,j,k)   )/( 0.5d0*( hy(j+1)+hy(j)) )
           Dmuy(i,j,k)=( u(i,j,k)  -u(i,j-1,k) )/( 0.5d0*( hy(j)+hy(j-1)) )
           Dpuz(i,j,k)=( u(i,j,k+1)-u(i,j,k)   )/( 0.5d0*( hz(k+1)+hz(k)) )
           Dmuz(i,j,k)=( u(i,j,k)  -u(i,j,k-1) )/( 0.5d0*( hz(k)+hz(k-1)) )

           Dpvx(i,j,k)=( v(i+1,j,k)-v(i,j,k)   )/( 0.5d0*( hx(i+1)+hx(i)) )
           Dmvx(i,j,k)=( v(i,j,k)  -v(i-1,j,k) )/( 0.5d0*( hx(i)+hx(i-1)) )
           Dpvy(i,j,k)=( v(i,j+1,k)-v(i,j,k)   )/( hy(j)   )
           Dmvy(i,j,k)=( v(i,j,k)  -v(i,j-1,k) )/( hy(j-1) )
           Dpvz(i,j,k)=( v(i,j,k+1)-v(i,j,k)   )/( 0.5d0*( hz(k+1)+hz(k)) )
           Dmvz(i,j,k)=( v(i,j,k)  -v(i,j,k-1) )/( 0.5d0*( hz(k)+hz(k-1)) )

           DpWx(i,j,k)=( W(i+1,j,k)-W(i,j,k)   )/( 0.5d0*( hx(i+1)+hx(i)) )
           DmWx(i,j,k)=( W(i,j,k)  -W(i-1,j,k) )/( 0.5d0*( hx(i)+hx(i-1)) )
           DpWy(i,j,k)=( W(i,j+1,k)-W(i,j,k)   )/( 0.5d0*( hy(j+1)+hy(j)) )
           DmWy(i,j,k)=( W(i,j,k)  -W(i,j-1,k) )/( 0.5d0*( hy(j)+hy(j-1)) )
           DpWz(i,j,k)=( W(i,j,k+1)-W(i,j,k)   )/( hz(k)   )
           DmWz(i,j,k)=( W(i,j,k)  -W(i,j,k-1) )/( hz(k-1) )

         end do 
       end do 
     end if

     if (e==(nx-1)) then
       i=nx
       do j=0,ny 
         do k=0,nz
          
           Dpux(i,j,k)=( u(i+1,j,k)-u(i,j,k)   )/( hx(i)   )
           Dmux(i,j,k)=( u(i,j,k)  -u(i-1,j,k) )/( hx(i-1) )
           Dpuy(i,j,k)=( u(i,j+1,k)-u(i,j,k)   )/( 0.5d0*( hy(j+1)+hy(j)) )
           Dmuy(i,j,k)=( u(i,j,k)  -u(i,j-1,k) )/( 0.5d0*( hy(j)+hy(j-1)) )
           Dpuz(i,j,k)=( u(i,j,k+1)-u(i,j,k)   )/( 0.5d0*( hz(k+1)+hz(k)) )
           Dmuz(i,j,k)=( u(i,j,k)  -u(i,j,k-1) )/( 0.5d0*( hz(k)+hz(k-1)) )

           Dpvx(i,j,k)=( v(i+1,j,k)-v(i,j,k)   )/( 0.5d0*( hx(i+1)+hx(i)) )
           Dmvx(i,j,k)=( v(i,j,k)  -v(i-1,j,k) )/( 0.5d0*( hx(i)+hx(i-1)) )
           Dpvy(i,j,k)=( v(i,j+1,k)-v(i,j,k)   )/( hy(j)   )
           Dmvy(i,j,k)=( v(i,j,k)  -v(i,j-1,k) )/( hy(j-1) )
           Dpvz(i,j,k)=( v(i,j,k+1)-v(i,j,k)   )/( 0.5d0*( hz(k+1)+hz(k)) )
           Dmvz(i,j,k)=( v(i,j,k)  -v(i,j,k-1) )/( 0.5d0*( hz(k)+hz(k-1)) )

           DpWx(i,j,k)=( W(i+1,j,k)-W(i,j,k)   )/( 0.5d0*( hx(i+1)+hx(i)) )
           DmWx(i,j,k)=( W(i,j,k)  -W(i-1,j,k) )/( 0.5d0*( hx(i)+hx(i-1)) )
           DpWy(i,j,k)=( W(i,j+1,k)-W(i,j,k)   )/( 0.5d0*( hy(j+1)+hy(j)) )
           DmWy(i,j,k)=( W(i,j,k)  -W(i,j-1,k) )/( 0.5d0*( hy(j)+hy(j-1)) )
           DpWz(i,j,k)=( W(i,j,k+1)-W(i,j,k)   )/( hz(k)   )
           DmWz(i,j,k)=( W(i,j,k)  -W(i,j,k-1) )/( hz(k-1) )

         end do 
       end do 
     end if





     return
   end subroutine  


   subroutine Advection(dmux,dpux,dmuy,dpuy,dmuz,dpuz, &
                      & dmvx,dpvx,dmvy,dpvy,dmvz,dpvz, &
                      & dmwx,dpwx,dmwy,dpwy,dmwz,dpwz, &
                      & advectu,advectv,advectw,s,e    )
      
     use M_MPI_General_3D, only: nx,ny,nz,x,y,z,u,v,w
     implicit none 
 
     integer,intent(in)                                       ::   s,e
     Double precision,dimension (s-1:e+1,0:ny,0:nz) ,intent(in)        ::   dpux,dmux,dpuy,dmuy,dpuz,dmuz
     Double precision,dimension (s-1:e+1,0:ny,0:nz) ,intent(in)        ::   dpvx,dmvx,dpvy,dmvy,dpvz,dmvz
     Double precision,dimension (s-1:e+1,0:ny,0:nz) ,intent(in)        ::   dpwx,dmwx,dpwy,dmwy,dpwz,dmwz
     Double precision,dimension (s:e,1:ny-1,1:nz-1) ,intent(out)       ::   advectu,advectv,advectw

     Double precision ux,uy,uz,vx,vy,vz,wx,wy,wz
     integer i,j,k


     do i=s,e 
       do j=1,ny-1 
         do k=1,nz-1     !!A!!


 
          if (u(i,j,k).gt.0.0) then

            if (  abs(  Dmux(i,j,k)-Dmux(i-1,j,k) ).lt.abs(  Dpux(i,j,k)-Dpux(i-1,j,k) )   ) then
              ux=Dmux(i,j,k)+  0.5*(  Dmux(i,j,k)-Dmux(i-1,j,k)  ) 
            else
              ux=Dmux(i,j,k)+  0.5*(  Dpux(i,j,k)-Dpux(i-1,j,k)  )
            end if 
    
          else

            if (  abs(  Dmux(i+1,j,k)-Dmux(i,j,k) ).lt.abs(  Dpux(i+1,j,k)-Dpux(i,j,k) )   ) then
              ux=Dpux(i,j,k)-  0.5*(  Dmux(i+1,j,k)-Dmux(i,j,k)  ) 
            else
              ux=Dpux(i,j,k)-  0.5*(  Dpux(i+1,j,k)-Dpux(i,j,k)  )
            end if 

          end if 

!!U2
          if (  (0.25*( V(i,j,k)+V(i,j+1,k)+v(i-1,j,k)+v(i-1,j+1,k) )).gt.0.0) then

            if (  abs(  DmuY(i,j,k)-DmuY(i,j-1,k) ).lt.abs(  DpuY(i,j,k)-DpuY(i,j-1,k) )   ) then
              uY=DmuY(i,j,k)+  0.5*(  DmuY(i,j,k)-DmuY(i,j-1,k)  ) 
            else
              uY=DmuY(i,j,k)+  0.5*(  DpuY(i,j,k)-DpuY(i,j-1,k)  )
            end if 
    
          else

            if (  abs(  DmuY(i,j+1,k)-DmuY(i,j,k) ).lt.abs(  DpuY(i,j+1,k)-DpuY(i,j,k) )   ) then
              uY=DpuY(i,j,k)-  0.5*(  DmuY(i,j+1,k)-DmuY(i,j,k)  ) 
            else
              uY=DpuY(i,j,k)-  0.5*(  DpuY(i,j+1,k)-DpuY(i,j,k)  )
            end if 

          end if 
 
 
          if (  (0.25*( W(i,j,k)+W(i-1,j,k)+W(i,j,k+1)+W(i-1,j,k+1) )).gt.0.0) then
  
            if (  abs(  DmuZ(i,j,k)-DmuZ(i,j,k-1) ).lt.abs(  DpuZ(i,j,k)-DpuZ(i,j,k-1) )   ) then
              uZ=DmuZ(i,j,k)+  0.5*(  DmuZ(i,j,k)-DmuZ(i,j,k-1)  ) 
            else
              uZ=DmuZ(i,j,k)+  0.5*(  DpuZ(i,j,k)-DpuZ(i,j,k-1)  )
            end if
    
          else 
    
            if (  abs(  Dmuz(i,j,k+1)-Dmuz(i,j,k) ).lt.abs(  DpuZ(i,j,k+1)-DpuZ(i,j,k) )   ) then
              uZ=DpuZ(i,j,k)-  0.5*(  DmuZ(i,j,k+1)-DmuZ(i,j,k)  ) 
            else
              uZ=DpuZ(i,j,k)-  0.5*(  DpuZ(i,j,k+1)-DpuZ(i,j,k)  )
            end if 
    
          end if 
    

          advectu(i,j,k)=-(                                      u(i,j,k)*uX+ &
          &          0.25*( V(i,j,k)+V(i,j+1,k)+v(i-1,j,k)+v(i-1,j+1,k) )*uY+ &
          &          0.25*( W(i,j,k)+W(i-1,j,k)+W(i,j,k+1)+W(i-1,j,k+1) )*uZ  )



        end do 
      end do 
    end do 
  


    do i=s,e 
      do j=1,ny-1
        do k=1,nz-1 !!A!!
 
 
          if (  (0.25*( u(i,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i+1,j-1,k) )).gt.0.0) then

            if (  abs(  Dmvx(i,j,k)-Dmvx(i-1,j,k) ).lt.abs(  Dpvx(i,j,k)-Dpvx(i-1,j,k) )   ) then
              vx=Dmvx(i,j,k)+  0.5*(  Dmvx(i,j,k)-Dmvx(i-1,j,k)  ) 
            else
              vx=Dmvx(i,j,k)+  0.5*(  Dpvx(i,j,k)-Dpvx(i-1,j,k)  )
            end if 
    
          else

            if (  abs(  Dmvx(i+1,j,k)-Dmvx(i,j,k) ).lt.abs(  Dpvx(i+1,j,k)-Dpvx(i,j,k) )   ) then
              vx=Dpvx(i,j,k)-  0.5*(  Dmvx(i+1,j,k)-Dmvx(i,j,k)  ) 
            else
              vx=Dpvx(i,j,k)-  0.5*(  Dpvx(i+1,j,k)-Dpvx(i,j,k)  )
            end if 

          end if 


          if ( V(i,j,k).gt.0.0) then

            if (  abs(  DmvY(i,j,k)-DmvY(i,j-1,k) ).lt.abs(  DpvY(i,j,k)-DpvY(i,j-1,k) )   ) then
              vY=DmvY(i,j,k)+  0.5*(  DmvY(i,j,k)-DmvY(i,j-1,k)  ) 
            else
              vY=DmvY(i,j,k)+  0.5*(  DpvY(i,j,k)-DpvY(i,j-1,k)  )
            end if 
    
          else

            if (  abs(  DmvY(i,j+1,k)-DmvY(i,j,k) ).lt.abs(  DpvY(i,j+1,k)-DpvY(i,j,k) )   ) then
              vY=DpvY(i,j,k)-  0.5*(  DmvY(i,j+1,k)-DmvY(i,j,k)  ) 
            else
              vY=DpvY(i,j,k)-  0.5*(  DpvY(i,j+1,k)-DpvY(i,j,k)  )
            end if 

          end if 
 
 !!V3
          if (  (0.25*( W(i,j,k)+W(i,j-1,k)+W(i,j,k+1)+W(i,j-1,k+1) )).gt.0.0) then

            if (  abs(  DmvZ(i,j,k)-DmvZ(i,j,k-1) ).lt.abs(  DpvZ(i,j,k)-DpvZ(i,j,k-1) )   ) then
              vZ=DmvZ(i,j,k)+  0.5*(  DmvZ(i,j,k)-DmvZ(i,j,k-1)  ) 
            else
              vZ=DmvZ(i,j,k)+  0.5*(  DpvZ(i,j,k)-DpvZ(i,j,k-1)  )
            end if 
    
          else

            if (  abs(  DmvZ(i,j,k+1)-DmvZ(i,j,k) ).lt.abs(  DpvZ(i,j,k+1)-DpvZ(i,j,k) )   ) then
              vZ=DpvZ(i,j,k)-  0.5*(  DmvZ(i,j,k+1)-DmvZ(i,j,k)  ) 
            else
              vZ=DpvZ(i,j,k)-  0.5*(  DpvZ(i,j,k+1)-DpvZ(i,j,k)  )
            end if 

          end if 


          advectv(i,j,k)=-(0.25*( u(i,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i+1,j-1,k) )*vX+   &
          &                                                            V(i,j,k)*vY+   &
          &                0.25*( W(i,j,k)+W(i,j,k+1)+W(i,j-1,k)+W(i,j-1,k+1) )*vZ    )



        end do 
      end do 
    end do


  

    do i=s,e 
      do j=1,ny-1  
        do k=1,nz-1 !!A!!
 
 
          if (  (0.25*( u(i,j,k)+u(i+1,j,k)+u(i,j,k-1)+u(i+1,j,k-1) )).gt.0.0) then

            if (  abs(  DmWx(i,j,k)-DmWx(i-1,j,k) ).lt.abs(  DpWx(i,j,k)-DpWx(i-1,j,k) )   ) then
              Wx=DmWx(i,j,k)+  0.5*(  DmWx(i,j,k)-DmWx(i-1,j,k)  ) 
            else
              Wx=DmWx(i,j,k)+  0.5*(  DpWx(i,j,k)-DpWx(i-1,j,k)  )
            end if 
    
          else

            if (  abs(  DmWx(i+1,j,k)-DmWx(i,j,k) ).lt.abs(  DpWx(i+1,j,k)-DpWx(i,j,k) )   ) then
              Wx=DpWx(i,j,k)-  0.5*(  DmWx(i+1,j,k)-DmWx(i,j,k)  ) 
            else
              Wx=DpWx(i,j,k)-  0.5*(  DpWx(i+1,j,k)-DpWx(i,j,k)  )
            end if 

          end if 


          if (  (0.25*( V(i,j,k)+V(i,j+1,k)+v(i,j,k-1)+v(i,j+1,k-1) )).gt.0.0) then

            if (  abs(  DmWY(i,j,k)-DmWY(i,j-1,k) ).lt.abs(  DpWY(i,j,k)-DpWY(i,j-1,k) )   ) then
              WY=DmWY(i,j,k)+  0.5*(  DmWY(i,j,k)-DmWY(i,j-1,k)  ) 
            else
              WY=DmWY(i,j,k)+  0.5*(  DpWY(i,j,k)-DpWY(i,j-1,k)  )
            end if 
    
          else

            if (  abs(  DmWY(i,j+1,k)-DmWY(i,j,k) ).lt.abs(  DpWY(i,j+1,k)-DpWY(i,j,k) )   ) then
              WY=DpWY(i,j,k)-  0.5*(  DmWY(i,j+1,k)-DmWY(i,j,k)  ) 
            else
              WY=DpWY(i,j,k)-  0.5*(  DpWY(i,j+1,k)-DpWY(i,j,k)  )
            end if 

          end if 

 
          if ( W(i,j,k).gt.0.0) then

            if (  abs(  DmWZ(i,j,k)-DmWZ(i,j,k-1) ).lt.abs(  DpWZ(i,j,k)-DpWZ(i,j,k-1) )   ) then
              WZ=DmWZ(i,j,k)+  0.5*(  DmWZ(i,j,k)-DmWZ(i,j,k-1)  ) 
            else
              WZ=DmWZ(i,j,k)+  0.5*(  DpWZ(i,j,k)-DpWZ(i,j,k-1)  )
            end if 
    
          else

            if (  abs(  DmWZ(i,j,k+1)-DmWZ(i,j,k) ).lt.abs(  DpWZ(i,j,k+1)-DpWZ(i,j,k) )   ) then
              WZ=DpWZ(i,j,k)-  0.5*(  DmWZ(i,j,k+1)-DmWZ(i,j,k)  ) 
            else
              WZ=DpWZ(i,j,k)-  0.5*(  DpWZ(i,j,k+1)-DpWZ(i,j,k)  )
            end if 

           end if 

           advectW(i,j,k)=-(0.25*( u(i,j,k)+u(i+1,j,k)+u(i,j,k-1)+u(i+1,j,k-1) )*WX+ &
           &                0.25*( V(i,j,k)+V(i,j+1,k)+v(i,j,k-1)+v(i,j+1,k-1) )*WY  &
           &                                                           +W(i,j,k)*WZ  )
 


         end do 
       end do 
     end do 


     return 

   end subroutine 


   subroutine Viscoset(miuv,ro,Tx,Ty,Tz,s,e)
     use M_MPI_General_3D, only: ny,nz,x,y,z,u,v,w,hx,hy,hz

     implicit none 
     integer,intent(in)                                      ::   s,e
     Double precision,dimension (s-1:e+1,0:ny,0:nz),intent(in)        ::   ro,miuv
     Double precision,dimension(s:e,1:ny-1,1:nz-1),intent(out)        ::   Tx,Ty,Tz

     Double precision Txxr,Txxl,Tyxd,Tyxu,Tzxf,Tzxb
     Double precision Tyyu,Tyyd,Txyr,Txyl,Tzyf,Tzyb
     Double precision Tzzf,Tzzb,Tyzd,Tyzu,Txzr,Txzl
     integer i,j,k


     Do i=s,e 
       Do j=1, ny-1 
         do k=1, nz-1


           Txxr=2*miuv(i,j,k)  *( U(i+1,j,K)-U(i,j,K)   )/hx(i)    !( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )
           Txxl=2*miuv(i-1,j,k)*( U(i,j,K)  -U(i-1,j,K) )/hx(i-1)  !( 0.5*(x(i)  +x(i-1))-0.5*(x(i-1)+x(i-2)) )

           Tyxu=0.25*( miuv(i,j,k)+miuv(i,j+1,k)+miuv(i-1,j,k)+miuv(i-1,j+1,k) )*&
           &(      ( U(i,j+1,K)-U(i,j,K)     )/( Y(j+1)-y(j) )+ ( V(i,j+1,K)-V(i-1,j+1,K) )/( x(i)-x(i-1) )     )
                 
           Tyxd=0.25*( miuv(i,j-1,k)+miuv(i,j,k)+miuv(i-1,j-1,k)+miuv(i-1,j,k) )*&
           &(      ( U(i,j,K)-U(i,j-1,K)     )/( y(j)-y(j-1) )+ ( V(i,j,K)  -V(i-1,j,K)   )/( x(i)-x(i-1) )     )

           Tzxf=0.25*( miuv(i,j,k)+miuv(i,j,k+1)+miuv(i-1,j,k)+miuv(i-1,j,k+1) )*&
           &(      ( U(i,j,K+1)-U(i,j,K)     )/( z(k+1)-z(k) )+ ( W(i,j,K+1)-W(i-1,j,K+1) )/( x(i)-x(i-1) )     )
  
           Tzxb=0.25*( miuv(i,j,k-1)+miuv(i,j,k)+miuv(i-1,j,k-1)+miuv(i-1,j,k) )*&
           &(      ( U(i,j,K)-U(i,j,K-1)     )/( z(k)-z(k-1) )+ ( W(i,j,K)  -W(i-1,j,K)   )/( x(i)-x(i-1) )     )

           Tx(i,j,k)= (  (Txxr-Txxl)/(x(i)-x(i-1))  + &
           &             (Tyxu-Tyxd)/hy(j)          + &
           &             (Tzxf-Tzxb)/hz(k)        )/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )
 
   
         end do 
       end do 
     end do  


     Do j=1,ny-1
       Do i=s,e
         do k=1, nz-1


           Tyyu=2*miuv(i,j,k)   *( V(i,j+1,k)-V(i,j,k)   )/hy(j)
           Tyyd=2*miuv(i,j-1,k)*( V(i,j,k)  -V(i,j-1,k) )/hy(j-1)
 
           Txyr=0.25*( miuv(i,j,k)+miuv(i,j-1,k)+miuv(i+1,j,k)+miuv(i+1,j-1,k) )*&
           & (     ( V(i+1,j,k)-V(i,j,k)     )/( x(i+1)-x(i) )+ ( U(i+1,j,k)  -U(i+1,j-1,k) )/( y(j)-y(j-1)  )   )
  
           Txyl=0.25*( miuv(i-1,j,k)+miuv(i-1,j-1,k)+miuv(i,j,k)+miuv(i,j-1,k) )*&
           & (     ( V(i,j,k)-V(i-1,j,k)     )/( x(i)-x(i-1) )+ ( U(i,j,k)    -U(i,j-1,k)   )/( y(j)-y(j-1)  )   )
  
           Tzyf=0.25*( miuv(i,j,k)+miuv(i,j-1,k)+miuv(i,j,k+1)+miuv(i,j-1,k+1) )*&
           & (     ( v(i,j,k+1)-V(i,j,k)     )/( z(k+1)-z(k) )+ (W(i,j,k+1)  -W(i,j-1,k+1)  )/( y(j)-y(j-1)  )   )
 
           Tzyb=0.25*( miuv(i,j,k-1)+miuv(i,j-1,k-1)+miuv(i,j,k)+miuv(i,j-1,k) )*&                         
           & (     ( v(i,j,k)-V(i,j,k-1)     )/( z(k)-z(k-1) )+ (W(i,j,k)    -W(i,j-1,k)    )/( y(j)-y(j-1)  )   )  
  

           Ty(i,j,k)= (  (Txyr-Txyl)/hx(i)         + &
                      &  (Tyyu-Tyyd)/(y(j)-y(j-1)) + &
                      &  (Tzyf-Tzyb)/hz(k)         )/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )
 
   
         end do 
       end do 
     end do 


     Do j=1,ny-1
       Do i=s,e 
         do k=1, nz-1


           Tzzf=2*miuv(i,j,k)  * ( w(i,j,k+1)-w(i,j,k)   )/hz(k)
           Tzzb=2*miuv(i,j,k-1)* ( w(i,j,k)  -w(i,j,k-1) )/hz(k-1)

           Txzr=0.25*( miuv(i,j,k)+miuv(i+1,j,k)+miuv(i,j,k-1)+miuv(i+1,j,k-1) )*&
           &(      ( U(i+1,j,K)-U(i+1,j,K-1)     )/( z(k)-z(k-1) )+ ( W(i+1,j,K)-W(i,j,K) )/( x(i+1)-x(i) )     )
  
           Txzl=0.25*( miuv(i-1,j,k)+miuv(i,j,k)+miuv(i-1,j,k-1)+miuv(i,j,k-1) )*&
           &(      ( U(i,j,K)  -U(i,j,K-1)       )/( z(k)-z(k-1) )+ ( W(i,j,K)-W(i-1,j,K) )/( x(i)-x(i-1) )     )
   
           Tyzu=0.25*( miuv(i,j,k)+miuv(i,j+1,k)+miuv(i,j,k-1)+miuv(i,j+1,k-1) )*&
           &(      ( V(i,j+1,K) -V(i,j+1,K-1)    )/( z(k)-z(k-1) )+ ( W(i,j+1,K)-W(i,j,K) )/( y(j+1)-y(j) )     )

           Tyzd=0.25*( miuv(i,j-1,k)+miuv(i,j,k)+miuv(i,j-1,k-1)+miuv(i,j,k-1) )*&
           &(      ( V(i,j,K)   -V(i,j,K-1)      )/( z(k)-z(k-1) )+ ( W(i,j,K)-W(i,j-1,K) )/( y(j)-y(j-1) )     )
   
           Tz(i,j,k)= (    (Txzr-Txzl)/hx(i)+ &
                      &    (Tyzu-Tyzd)/hy(j)+ &
                      &    (Tzzf-Tzzb)/(z(k)-z(k-1))       )/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )  
 

         end do 
       end do 
     end do 


     return 
   end subroutine 


  subroutine Proprety(pro,prodrop,proAir,s,e)
    use M_MPI_General_3D,             only: ny,nz
    use M_MPI_FrSrf_3D,               only: phi,eps
    implicit none 

    integer,intent(in)                                               ::  s,e
    dOUBLE PRECISION,intent(in)                                      ::  prodrop,proAir
    double precision,dimension (s-1:e+1,0:ny,0:nz),intent(out)       ::  pro
    
    double precision                                                 ::  HV
    integer i,j,k 

    do i=s-1,e+1
      do j=0,ny
        do k=0 ,nz

          pro(i,j,k)=proAir+ HV(-phi(i,j,k),eps)*(prodrop-proAir)
          ! pro(i,j,k)=prodrop
        end do
      end do
    end do


    return 

  end subroutine 
 
double precision function HV(phi_D,eps)

double precision            ::  phi_D,pi,eps


pi=3.141592654D0
if (phi_D.gt.eps) then
HV=1.d0
else if (phi_D.lt.-eps) then
HV=0.d0
else
HV=0.5d0*(  1.d0+ phi_D/eps +(1/pi)*dsin(pi*phi_D/eps)  )
end if

return
end function

  
subroutine Poisson_COF(ro,dt,Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,s,e)
     use M_MPI_General_3D, only: u,v,w,x,y,z,nx,ny,nz
     implicit none
      
     integer                              ,intent(in)     ::   s,e
     Double precision,dimension (s-1:e+1,0:ny,0:nz),intent(in)     ::   ro
     Double precision                              ,intent(in)     ::   dt
     Double precision,dimension (s:e,1:ny-1,1:nz-1),intent(out)    ::   Apx,Amx,Apy,Amy,Apz,Amz,AP,Q

     integer                                              ::   i,j,k


     do i=s,e 
       do j=1,ny-1
         do k=1,nz-1



           if (i==nx-1) then
              Apx(i,j,k)=0
           else
              Apx(i,j,k)=(1)/( (x(i+1)-x(i))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i+1,j,k))  )
            end if 
 
            if (i==1)then 
              Amx(i,j,k)=0       
            else 
              Amx(i,j,k)=(1)/( (x(i)-x(i-1))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )
             end if 
 
             if (j==ny-1) then
               Apy(i,j,k)=0
             else
               Apy(i,j,k)=(1)/( (y(j+1)-y(j))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j+1,k))  )
             end if 
 
             if (j==1) then
               Amy(i,j,k)=0
             else
               Amy(i,j,k)=(1)/( (y(j)-y(j-1))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )
             end if 

             if (k==nz-1) then
               Apz(i,j,k)=0
             else
               Apz(i,j,k)=(1)/( (z(k+1)-z(k))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k+1))  )  
             end if 
 
             if (k==1) then
               Amz(i,j,k)=0   
             else
               Amz(i,j,k)=(1)/( (z(k)-z(k-1))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )
             end if 


             AP(i,j,k)=-(   Apx(i,j,k)+ &
                          & Amx(i,j,k)+ &
                          & Apy(i,j,k)+ &
                          & Amy(i,j,k)+ &
                          & Apz(i,j,k)+ &
                          & Amz(i,j,k) ) 
  
             Q(I,J,K)=(   (U(I+1,J,k)-U(I,J,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
                       &  (V(I,J+1,k)-V(I,J,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+& 
                       &  (W(I,j,k+1)-W(I,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )   )/dt 



           end do 
         end do 
       end do

       return

     end subroutine 



subroutine Poisson (Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,beta,pdif,p,s,e)
  use M_MPI_General_3D,           only: nx,ny,nz
  implicit none 
 
  integer                              ,intent(in)             ::   s,e
  Double precision,dimension (s:e,1:ny-1,1:nz-1),intent(in)             ::   Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q
  Double precision,dimension (s-1:e+1,0:ny,0:nz),intent(inout)          ::   p
  Double precision                              ,intent(in)             ::   pdif,beta

  Double precision,dimension (s-1:e+1,0:ny,0:nz)                        ::   pold
  Double precision                                                      ::   maxp
  integer                                                      ::   i,j,k,number2


  number2=0
  maxp=1000.0
  do while (maxp.gt.pdif.AND.number2.lt.2500) 


    maxp=0 
    number2=number2+1
 
    do i=s,e 
      do j=1,ny-1 
        do k=1, nz-1

          pold(i,j,k)=p(i,j,k)
 
          p(i,j,k)=beta*( ( Apx(i,j,k)*P(I+1,j,k)+Amx(i,j,k)*P(I-1,j,k)+Apy(i,j,k)*P(I,j+1,k) &
          &                +Amy(i,j,k)*P(I,j-1,k)+Apz(i,j,k)*P(I,j,k+1)+Amz(i,j,k)*P(I,j,k-1)-Q(I,j,k)  )/(-AP(i,j,k))  ) &
          & +(1-beta)*p(i,j,k)

          if (abs( pold(i,j,k)-p(i,j,k) ).gt.maxp) then 
            maxp=abs( pold(i,j,k)-p(i,j,k) )
          end if
 
        end do 
      end do 
    end do 


  end do  !!while !
  print *,maxp,number2 


  if (s==1) then  
    do j=1,ny-1
      do k=1,nz-1 
        p(0,j,k)=p(1,j,k)
      end do 
    end do 
  end if 
  
  if (e==(nx-1)) then 
    do j=1,ny-1
      do k=1,nz-1 
        p(nx,j,k)=p(nx-1,j,k)
      end do
    end do 
  end if 


  do i=s,e
    do k=1,nz-1 
      p(i,0,k)=p(i,1,k)
      p(i,ny,k)=p(i,ny-1,k)
    end do
  end do 


  do i=s,e
    do j=1,ny-1
      p(i,j,0)=p(i,j,1)
      p(i,j,nz)=p(i,j,nz-1)
    end do 
  end do 


  return 

end subroutine  


   




subroutine PrintData(p,ro,plot,tp,TStepE,dt,MyRank,numprocs,s,e,comm1d)

!Parallel
  use mpi
  use M_MPI_General_3D, only: u,v,w,x,y,z,nx,ny,nz,Parallel
  use M_MPI_FrSrf_3D,   only: phi
!Parallel
  implicit none

  integer,intent(in)                                               :: s,e
  double precision,dimension (s-1:e+1,0:ny,0:nz),intent(in)        :: p,ro
  dOUBLE PRECISION,intent(in)                                      :: dt
  integer,intent(in)                                               :: MyRank,numprocs,comm1d
  integer,intent(in)                                               :: tp,plot,TStepE

  double precision,dimension(:,:,:),allocatable                    :: us,vs,ws,phis
  !double precision,dimension ((e-s+1)*(ny-1)*(nz-1))                      :: Upak,VPak,WPak,PPak,roPak
  double precision,dimension (:),allocatable                       :: Upak2,VPak2,WPak2,PPak2,roPak2,phiPak2
  double precision,dimension (:),allocatable                       :: UTo,VTo,WTo,PTo,roTo,phiTo
  double precision,dimension (1:nx-1)                              :: XTo
  double precision                                                 :: ZFreeG,ZFreeG_cpu
  
  integer                                                          :: status(MPI_STATUS_SIZE), ierr
  integer                                                          :: i,j,k,kk,kkk
  integer,save                                                     :: count=100000

    count=count+1
    if (count.ge.plot)then

      allocate ( us(s-1:e+1,0:ny,0:nz),vs(s-1:e+1,0:ny,0:nz),ws(s-1:e+1,0:ny,0:nz),phis(s-1:e+1,0:ny,0:nz) )
      allocate ( Upak2((nx-1)*(ny-1)*(nz-1)), VPak2((nx-1)*(ny-1)*(nz-1)),  WPak2((nx-1)*(ny-1)*(nz-1)), &
           &     PPak2((nx-1)*(ny-1)*(nz-1)),roPak2((nx-1)*(ny-1)*(nz-1)),phiPak2((nx-1)*(ny-1)*(nz-1)) )

      allocate ( UTo((nx-1)*(ny-1)*(nz-1)), VTo((nx-1)*(ny-1)*(nz-1)),  WTo((nx-1)*(ny-1)*(nz-1)), &
           &     PTo((nx-1)*(ny-1)*(nz-1)),roTo((nx-1)*(ny-1)*(nz-1)),phiTo((nx-1)*(ny-1)*(nz-1)) )
      count=0
      print*,"data is written"
      call flush(6)
  
      us  (s-1:e+1,0:ny,0:nz)=u(s-1:e+1,0:ny,0:nz)
      vs  (s-1:e+1,0:ny,0:nz)=v(s-1:e+1,0:ny,0:nz)
      ws  (s-1:e+1,0:ny,0:nz)=w(s-1:e+1,0:ny,0:nz)
      phis(s-1:e+1,0:ny,0:nz)=phi(s-1:e+1,0:ny,0:nz)

      Upak2(:)=0 
      VPak2(:)=0
      WPak2(:)=0 
      PPak2(:)=0 
      roPak2(:)=0
      phiPak2(:)=0 
      do k=1,nz-1
        kk=(k-1)*(e-s+1)*(ny-1)
        do j=1,ny-1
          do i=s,e
            kk=kk+1
            kkk=kk+ (myRank)*(e-s+1) + floor( real(kk-1)/(e-s+1) )*(numprocs-1)*(e-s+1)
            roPak2(kkk)=ro(i,j,k)
            PPak2(kkk) =p(i,j,k)
            UPak2(kkk) =us(i,j,k)
            VPak2(kkk) =vs(i,j,k)
            WPak2(kkk) =ws(i,j,k)
            phiPak2(kkk) =phis(i,j,k)
          end do
        end do
      end do

      !if (s==1) then
      !  print*,p(10,20,38),"point 10"
      !end if 

      !if (e==(nx-1)) then
      !  print*,p(30,20,38),"point 30"
      !end if 

      if (Parallel) then
        
        call MPI_GATHER (X(s)       ,e-s+1                ,MPI_dOUBLE_PRECISION ,XTo(1)   ,e-s+1          ,MPI_dOUBLE_PRECISION,0,comm1d,ierr )  !! I don't know if the receive buffer should be one or not
        call MPI_Allreduce(upak2(1)   ,uTo(1)   ,(nx-1)*(ny-1)*(nz-1),MPI_dOUBLE_PRECISION,MPI_sum,comm1d,ierr)
        call MPI_Allreduce(vpak2(1)   ,vTo(1)   ,(nx-1)*(ny-1)*(nz-1),MPI_dOUBLE_PRECISION,MPI_sum,comm1d,ierr)
        call MPI_Allreduce(wpak2(1)   ,wTo(1)   ,(nx-1)*(ny-1)*(nz-1),MPI_dOUBLE_PRECISION,MPI_sum,comm1d,ierr)
        call MPI_Allreduce(phipak2(1) ,phiTo(1) ,(nx-1)*(ny-1)*(nz-1),MPI_dOUBLE_PRECISION,MPI_sum,comm1d,ierr)
        call MPI_Allreduce(ropak2(1) ,roTo(1)   ,(nx-1)*(ny-1)*(nz-1),MPI_dOUBLE_PRECISION,MPI_sum,comm1d,ierr)
        call MPI_Allreduce(ppak2(1)   ,PTo(1)   ,(nx-1)*(ny-1)*(nz-1),MPI_dOUBLE_PRECISION,MPI_sum,comm1d,ierr)

      end if

      if (MyRank.eq.0) then

        kk=0
        write(35,*) 'ZONE T="',tp,'" SOLUTIONTIME=',tp*dt,'i=',nx-1,' j=',ny-1,' k=',nz-1 
        call flush(35)
        do k=1,nz-1
          Do j=1,ny-1
            Do i=1,nx-1
              kk=kk+1
              
              write(35,135) xTo(i),Y(j),Z(k),UTo(kk),VTo(kk),WTo(kk),pTo(kk),roTo(kk),phiTo(kk)
               !write(35,135) xTo(i),Y(j),UTo(kk),VTo(kk),WTo(kk),pTo(kk),roTo(kk)
          !    write(35,135) xTo(i),Y(j),Z(k),UTo(kk),roTo(kk)
              call flush(35)
            end do
          end do
        end do
 
      end if



   end if
   135  format (9(1x,e15.7))
  !  135  format (7(1x,e15.7))


  ZFreeG_cpu=0
  if ( s.le.int((nx-1)/2).AND.e.ge.int((nx-1)/2) ) then

    print*,"MyRank for free surface is",MyRank 
    do k=1,nz-1
      if ( (Phi(int((nx-1)/2),int((ny-1)/2),k)*Phi(int((nx-1)/2),int((ny-1)/2),k+1)).le.0 ) then
        ZFreeG_cpu=Z(k)
        exit
      end if
    end do
    print*, ZFreeG_cpu
    
  endif
 
  if (Parallel) then
    call MPI_AllReduce( ZFreeG_cpu,ZFreeG,1,MPI_DOUBLE_PRECISION,MPI_sum,comm1d, ierr )
  end if


  if (MyRank.eq.0) then
    write (65,165) tp*dt,',',ZFreeG
  end if
  165  format (2(1x,e15.7,1A))



    if (tp==tstepE) then

      print *,"I am in last time step cavity write"
      if (myRank==0) then 

        Do k=1,nz-1
          write(45,145) UTo( (nx-1)*(ny-1)*(k-1)+ (nx-1)*int((ny-1)/2)+int((nx-1)/2) )
          call flush(45)
        end do

        Do i=1,nx-1
          write(55,155) WTo( (nx-1)*(ny-1)*(int((nz-1)/2))+ (nx-1)*(int((ny-1)/2)-1)+i)
          call flush(55)
        end do

      end if 

    end if
    145  format (1(1x,e15.7))
    155  format (1(1x,e15.7))



return

end subroutine 









     







