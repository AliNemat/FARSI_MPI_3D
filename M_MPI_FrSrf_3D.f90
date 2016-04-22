module M_MPI_FrSrf_3D

  use mpi

  implicit none
  double precision,dimension (:,:,:),allocatable       ::   phi
  double precision                                     ::   eps
  double precision                                     ::   dtau
  logical, parameter                                   ::   TestFrsr=.false.


  contains


  subroutine Levelset_Ini(MyRank,comm1d,s,e)

    use M_MPI_General_3D, only: x,y,z,ZFree,nx,ny,nz,parallel,Lx,pi
    implicit none

    integer                                     ,intent(in)          ::  s,e
    integer                                     ,intent(in)          ::  MyRank
    integer,intent(in)                                               ::  comm1d

    integer                                                          ::  ierr
    double precision                                                 ::  avgh_cpu,avgh,bb
    double precision                                                 ::  phiges,cellX,Xipr
    integer                                                          ::  i,j,k,ipr,jpr

    allocate( Phi(s-2:e+2,-1:ny+1,-1:nz+1) )
    do i=s-2,e+2
      do j=-1,ny+1
        do k=-1,nz+1
          phi(i,j,k)=z(k)-zfree
        end do
     end do 
    end do

   eps=3*( z(int(nz/2))-z(int(nz/2)-1) )
   print*,eps, "eps equal to from Cpu",MyRank


    avgh_cpu=1000 !! large number!! 
    do i=s,e 
      do j=1, ny-1 
        do k=1,nz-1

          bb=min ( (x(i+1)-x(i)),(y(j+1)-y(j)),(z(k+1)-z(k)) )
          if (bb.lt.avgh_cpu) then
            avgh_cpu=bb
          end if

        end do 
      end do
    end do

    if (Parallel) then
      call MPI_Allreduce( avgh_cpu  ,avgh  ,1 ,MPI_Double_Precision, MPI_Min, comm1d, ierr )
    endif

    dtau=0.5*avgh
    !! for validation purpose 
    !!!!!!!!!!!!!!!!!! for free surface testing with Wu et al. (2001)
    !only for uniform grid

    !CellX=Lx/(nx-1)
    if (TestFrsr) then
      do k=-1,nz+1 
        do j=-1,ny+1 
          do i=s-2,e+2

            phi(i,j,k)=10000
          !do ipr=1,nx-1
          !  Xipr=CellX*(ipr-0.5)
          !  phiges=sqrt(   ( x(i)-xipr )**2  +  (  z(k)-( ZFree+0.05*cos(pi*xipr) )  )**2   )

            do jpr=1,ny-1
              phiges=sqrt(   ( Y(j)-Y(jpr) )**2  +  (  z(k)-( ZFree+0.05*cos(pi*Y(jpr)) )  )**2   )
              if (phiges.lt.phi(i,j,k) ) then 
                phi(i,j,k)=phiges
              end if 
            end do 
      
          end do 
        end do 
      end do 
 
      do k=-1,nz+1 
        do j=-1,ny+1 
          do i=s-2,e+2

            if ( z(k).lt.( ZFree+0.05*cos(pi*(y(j)) ) ) ) then 
              phi(i,j,k)=-phi(i,j,k)
            end if 

          end do 
        end do
      end do  
    end if 


    return

  end subroutine


  subroutine levelset(dt,s,e,comm1d, nbr_left,  nbr_right)

    use M_MPI_General_3D, only:nx,ny,nz,x,y,z,u,v,w,parallel
    use M_MPI_Exch_3D

    implicit none
    integer,intent(in)                                             ::   s,e
    integer,intent(in)                                             ::   comm1d,nbr_left,nbr_right
    double precision,intent(in)                                    ::   dt

    double precision,dimension (:,:,:),allocatable                 ::   phiold2
    double precision,DIMENSION (:,:,:),allocatable                 ::   Dpphix,Dmphix,Dpphiy,Dmphiy,Dpphiz,Dmphiz
    double precision,DIMENSION (:,:,:),allocatable                 ::   phix,phiy,phiz,Lphin
    double precision                                               ::   Lphis
    integer                                                        ::   ierr
    integer i,j,k,kk
   
   print*, "I am in level set 1"

    allocate( Phiold2(s-2:e+2,-1:ny+1,-1:nz+1) )
    allocate(Dpphix(s-1:e+1,0:ny,0:nz),Dmphix(s-1:e+1,0:ny,0:nz),Dpphiy(s-1:e+1,0:ny,0:nz),Dmphiy(s-1:e+1,0:ny,0:nz), &
       &     Dpphiz(s-1:e+1,0:ny,0:nz),Dmphiz(s-1:e+1,0:ny,0:nz)  )
    allocate(phix(s:e,1:ny-1,1:nz-1),phiy(s:e,1:ny-1,1:nz-1),phiz(s:e,1:ny-1,1:nz-1),Lphin(s:e,1:ny-1,1:nz-1))


    phiold2(s-2:e+2,-1:ny+1,-1:nz+1)=phi(s-2:e+2,-1:ny+1,-1:nz+1)


    do kk=1,2  !!prediction correction method!!

      do i=s,e
        do j=0,ny
          do  k=0,nz
            Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i)  )
            Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1)  )
            Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j)  )
            Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1)  )
            Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k)  )
            Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1)  )
          end do 
        end do 
      end do
   
    if (s.eq.1) then
      i=0
      do j=0,ny
        do  k=0,nz
          Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i)  )
          Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1)  )
          Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j)  )
          Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1)  )
          Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k)  )
          Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1)  )
        end do 
      end do
    end if

    if (e.eq.(nx-1)) then
      i=nx
      do j=0,ny
        do  k=0,nz
          Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i)  )
          Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1)  )
          Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j)  )
          Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1)  )
          Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k)  )
          Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1)  )
        end do 
      end do
    end if

    if (Parallel) then

      call Dy_A_Exch(Dmphix,nx,ny,nz, s, e, comm1d, nbr_left,  nbr_right)
      call Dy_A_Exch(Dpphix,nx,ny,nz, s, e, comm1d, nbr_left, nbr_right)

    end if
 


    do i=s,e 
      do j=1,ny-1 
        do k=1,nz-1   !!A!!



          if (0.5*( u(i,j,k)+u(i+1,j,k) ).gt.0.0) then

            if (  abs(  Dmphix(i,j,k)-Dmphix(i-1,j,k) ).lt.abs(  Dpphix(i,j,k)-Dpphix(i-1,j,k) )   ) then
              phix(i,j,k)=Dmphix(i,j,k)+  0.5*(  Dmphix(i,j,k)-Dmphix(i-1,j,k)  )
            else
              phix(i,j,k)=Dmphix(i,j,k)+  0.5*(  Dpphix(i,j,k)-Dpphix(i-1,j,k)  )
            end if

          else

            if (  abs(  Dmphix(i+1,j,k)-Dmphix(i,j,k) ).lt.abs(  Dpphix(i+1,j,k)-Dpphix(i,j,k) )   ) then
              phix(i,j,k)=Dpphix(i,j,k)-  0.5*(  Dmphix(i+1,j,k)-Dmphix(i,j,k)  )
            else
              phix(i,j,k)=Dpphix(i,j,k)-  0.5*(  Dpphix(i+1,j,k)-Dpphix(i,j,k)  )
            end if

          end if


          if (0.5*( V(i,j,k)+V(i,j+1,k) ).gt.0.0) then

            if (  abs(  DmphiY(i,j,k)-DmphiY(i,j-1,k) ).lt.abs(  DpphiY(i,j,k)-DpphiY(i,j-1,k) )   ) then
              phiY(i,j,k)=DmphiY(i,j,k)+  0.5*(  DmphiY(i,j,k)-DmphiY(i,j-1,k)  )
            else
              phiY(i,j,k)=DmphiY(i,j,k)+  0.5*(  DpphiY(i,j,k)-DpphiY(i,j-1,k)  )
            end if

          else

            if (  abs(  DmphiY(i,j+1,k)-DmphiY(i,j,k) ).lt.abs(  DpphiY(i,j+1,k)-DpphiY(i,j,k) )   ) then
              phiY(i,j,k)=DpphiY(i,j,k)-  0.5*(  DmphiY(i,j+1,k)-DmphiY(i,j,k)  )
            else
              phiY(i,j,k)=DpphiY(i,j,k)-  0.5*(  DpphiY(i,j+1,k)-DpphiY(i,j,k)  )
            end if

          end if


          if (0.5*( W(i,j,k)+W(i,j,k+1) ).gt.0.0) then

            if (  abs(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1) ).lt.abs(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1) )   ) then
              phiZ(i,j,k)=DmphiZ(i,j,k)+  0.5*(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1)  )
            else
              phiZ(i,j,k)=DmphiZ(i,j,k)+  0.5*(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1)  )
            end if

          else

            if (  abs(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k) ).lt.abs(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k) )   ) then
              phiZ(i,j,k)=DpphiZ(i,j,k)-  0.5*(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k)  )
            else
              phiZ(i,j,k)=DpphiZ(i,j,k)-  0.5*(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k)  )
            end if

          end if



        end do 
      end do 
    end do   !!A!!


    if (kk.eq.1) then

      do i=s,e 
        do j=1,ny-1 
          do k=1,nz-1

            Lphin(i,j,k)=   (-0.5)*( u(i,j,k)+u(i+1,j,k) )* phix(i,j,k)+&
                       &    (-0.5)*( V(i,j,k)+V(i,j+1,k) )* phiy(i,j,k)+&
                       &    (-0.5)*( W(i,j,k)+W(i,j,k+1) )* phiz(i,j,k)

            phi(i,j,k)=phiold2(i,j,k)+ &
                     & dt*( (-0.5)*( u(i,j,k)+u(i+1,j,k) )* phix(i,j,k)+&
                     &      (-0.5)*( V(i,j,k)+V(i,j+1,k) )* phiy(i,j,k)+&
                     &      (-0.5)*( W(i,j,k)+W(i,j,k+1) )* phiz(i,j,k) )

          end do 
        end do 
      end do

    end if


    if (kk.eq.2) then

      do i=s,e 
        do j=1,ny-1 
          do k=1,nz-1
!phi(i,j)=phi(i,j)+dt*( -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
 !                    & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j)  )

            Lphis=   (-0.5)*( u(i,j,k)+u(i+1,j,k) )* phix(i,j,k)+&
            &        (-0.5)*( V(i,j,k)+V(i,j+1,k) )* phiy(i,j,k)+&
                   & (-0.5)*( W(i,j,k)+W(i,j,k+1) )* phiz(i,j,k)

            phi(i,j,k)=phiold2(i,j,k)+0.5*dt*(Lphis+Lphin(i,j,k) )
          end do 
        end do 
      end do

    end if

    call Levelset_bound(s,e)

    if (Parallel) then
      call UVPhi_Exch(phi,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
    end if

  end do  !!prediction corection method!!                      



  RETURN
  END subroutine 

  subroutine reinitialize(s,e,comm1d, nbr_left,  nbr_right)

    use M_MPI_General_3D, only:nx,ny,nz,x,y,z,pi,parallel
    use M_MPI_Exch_3D

    implicit none
    integer,intent(in)                                             ::   s,e
    integer,intent(in)                                             ::   comm1d,nbr_left,nbr_right

    double precision, dimension(:,:,:),allocatable                 ::   phiold_Orig,phin,lphin
    double precision, dimension(:,:,:),allocatable                 ::   Dpphix,Dmphix,Dpphiy,Dmphiy,Dpphiz,Dmphiz
    double precision                                               ::   Lphis
    double precision                                               ::   phixm,phiym,phizm,phixp,phiyp
    double precision                                               ::   phizp,phixr,phiyr,phizr,sphi,ss

    integer                                                        ::   ierr
    integer i,j,k,kk,m

    allocate( Phiold_Orig(s-2:e+2,-1:ny+1,-1:nz+1),Phin(s-2:e+2,-1:ny+1,-1:nz+1),lphin(s-2:e+2,-1:ny+1,-1:nz+1) )
    allocate( Dpphix(s-1:e+1,0:ny,0:nz),Dmphix(s-1:e+1,0:ny,0:nz),Dpphiy(s-1:e+1,0:ny,0:nz),Dmphiy(s-1:e+1,0:ny,0:nz), &
       &      Dpphiz(s-1:e+1,0:ny,0:nz),Dmphiz(s-1:e+1,0:ny,0:nz)  )


   phiold_Orig(s-2:e+2,-1:ny+1,-1:nz+1)=phi(s-2:e+2,-1:ny+1,-1:nz+1)
   do kk=1,3 !? 






     do m=1,2   !!runge kutta!!





       do i=s,e 
         do j=0,ny 
           do k=0,nz
             Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i) )
             Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1) )
             Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j) )
             Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1) )
             Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k) )
             Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1) )
           end do 
         end do 
       end do
      
       if (s==1) then
         i=0 
         do j=0,ny 
           do k=0,nz
             Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i) )
             Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1) )
             Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j) )
             Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1) )
             Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k) )
             Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1) )
           end do 
         end do
       end if
  
       if (e==(nx-1)) then
         i=nx 
         do j=0,ny 
           do k=0,nz
             Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i) )
             Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1) )
             Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j) )
             Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1) )
             Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k) )
             Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1) )
           end do 
         end do
       end if  

       if (Parallel) then
         call Dy_A_Exch(Dmphix,nx,ny,nz, s, e, comm1d, nbr_left,  nbr_right)
         call Dy_A_Exch(Dpphix,nx,ny,nz, s, e, comm1d, nbr_left,  nbr_right)
       end if



       do i=s,e 
         do j=1,ny-1
           do k=1,nz-1 !!A!!  !!A!!

             !sphi=phiold(i,j)/(  sqrt(phiold(i,j)*phiold(i,j)+h*h)  )
             !sphi=sgn(phiold)

             if (phiold_Orig(i,j,k).gt.eps) then
               sphi=1.0
             else if (phiold_Orig(i,j,k).lt.-eps) then
               sphi=-1.0
             else
               sphi=phiold_Orig(i,j,k)/eps -(1/pi)*sin(pi*phiold_Orig(i,j,k)/eps)
             end if



             if (  abs(  Dmphix(i,j,k)-Dmphix(i-1,j,k) ).lt.abs(  Dpphix(i,j,k)-Dpphix(i-1,j,k) )   ) then
               phixm=Dmphix(i,j,k)+  0.5*(  Dmphix(i,j,k)-Dmphix(i-1,j,k)  )
             else
               phixm=Dmphix(i,j,k)+  0.5*(  Dpphix(i,j,k)-Dpphix(i-1,j,k)  )
             end if

             if (  abs(  Dmphix(i+1,j,k)-Dmphix(i,j,k) ).lt.abs(  Dpphix(i+1,j,k)-Dpphix(i,j,k) )   ) then
               phixp=Dpphix(i,j,k)-  0.5*(  Dmphix(i+1,j,k)-Dmphix(i,j,k)  )
             else
               phixp=Dpphix(i,j,k)-  0.5*(  Dpphix(i+1,j,k)-Dpphix(i,j,k)  )
             end if


             if (  abs(  DmphiY(i,j,k)-DmphiY(i,j-1,k) ).lt.abs(  DpphiY(i,j,k)-DpphiY(i,j-1,k) )   ) then
               phiYm=DmphiY(i,j,k)+  0.5*(  DmphiY(i,j,k)-DmphiY(i,j-1,k)  )
             else
               phiYm=DmphiY(i,j,k)+  0.5*(  DpphiY(i,j,k)-DpphiY(i,j-1,k)  )
             end if

             if (  abs(  DmphiY(i,j+1,k)-DmphiY(i,j,k) ).lt.abs(  DpphiY(i,j+1,k)-DpphiY(i,j,k) )   ) then
               phiYp=DpphiY(i,j,k)-  0.5*(  DmphiY(i,j+1,k)-DmphiY(i,j,k)  )
             else
               phiYp=DpphiY(i,j,k)-  0.5*(  DpphiY(i,j+1,k)-DpphiY(i,j,k)  )
             end if


             if (  abs(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1) ).lt.abs(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1) )   ) then
               phiZm=DmphiZ(i,j,k)+  0.5*(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1)  )
             else
               phiZm=DmphiZ(i,j,k)+  0.5*(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1)  )
             end if

             if (  abs(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k) ).lt.abs(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k) )   ) then
               phiZp=DpphiZ(i,j,k)-  0.5*(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k)  )
             else
               phiZp=DpphiZ(i,j,k)-  0.5*(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k)  )
             end if



            if (sphi*phixp.ge.0.0.AND.sphi*phixm.ge.0.0) then
              phixr=phixm
            else if (sphi*phixp.le.0.0.AND.sphi*phixm.le.0.0) then
              phixr=phixp
            else if (sphi*phixp.gt.0.0.AND.sphi*phixm.lt.0.0) then
              phixr=0.0
            else if (sphi*phixp.lt.0.0.AND.sphi*phixm.gt.0.0) then
              ss=sphi*( abs(phixp)-abs(phixm) )/(phixp-phixm)
              if (ss.gt.0.0) then
                phixr=phixm
              else
                phixr=phixp
             end if
            end if

            if (sphi*phiyp.ge.0.0.AND.sphi*phiym.ge.0.0) then
              phiyr=phiym
            else if (sphi*phiyp.le.0.0.AND.sphi*phiym.le.0.0) then
              phiyr=phiyp
            else if (sphi*phiyp.gt.0.0.AND.sphi*phiym.lt.0.0) then
              phiyr=0.0
            else if (sphi*phiyp.lt.0.0.AND.sphi*phiym.gt.0.0) then
              ss=sphi*( abs(phiyp)-abs(phiym) )/(phiyp-phiym)
              if (ss.gt.0.0) then
                phiyr=phiym
              else
                phiyr=phiyp
              end if
            end if

            if (sphi*phizp.ge.0.0.AND.sphi*phizm.ge.0.0) then
              phizr=phizm
            else if (sphi*phizp.le.0.0.AND.sphi*phizm.le.0.0) then
              phizr=phizp
            else if (sphi*phizp.gt.0.0.AND.sphi*phizm.lt.0.0) then
              phizr=0.0
            else if (sphi*phizp.lt.0.0.AND.sphi*phizm.gt.0.0) then
              ss=sphi*( abs(phizp)-abs(phizm) )/(phizp-phizm)
              if (ss.gt.0.0) then
                phizr=phizm
              else
                phizr=phizp
              end if
            end if


            if (m.eq.1) then
              lphin(i,j,k)=sphi*(  1.0-sqrt(phiyr*phiyr+phixr*phixr+phizr*phizr)  )
              phin(i,j,k)=phi(i,j,k)
              phi(i,j,k)=phi(i,j,k)+dtau*sphi*(  1.0-sqrt(phiyr*phiyr+phixr*phixr+phizr*phizr)  )
            end if

            if (m.eq.2) then
              lphis   =sphi*(  1.0-sqrt(phiyr*phiyr+phixr*phixr+phizr*phizr)  )
              phi(i,j,k)=phin(i,j,k)+0.5*dtau*(  lphis+lphin(i,j,k)  )
            end if




          end do 
        end do 
      end do

      call Levelset_bound(s,e)

      if (Parallel) then
        call UVPhi_Exch(phi,nx,ny,nz,s, e, comm1d, nbr_left,  nbr_right)
      end if





   end do






 end do 

 return
END subroutine 





  subroutine Levelset_bound(s,e)
    use M_MPI_General_3D, only: nx,ny,nz

    implicit none
    integer,intent(in)              :: s,e
    integer i,j,k


    if (s==1) then 
      do j=1,ny-1 
        do k=1,nz-1
          phi(0,j,k)=2*phi(1,j,k)-phi(2,j,k)
          phi(-1,j,k)=3*phi(1,j,k)-2*phi(2,j,k)
        end do
      end do
    end if
  
    if (e==(nx-1)) then 
      do j=1,ny-1
        do k=1,nz-1   
          phi(nx,j,k)=2*phi(nx-1,j,k)-phi(nx-2,j,k)
          phi(nx+1,j,k)=3*phi(nx-1,j,k)-2*phi(nx-2,j,k)
        end do 
      end do
    end if


    
    do i=s,e 
      do k=1,nz-1
        phi(i,0,k)=2*phi(i,1,k)-phi(i,2,k)
        phi(i,-1,k)=3*phi(i,1,k)-2*phi(i,2,k)
        phi(i,ny,k)=2*phi(i,ny-1,k)-phi(i,ny-2,k)
        phi(i,ny+1,k)=3*phi(i,ny-1,k)-2*phi(i,ny-2,k)
      end do 
    end do

    if (s==1) then
      do i=-1,0 
        do k=1,nz-1
          phi(i,0,k)=2*phi(i,1,k)-phi(i,2,k)
          phi(i,-1,k)=3*phi(i,1,k)-2*phi(i,2,k)
          phi(i,ny,k)=2*phi(i,ny-1,k)-phi(i,ny-2,k)
          phi(i,ny+1,k)=3*phi(i,ny-1,k)-2*phi(i,ny-2,k)
        end do 
      end do
    end if
 
    if (e==(nx-1)) then
      do i=nx,nx+1 
        do k=1,nz-1
          phi(i,0,k)=2*phi(i,1,k)-phi(i,2,k)
          phi(i,-1,k)=3*phi(i,1,k)-2*phi(i,2,k)
          phi(i,ny,k)=2*phi(i,ny-1,k)-phi(i,ny-2,k)
          phi(i,ny+1,k)=3*phi(i,ny-1,k)-2*phi(i,ny-2,k)
        end do 
      end do
    end if 
  

 
    do i=s,e 
      do j=-1 , ny+1
        phi(i,j,0)=2*phi(i,j,1)-phi(i,j,2)
        phi(i,j,-1)=3*phi(i,j,1)-2*phi(i,j,2)
        phi(i,j,nz)=2*phi(i,j,nz-1)-phi(i,j,nz-2)
        phi(i,j,nz+1)=3*phi(i,j,nz-1)-2*phi(i,j,nz-2)
      end do 
    end do

    if (s==1) then 
      do i=-1,0 
        do j=-1 , ny+1
          phi(i,j,0)=2*phi(i,j,1)-phi(i,j,2)
          phi(i,j,-1)=3*phi(i,j,1)-2*phi(i,j,2)
          phi(i,j,nz)=2*phi(i,j,nz-1)-phi(i,j,nz-2)
          phi(i,j,nz+1)=3*phi(i,j,nz-1)-2*phi(i,j,nz-2)
        end do 
      end do
    end if

    if (e==(nx-1)) then
      do i=nx-1,nx
        do j=-1 , ny+1
          phi(i,j,0)=2*phi(i,j,1)-phi(i,j,2)
          phi(i,j,-1)=3*phi(i,j,1)-2*phi(i,j,2)
          phi(i,j,nz)=2*phi(i,j,nz-1)-phi(i,j,nz-2)
          phi(i,j,nz+1)=3*phi(i,j,nz-1)-2*phi(i,j,nz-2)
        end do 
      end do
    end if 

    return

  end subroutine 












end module
