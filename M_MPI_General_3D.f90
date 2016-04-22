module M_MPI_General_3D

  implicit none 
  save
  logical     , parameter                           :: parallel=.true.
  integer     , parameter                           :: nx=241, ny=101, nz=151   !nx=191, ny=61, nz=91  !nx=211, ny=60, nz=71 !nx=191,ny=60  ,nz=91 !  ny=60  nz=91
  Double precision     , parameter                           :: pi=3.141592653
  Double precision     , parameter                           :: Froude=(1.0/12.0)
  double precision     , parameter                           :: Landa=Froude*156 
  Double precision,dimension (:,:,:),allocatable        ::   u,v,w
  Double precision,dimension (:),allocatable            ::   x,hx
  Double precision,dimension (-1:ny+1)                  ::   y,hy
  Double precision,dimension (-1:nz+1)                  ::   z,hz

  Double precision,dimension (1,3)                      ::   Bp
  Double precision,dimension (3)                        ::   Chv

  Double precision                                      ::   lx,ly,lz,dt,gx,gy,gz
  Double precision                                      ::   beta,maxSOR,pdif,toler,maxError
  Double precision                                      ::   miudrop,rodrop,miuAir,roAir
  Double precision                                      ::   ZFree


  integer                                      ::   tstepE
  integer                                      ::   maxIteration,plot


  contains
   

    Subroutine General_Constant_Ini(s,e)
      implicit none
      integer,intent(in)           ::s,e
  
 
      Include "Par_Constant_General.txt"    

      allocate (u(s-2:e+2,-1:ny+1,-1:nz+1),v(s-2:e+2,-1:ny+1,-1:nz+1),w(s-2:e+2,-1:ny+1,-1:nz+1))
      allocate (x(s-2:e+2),hx(s-2:e+2))
      u(s-2:e+2,-1:ny+1,-1:nz+1)=0
      v(s-2:e+2,-1:ny+1,-1:nz+1)=0
      w(s-2:e+2,-1:ny+1,-1:nz+1)=0

      return 
    end subroutine 




end module
