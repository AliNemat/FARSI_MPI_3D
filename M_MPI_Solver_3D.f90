Module M_MPI_Solver_3D
  implicit none 
  private                                                  ::    DIFF
  contains

subroutine Poisson_Solver_Parallel_3D(Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,pdif,beta,P,s,e,comm1d,nbr_left,nbr_right)

!Parallel
    use mpi
!Parallel
    use M_MPI_General_3D, only: nx,ny,nz
    implicit none  
    integer,intent(in)                                                 ::    s,e
    integer,intent(in)                                                 ::    comm1d,nbr_left,nbr_right
    double precision,intent(in)                                        ::    pdif,beta
    double precision,dimension (s:e,1:ny-1,1:nz-1)  ,intent(in)        ::    Apx,Amx,Apy,Amy,Apz,Amz,AP,Q
    double precision,dimension (s-1:e+1,0:ny,0:nz)  ,intent(inout)     ::    p

    double precision                                                   ::    diffnorm,dwork
    double precision,dimension (s-1:e+1,0:ny,0:nz)                     ::    p1,p2
    integer i,j,k,it
    integer                                                            ::    status(MPI_STATUS_SIZE), ierr

    print*,"I am in Poisson parallel solver"
    p1(s-1:e+1,0:ny,0:nz)=p(s-1:e+1,0:ny,0:nz)
    p2(s-1:e+1,0:ny,0:nz)=p(s-1:e+1,0:ny,0:nz)
    do it=1,10000


      call Sweep1D   (Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,beta,p1,p2,s,e,ny,nz)
      call PTemp_Exch( p2, s, e,ny,nz, comm1d, nbr_left, nbr_right )  
      call Sweep1D   (Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,beta,p2,p1,s,e,ny,nz)  
      call PTemp_Exch( p1, s, e,ny,nz, comm1d, nbr_left, nbr_right )   
  
       dwork = DIFF   ( p1, p2, s, e,ny,nz )
    !   print*, dwork,"difference captured"
      call MPI_allreduce( dwork, diffnorm, 1, MPI_DOUBLE_PRECISION,MPI_SUM, comm1d, ierr )
 
      ! diffnorm=dwork

       if (diffnorm .lt. pdif) then
         exit
       end if


    end do 

    if (s.eq.1) then
      print *, 'Converged after ', 2*it, 'Iterations','diff',diffnorm
    end if

    p(s-1:e+1,0:ny,0:nz)=p1(s-1:e+1,0:ny,0:nz)

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

  print*, "I am finished poisson "
     return
   end subroutine 




  Subroutine Sweep1D (Amx,Apx,Amy,Apy,Amz,Apz,Ap,Q,beta,pold,pnew,s,e,ny1,nz1)  
    implicit none

    integer,intent(in)                                                         ::    s,e,ny1,nz1
    double precision,intent(in)                                                ::    beta
    double precision,dimension (s:e,1:ny1-1,1:nz1-1)    ,intent(in)              ::    Apx,Amx,Apy,Amy,Apz,Amz,AP,Q
    double precision,dimension (s-1:e+1,0:ny1,0:nz1)    ,intent(in)              ::    pold
    double precision,dimension (s-1:e+1,0:ny1,0:nz1)    ,intent(out)             ::    pnew


    integer i,j,k

    do i=s,e
      do j=1,ny1-1
        do k=1,nz1-1 
                              
          pnew(i,j,k)=beta*( ( Apx(i,j,k)*Pold(I+1,j,k)+Amx(i,j,k)*Pold(I-1,j,k)+&
                               Apy(i,j,k)*Pold(I,j+1,k)+Amy(i,j,k)*Pold(I,j-1,k)+&
                               Apz(i,j,k)*Pold(I,j,k+1)+Amz(i,j,k)*Pold(I,j,k-1)-&
                               Q(I,j,k) ) /(-AP(i,j,k))  ) &
                & +(1-beta) *pold(i,j,k)

        end do 
      end do 
    end do



    return

  end subroutine


subroutine PTemp_Exch( pex,s, e,ny1,nz1, comm1d, nbr_left, nbr_right )
!Parallel
  use mpi
!Parallel
  implicit none
  integer,intent(in)                                                  :: s, e,ny1,nz1
  double precision,dimension(s-1:e+1,0:ny1,0:nz1),intent(inout)       :: Pex
  integer,intent(in)                                                  :: comm1d, nbr_left, nbr_right
  integer                                                             :: status(MPI_STATUS_SIZE), ierr,isize,itemp,idist
  logical,save                                                        :: First_Time=.true.
  integer,save                                                            :: FacexP

   if (first_time) then

      first_time=.false.
      CALL MPI_TYPE_EXTENT(MPI_double_precision, isize, ierr)
      CALL MPI_TYPE_VECTOR(ny1+1, 1,e-s+3, MPI_double_precision, itemp, ierr)
      idist = (e-s+3) * (ny1+1) * isize
      CALL MPI_TYPE_HVECTOR(nz1+1, 1, idist, itemp, FaceXp, ierr)
      CALL MPI_TYPE_COMMIT(FaceXp, ierr)
    print*,"Derived data type for poisson solver created"
    end if 

     call MPI_SENDRECV(  &
     &            pex(e,0,0)  , 1, FaceXp, nbr_right,    0, &
     &            pex(s-1,0,0), 1, FaceXp, nbr_left,     0, &
     &            comm1d, status, ierr )

     call MPI_SENDRECV( &
     &            pex(s,0,0),   1, FaceXp, nbr_left,     1, &
     &            pex(e+1,0,0), 1, FaceXp, nbr_right,    1, &
     &            comm1d, status, ierr )

    !print*,"1st Data exchange in Poisson solver succesfull "
    !  CALL MPI_TYPE_Free(FaceXp, ierr)
  return


end subroutine 





  double precision function DIFF( a, b, s, e,ny1,nz1 )
      integer s,e,ny1,nz1
      double precision a(s-1:e+1,0:ny1,0:nz1), b(s-1:e+1,0:ny1,0:nz1)
      double precision sump
      integer i, j,k

      sump = 0.0d0
      do i=s,e
        do j=1,ny1-1
          do k=1,nz1-1
            sump = sump + (a(i,j,k) - b(i,j,k)) ** 2
          end do
        end do
      end do 

      diff = sump
      
     ! print*, "I am in the summation function"
    return

  end function


end module 

