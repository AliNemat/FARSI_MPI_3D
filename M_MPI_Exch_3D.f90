
Module M_MPI_Exch_3D
use mpi
implicit none

contains


!=================================================================================================
!=================================================================================================
!=================================================================================================
! creates new data type for passing non-contiguous 
!-------------------------------------------------------------------------------------------------
!SUBROUTINE Create_Face_DataType(nx1,ny1,nz1, inewtype)
!  implicit none
!  integer ::  inewtype,isize, ierr, itemp, idist

!  CALL MPI_TYPE_EXTENT(MPI_real, isize, ierr)

!  CALL MPI_TYPE_VECTOR(ny+3, 1,e-s+1, MPI_real, itemp, ierr)

!  idist = (nx+3) * (ny+3) * isize
!  CALL MPI_TYPE_HVECTOR(nz+3, 1, idist, itemp, inewtype, ierr)
!  CALL MPI_TYPE_COMMIT(inewtype, ierr)
!END subroutine 

subroutine Dy_A_Exch(Dy_A,nx1,ny1,nz1,s, e, comm1d, nbr_left, nbr_right)

  implicit none
  integer,intent(in)                                                 ::  nx1,ny1,nz1,s,e
  Double precision,dimension(s-1:e+1, 0:ny1,0:nz1),intent(inout)         ::  Dy_A
  integer,intent(in)                                                 ::  comm1d, nbr_left, nbr_right
  integer                                                            ::  status(MPI_STATUS_SIZE), ierr
  integer                                                            ::  isize,  itemp, idist
  logical,save                                                       ::  First_Time=.true.
  integer,save                                                       :: FaceX


!  print*, "I am in adevction exch"

  if (First_Time) then 
    First_time=.false.
    CALL MPI_TYPE_EXTENT(MPI_double_precision, isize, ierr)

    CALL MPI_TYPE_VECTOR(ny1+1, 1,e-s+3, MPI_double_precision, itemp, ierr)

    idist = (e-s+3) * (ny1+1) * isize
    CALL MPI_TYPE_HVECTOR(nz1+1, 1, idist, itemp, FaceX, ierr)
    CALL MPI_TYPE_COMMIT(FaceX, ierr)
  !CALL MPI_TYPE_COMMIT(itemp, ierr)
  end if 
      
     call MPI_SENDRECV( &
     &            Dy_A(e,0,0)  ,1, FaceX, nbr_right,      0, &
     &            Dy_A(s-1,0,0)  ,1,FaceX, nbr_left,     0, &
     &            comm1d, status, ierr )

     call MPI_SENDRECV( &
     &            Dy_A(s,0,0)  ,1, FaceX, nbr_left,      1, &
     &            Dy_A(e+1,0,0),1, FaceX, nbr_right,     1, &
     &            comm1d, status, ierr )


  !  CALL MPI_TYPE_Free(FaceX, ierr)
  return
end subroutine


subroutine UVPhi_Exch(UVPhi,nx1,ny1,nz1,s, e, comm1d, nbr_left, nbr_right)

  implicit none
  integer,intent(in)                                                 ::  nx1,ny1,nz1,s,e
  Double precision,dimension(s-2:e+2, -1:ny1+1,-1:nz1+1),intent(inout)         ::  UVPhi
  integer,intent(in)                                                 ::  comm1d, nbr_left, nbr_right
  integer                                                            ::  status(MPI_STATUS_SIZE), ierr

  integer                                                            ::   isize,  itemp, idist
  logical,save                                                       ::  First_Time=.true.
  integer,save                                                            :: FaceX2

  !print*, "I am in velocity exch"

  if (First_Time) then 
    First_time=.false.
    CALL MPI_TYPE_EXTENT(MPI_double_precision, isize, ierr)

    CALL MPI_TYPE_VECTOR(ny1+3, 1,e-s+5, MPI_double_precision, itemp, ierr)

    idist = (e-s+5) * (ny1+3) * isize
    CALL MPI_TYPE_HVECTOR(nz1+3, 1, idist, itemp, FaceX2, ierr)
    CALL MPI_TYPE_COMMIT(FaceX2, ierr)
    !CALL MPI_TYPE_COMMIT(itemp, ierr)
    !print*, FaceX, "Facex"
      !  if (e==(nx-1)) then
          ! if (s==1) then
            !print*,Dy_A(3,0,0),s,e
            ! print*,dmux(3,3,0),s,e
           !end if

           !if (s==1) then
           !if (e==(nx1-1)) then
            ! print*,dy_A(4,0,0),s,e
            ! print*,dmux(4,3,0),s,e
           !end if

   end if    
     call MPI_SENDRECV( &
     &            UVPhi(e,-1,-1)  ,1, FaceX2, nbr_right,      0, &
     &            UVPhi(s-1,-1,-1)  ,1,FaceX2, nbr_left,     0, &
     &            comm1d, status, ierr )

     call MPI_SENDRECV( &
     &            UVPhi(s,-1,-1)  ,1, FaceX2, nbr_left,      1, &
     &            UVPhi(e+1,-1,-1),1, FaceX2, nbr_right,     1, &
     &            comm1d, status, ierr )


   ! CALL MPI_TYPE_Free(FaceX2, ierr)
  return
end subroutine



subroutine Xhx_Exch(Xhx, s, e, comm1d, nbr_left, nbr_right)
  implicit none
  integer,intent(in)                                             :: s, e
  double precision,dimension (s-2:e+2),intent(inout)             :: Xhx
  integer,intent(in)                                             :: comm1d, nbr_left, nbr_right
  integer                                                        :: status(MPI_STATUS_SIZE), ierr


     call MPI_SENDRECV( &
     &            Xhx(e)  ,1, MPI_DOUBLE_PRECISION, nbr_right,    0, &
     &            Xhx(s-1),1, MPI_DOUBLE_PRECISION, nbr_left, 0, &
     &            comm1d, status, ierr )

     call MPI_SENDRECV( &
     &            Xhx(s)  ,1, MPI_DOUBLE_PRECISION, nbr_left, 1, &
     &            Xhx(e+1),1, MPI_DOUBLE_PRECISION, nbr_right, 1, &
     &            comm1d, status, ierr )


  return
end subroutine

subroutine XE_Exch(XE, s, e, comm1d, nbr_left, nbr_right)
  implicit none
  integer,intent(in)                                             :: s, e
  double precision,dimension (s-1:e),intent(inout)               :: XE
  integer,intent(in)                                             :: comm1d, nbr_left, nbr_right
  integer                                                        :: status(MPI_STATUS_SIZE), ierr


     call MPI_SENDRECV( &
     &            XE(e)  ,1, MPI_DOUBLE_PRECISION, nbr_right,    0, &
     &            XE(s-1),1, MPI_DOUBLE_PRECISION, nbr_left ,    0, &
     &            comm1d, status, ierr )

   return

  end subroutine 


end module 
