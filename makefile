

NICE    = ifort
Gfort   = gfortran44
FCOMP   = mpifort
F77     = gfortran
CCOMP   = mpicc
CXXCOMP = mpicxx

HYPRE_DIR = /usr/local/lib #HYPRE_DIR = /afs/crc.nd.edu/x86_64_linux/scilib/hypre/2.0.0/intel
HYPRE_LIBS =  -L$(HYPRE_DIR)-lHYPRE #-lHYPRE_struct_ls -lHYPRE_struct_mv -lHYPRE_krylov -lHYPRE_utilities -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv

FOPTS   = #-shared-intel -mcmodel=medium
CFLAGS  = -O3
#-----------------------------------------------------------------------------------------------

EXE = Fsrf2_Strech_15May

OBJ =       M_MPI_General_3D.o   M_MPI_Exch_3D.o   M_MPI_Mesh_3D.o       M_MPI_Solver_3D.o   M_MPI_FrSrf_3D.o   Code_2Phase_3D.o 
INC =       M_MPI_General_3D.f90 M_MPI_Exch_3D.f90 M_MPI_Mesh_3D.f90     M_MPI_Solver_3D.f90 M_MPI_FrSrf_3D.f90 Code_2Phase_3D.f90 

$(EXE): $(OBJ)
	$(FCOMP)  -o $(EXE) $(OBJ) $(CFLAGS)

$(OBJ): $(INC)
        
%.o:	%.f90 
	$(FCOMP) $(FOPTS) -c $<  $(CFLAGS)

clean:
	rm -f *.o *.mod *__genmod.f90

