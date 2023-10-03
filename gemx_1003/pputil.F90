!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE pputil
!
      use mpi
      IMPLICIT NONE
  PRIVATE
  PUBLIC :: ppinit
!
  INTEGER, SAVE :: me, nvp,npp,GCLR,TCLR

  INTEGER, SAVE :: TUBE_COMM,GRID_COMM

CONTAINS
!======================================================================
!
 	SUBROUTINE ppinit(idproc,nproc,ntube,com1,com2)
!
     use mpi
     INTEGER, INTENT(IN) :: ntube
     INTEGER, INTENT(OUT) :: nproc
	INTEGER, INTENT(OUT) :: idproc,com1,com2
	INTEGER :: ierr,npp
!
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npp, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
	nproc = npp
	IDPROC = me
! 
	GCLR=INT(me/ntube)
	TCLR=MOD(me,ntube)
!
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,GCLR,&
     	 &	TCLR,GRID_COMM,ierr)
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,TCLR,&
         &	GCLR,TUBE_COMM,ierr)
!
	com1=TUBE_COMM
	com2=GRID_COMM
	nvp=npp/ntube
!
	END SUBROUTINE ppinit
END MODULE pputil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
