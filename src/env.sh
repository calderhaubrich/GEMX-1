module restore
module load nvhpc
export PETSC_DIR=/global/cfs/cdirs/mp118/software/petsc_new
export PETSC_ARCH=perlmutter-nvhpc
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH