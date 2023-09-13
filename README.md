# GEMX
Gyrokinetic simulation of the tokamak edge including the scrape-off layer

Please add any useful notes to this file. Place it anywhere and we can clean up later.

GEM_X-3.pdf is our equations or user manual.  Consider placing any notes on implementations in this overleaf document: https://www.overleaf.com/project/64456a5749a720661d4d1f0a

Versions:

Sept 13, 2023 
Have GK Poisson working and particle tracers in equilibrium B working.  Zhichen is working getting Alfven waves.  Field solve needs to be parallelized and optimized. Qiheng and Scott are cleaning up wrappers, adding diagnostic routines and .ipynb diagnostic scripts.

-----------------------Cheet Sheet stuff follows------------------------------------

tar cvf * - When taring up GEMX, do not include top directory, just the files and any subdirectories.

tar xvf file.tar

NERSC user names:
cqh2021 - Qiheng
u10198 - Yang
nonohing - Zhichen
jycheng - Junyi

give -u <username> file.tar
take -u <username> -a

source env.sh

sbatch submit.cmd
sacct
sqs
scancel $jobid
scancel -u $USER
