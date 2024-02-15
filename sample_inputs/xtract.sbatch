#!/bin/tcsh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------
#SBATCH -J xtr_0.0_1_1024
#SBATCH -o xtr_0.0_1_1024.o%j
#SBATCH -e xtr_0.0_1_1024.e%j
#SBATCH -p normal          # Submit to the 'normal' or 'development' queue
#SBATCH -t 48:00:00             # Run time (hh:m) - 1.5 hours
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --tasks-per-node 64
#SBATCH --mail-user=rabatin.b@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
# The next line is required if the user has more than one project
#SBATCH -A TG-AST140008      # Allocation name to charge job against

#USAGE 
#Fix the node count and time above.
#Make sure enzo is built and in the working directory.
#Make sure the parameter field on the ibrun  line.

set echo
mpirun main input_xtr_logK
end
