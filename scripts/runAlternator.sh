#!/bin/bash

#-----------------------------------------------------------------
# Pure MPI job , using 80 procs on 4 nodes ,
# with 20 procs per node and 1 thread per MPI task

#-----------------------------------------------------------------

#SBATCH --job-name=ACGeneratorTest # Job name
#SBATCH --output=ACGeneratorTest.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=ACGeneratorTest.%j.err # Stderr (%j expands to jobId)
#SBATCH --ntasks=80 # Total number of tasks
#SBATCH --nodes=4 # Total number of nodes requested
#SBATCH --ntasks-per-node=20 # Tasks per node
#SBATCH --cpus-per-task=1 # Threads per task(=1) for pure MPI
#SBATCH --mem=56000 # Memory per job in MB
#SBATCH -t 01:30:00 # Run time (hh:mm:ss) - (max 48h)
#SBATCH --partition=compute # Submit queue
#SBATCH -A ACGenerator # Accounting project

# Load any necessary modules

module load intel
module load intelmpi

export I_MPI_FABRICS=shm:dapl

# Launch the executable

srun dotnet test --filter "MGroup.FEM.Thermal.Tests.ElectricalEquipment.AlternatorTest.Test1"
