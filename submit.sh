#PBS -N experiment  # job name
#PBS -o output/  # output path
#PBS -j oe  # merge error messages with output
#PBS -l ncpus=48  # limit total number of cpus

export num_threads=4 # number of threads to use for each process
export num_processes=12 # total number of processes to use. num_treads * num_processes = ncpus.

# Set Gurobi environment variables
export GUROBI_HOME="/opt/gurobi901/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

# Copy and paste this block into all STORM submission scripts. Limits thread usage of various packages by num_threads.
export MKL_NUM_THREADS=$num_threads      # Limits Intel math kernel
export OPENBLAS_NUM_THREADS=$num_threads # Limits OPENBLAS, this is the most important one
export MC_CORES=$num_threads             # Limits some packages in R
export OMP_NUM_THREADS=$num_threads      # Limits OpenMP
export NUMEXPR_NUM_THREADS=$num_threads  # Limits NumExpr in python

cd ~/robust-mrcpsp
python parallel_experiment.py -instance_dir=~/robust-mrcpsp/instances/j10.mm -solve_method=compact_reformulation
  -Gamma=0 -time_limit=3600 -num_threads=$num_threads -num_processes=$num_processes
