#!/bin/bash
#SBATCH --account=PCON0003
#SBATCH --job-name=omega_test
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20  # Use up to 13 cores for parallel processing
#SBATCH --output=omega_test_%j.out  # Output file
#SBATCH --error=omega_test_%j.err   # Error file

# Load the necessary Python module
module load python/3.9-2022.05

# Copy necessary files to the local node's temporary directory
cp makeintegrand.py Constants.fits I_functions_class.py module1.py merge_fits.py $TMPDIR

# Change to the temporary directory
cd $TMPDIR

# Set x within the script
x=200  # Example value for x
XSTEP=0.01  # Define XSTEP value

# Calculate END_OMEGA based on the formula
END_OMEGA=$(echo "scale=2; 20.0 - ($x * $XSTEP)" | bc)

# Define the step size for omega_scale
STEP=0.1

# Generate points between 0.1 and END_OMEGA with the defined STEP
omega_scales=$(seq 0.1 $STEP $END_OMEGA)

# Run makeintegrand.py for each omega scale in parallel
for omega_scale in $omega_scales; do
    python -u makeintegrand.py $x $omega_scale &  # Run each task in the background
    
    # Manually check if there are 13 jobs running (since 13 cores are available)
    while (( $(jobs | wc -l) >= 20 )); do
        sleep 1  # Wait for 1 second before checking again
    done
done

# Wait for all jobs to complete
wait

# Run the merging script to merge the FITS files
python -u merge_fits.py $x

# Copy the merged FITS file back to the submission directory
cp merged_output_x${x}.fits $SLURM_SUBMIT_DIR




