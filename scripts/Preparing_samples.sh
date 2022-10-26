#!/bin/bash
#SBATCH -p himem
#SBATCH --mem=30G
#SBATCH -t 12:00:00
#SBATCH -J Preparing_CRPC_files_and_consensus_set
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -c 1

#cd $SLURM_SUBMIT_DIR

#----------------------------------------
# Load dependencies
#----------------------------------------

module load R/4.1.0
module load bedtools/2.30.0

mkdir peaks_progress bedfiles

###Update these two
narrowpeakpath="/mnt/work1/users/lupiengroup/People/ankita/CRPC/mapping/data/old-run/peaks"
narrowpeakextension=".filtered.narrowPeak.gz"

######

for f in "$narrowpeakpath"/*"$narrowpeakextension" ;
do
  filename=$(echo "${f##*/}" | sed 's/\'$narrowpeakextension'[^.]*$//')
  zcat "$f" | awk 'BEGIN {FS=OFS="\t"}{print $1,$2,$3}' > peaks_progress/"$filename".bed
  sortBed -i peaks_progress/"$filename".bed > bedfiles/"$filename".narrowPeak.sorted.bed
done
