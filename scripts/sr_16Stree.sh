#! /bin/bash

#SBATCH --time=04-00:00:00
#SBATCH --partition=defq
#SBATCH --mem=192GB
#SBATCH --mail-user=myemail@myemail.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --job-name=sr_16S
#SBATCH --comment=sr_16S

module load mamba
mamba activate sr_16S

DIR="/home/tmhagm8/scratch/fmtad_wabx/data"

mafft --maxiterate 1000 --genafpair --thread 64 "$DIR/rep-seqs.fasta" > "$DIR/rep-seqs-aln.fasta"

FastTree -nt -gtr "$DIR/rep-seqs-aln.fasta" > "$DIR/rep-seqs.tree"
