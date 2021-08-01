### Read mapping pipeline
## Modified from the Matz lab tagSeq pipeline

### Step 1 Concatenate the host and sym transcriptome
# Made dir Sarahs_scripts/final_cleaned_contigs/Mapping_collapsed

# Copied over final transcriptomes
cp ../host_transc/collapsed/O_arbuscula_transcriptome/* ./
cp ../sym_transc/collapsed/B_psygmophilum_transcriptome/* ./

# Concatenate both transcriptomes and respective files
cat B_psygmophilum_transcriptome.fasta >> O_arbuscula_transcriptome.fasta
cat B_psygmophilum_isogroup_to_genename.tab >> O_arbuscula_isogroup_to_genename.tab
cat B_psygmophilum_isogroup_to_GOterm.tab >> O_arbuscula_isogroup_to_GOterm.tab
cat B_psygmophilum_sequenceID_to_isogroup.tab >> O_arbuscula_sequenceID_to_isogroup.tab

# Rename files
mv O_arbuscula_transcriptome.fasta combined_transcriptome.fasta
mv O_arbuscula_isogroup_to_GOterm.tab combined_iso2go.tab
mv O_arbuscula_isogroup_to_genename.tab combined_iso2gene.tab
mv O_arbuscula_sequenceID_to_isogroup.tab combined_seq2iso.tab

### Step 2 concatenate raw forward and reverse reads across the two lanes per sample
# In raw reads Directory files are named OC[genet][symbiontstate]_[lane number]_[1 if forward, 2 if reverse].fq
# this does the equivalent of cat sample1_lane2_forwardreads.fq >> sample1_lane1_forwardreads.fq
ls -1 *_1.fq | while read line; do filename=$(echo $line | sed -E 's/(O\w\w\w).*/\1_R1.fq/'); cat $line >> $filename ; done
ls -1 *_2.fq | while read line; do filename=$(echo $line | sed -E 's/(O\w\w\w).*/\1_R2.fq/'); cat $line >> $filename ; done

# move to mapping folder
mv *_R[12]*.fq ../Sarah_scripts/final_cleaned_contigs/Mapping_collapsed/
# make a samples.txt list
ls -1 OC*_R1.fq | while read line; do echo $line | sed -E 's/(O\w\w\w).*/\1/' >> samples.txt; done

### Step 3 Trim adapaters and quality filter
# Going to use the same parameters from the transcriptome assembly
# this is the trim_reads.sh job script

#!/bin/bash
#SBATCH --partition=compute         # Queue selection
#SBATCH --job-name=trim_reads       # Job name
#SBATCH --mail-type=END             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hrivera@whoi.edu  # Where to send mail
#SBATCH --mem=32G                   # Job memory request
#SBATCH --time=02:00:00             # Time limit hrs:min:sec
#SBATCH --output=trim_reads_collapsed_%j.log  # Standard output/error
#SBATCH --cpus-per-task 12
#SBATCH --array=1-12

module load bio
module load anaconda/5.1
source activate trimmomatic

cd /vortexfs1/omics/tarrant/hrivera/DaviesData/Sarah_scripts/final_cleaned_contigs/Mapping_collapsed/
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)

trimmomatic PE -phred33 -threads 12 ./${LINE}_R1.fq ./${LINE}_R2.fq ./${LINE}_R1.P.qtrim ./${LINE}_R1.U.qtrim ./${LINE}_R2.P.qtrim ./${LINE}_R2.U.qtrim ILLUMINACLIP:/vortexfs1/omics/tarrant/hrivera/DaviesData/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25

### Step 4 Index full transcriptome
bowtie2-build combined_transcriptome.fasta combined_transcriptome

### Step 5 Map Reads (map_reads.sh)
#!/bin/bash
#SBATCH --partition=compute         # Queue selection
#SBATCH --job-name=map_reads       # Job name
#SBATCH --mail-type=END             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hrivera@whoi.edu  # Where to send mail
#SBATCH --mem=32G                   # Job memory request
#SBATCH --time=04:00:00             # Time limit hrs:min:sec
#SBATCH --output=maps_reads_collapsed_%j.log  # Standard output/error
#SBATCH --cpus-per-task 12
#SBATCH --array=1-12

module load bio
module load bowtie2

cd /vortexfs1/omics/tarrant/hrivera/DaviesData/Sarah_scripts/final_cleaned_contigs/Mapping_collapsed/
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt) # this prints the 1, 2nd up to array max line in the file samples.txt and assigns $LINE with that value

#bowtie2 --local -x ./combined_transcriptome -1 ./${LINE}_R1.P.qtrim -2 ./${LINE}_R2.P.qtrim -S ${LINE}.sam -q --no-hd --no-sq --no-unal -k 5 -p 12 >> ${LINE}_map
#decided to go with:
bowtie2 --local -x ./combined_transcriptome -1 ./trimmed_reads/${LINE}_R1.P.qtrim -2 ./trimmed_reads/${LINE}_R2.P.qtrim -S ${LINE}_best_alignment.sam -q --no-hd --no-sq --no-unal -p 12
# -no-hd suppress sam @ header lines
# -no-sq suppres sam @SQ header lines
# -no-unal don't print reads that didn't align
# -p designate number of threads for execution


### Step 6 Get counts (get_counts.sh)
#!/bin/bash
#SBATCH --partition=compute         # Queue selection
#SBATCH --job-name=get_counts      # Job name
#SBATCH --mail-type=END             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hrivera@whoi.edu  # Where to send mail
#SBATCH --mem=32G                   # Job memory request
#SBATCH --time=02:00:00             # Time limit hrs:min:sec
#SBATCH --output=get_sam_counts_%j.log  # Standard output/error
#SBATCH --cpus-per-task 1
#SBATCH --array=1-12

cd /vortexfs1/omics/tarrant/hrivera/DaviesData/Sarah_scripts/final_cleaned_contigs/Mapping_collapsed
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)

/vortexfs1/omics/tarrant/hrivera/DaviesData/Sarah_scripts/samcount.pl ${LINE}_best_alignment.sam combined_seq2iso.tab aligner=bowtie2 > ${LINE}_best_alignment.sam.counts

### Step 6 Clean counts (clean_counts.sh)
#!/bin/bash
#SBATCH --partition=compute         # Queue selection
#SBATCH --job-name=clean_counts      # Job name
#SBATCH --mail-type=END             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hrivera@whoi.edu  # Where to send mail
#SBATCH --mem=32G                   # Job memory request
#SBATCH --time=01:00:00             # Time limit hrs:min:sec
#SBATCH --output=clean_sam_counts_%j.log  # Standard output/error
#SBATCH --cpus-per-task 1
#SBATCH --array=1-12

cd /vortexfs1/omics/tarrant/hrivera/DaviesData/Sarah_scripts/final_cleaned_contigs/Mapping_collapsed/best_alignments
SAMP=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)

# Get only coral counts
cat ${SAMP}.sam.counts | grep 'Oarb' >> ${SAMP}_coral.sam.counts
# Get only sym counts
cat ${SAMP}.sam.counts | grep 'Sym' >> ${SAMP}_sym.sam.counts

# From here coral sam counts are good to go.
# The sym counts ended up not being usable all the attempts are described in alternate_pipeline.sh

# doing this locally:
ls -1 *alignment* | while read line; do SAMP=$(echo $line | grep -o -E  'OC[A-F][AS]'); cat $line | grep 'Oarb' >>${SAMP}_coral_sam.counts; cat $line | grep 'Sym' >> ${SAMP}_sym.sam.counts; done
# to get number of maps to coral and sym:
ls -1 *_coral* | while read line; do awk '{sum+=$2;} END {print sum}' $line; done
