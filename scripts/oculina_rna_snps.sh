# Directory layout
# Demultiplexed and trimmed reads by sample live in directory called trimmed_reads
# Contains 12 fastq files (named OC[A-F][AS].trim where the
# first letter variant is the coral colony ID and the second is whether it was an aposymbiotic (A) or symbiotic (S) branch)

#######################################
# To call SNPs we'll combine the apo and sym samples for each coral colony

cat OCAS.trim OCAA.trim >> OcA.trim
# repeated for each genet

### Re-map the concatenated reads to the just the host transcriptome, adding flags to make sure bowtie outputs the header lines and adding the read group info:
# SNP_Calling directory
# run as a slurm job script:
#####
#!/bin/bash
#SBATCH --partition=compute         # Queue selection
#SBATCH --job-name=merge-bam       # Job name
#SBATCH --mail-type=END             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hrivera@whoi.edu  # Where to send mail
#SBATCH --ntasks=1                  # Run X CPU
#SBATCH -N 1
#SBATCH --mem=28G                   # Job memory request
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=merge_bam_%j.log  # Standard output/error

cd /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_Calling/
module load bio
module load bowtie2

bowtie2 --local -p 20 -x /vortexfs1/scratch/hrivera/BU/Sarah_scripts/mapping/host_transcriptome_for_mapping/host_transcriptome -U /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcA.trim -S /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcA_forsnps.bam --no-unal --rg-id ColA -k 5 > /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/ColA
bowtie2 --local -p 20 -x /vortexfs1/scratch/hrivera/BU/Sarah_scripts/mapping/host_transcriptome_for_mapping/host_transcriptome -U /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcB.trim -S /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcB_forsnps.bam --no-unal --rg-id ColB -k 5 > /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/ColB
bowtie2 --local -p 20 -x /vortexfs1/scratch/hrivera/BU/Sarah_scripts/mapping/host_transcriptome_for_mapping/host_transcriptome -U /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcC.trim -S /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcC_forsnps.bam --no-unal --rg-id ColC -k 5 > /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/ColC
bowtie2 --local -p 20 -x /vortexfs1/scratch/hrivera/BU/Sarah_scripts/mapping/host_transcriptome_for_mapping/host_transcriptome -U /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcD.trim -S /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcD_forsnps.bam --no-unal --rg-id ColD -k 5 > /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/ColD
bowtie2 --local -p 20 -x /vortexfs1/scratch/hrivera/BU/Sarah_scripts/mapping/host_transcriptome_for_mapping/host_transcriptome -U /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcE.trim -S /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcE_forsnps.bam --no-unal --rg-id ColE -k 5 > /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/ColE
bowtie2 --local -p 20 -x /vortexfs1/scratch/hrivera/BU/Sarah_scripts/mapping/host_transcriptome_for_mapping/host_transcriptome -U /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcF.trim -S /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/OcF_forsnps.bam --no-unal --rg-id ColF -k 5 > /vortexfs1/omics/tarrant/hrivera/DaviesData/concat_reads/trimmed_reads/SNP_calling/ColF
#####

# Create a merged bam file with proper RG headers for GATK downstream analyses
samtools merge -rh rg.txt merged.bam *.sam

# index the input file
samtools sort merged.bam -o merged_sorted.bam -l 9 -@ 12
samtools index merged_sorted.bam -@ 12


# move the host transcriptome to the folder and make dictionary and index
# did this in interactive mode

gatk CreateSequenceDictionary -R host_newheader_transcriptome.fasta -O host_newheader_transcriptome.dict
samtools faidx host_newheader_transcriptome.fasta

## merge the reads
samtools mpileup -g -f host_newheader_transcriptome.fasta  merged_sorted.bam > merged_raw.bcf

# I still had to futz around with the header until I got it to work
# I used samtools view -H merged_sorted.bam > headers.txt
# to get the header info, then I had to add in proper header for each @RG tag
# I had also relabled all the files to be ColA.sam etc and remerged but there were still issues so I had to do this.
# Essentially using sed and piping I changed all the headers to
# @RG ID: ColA SN: ColA PL: Illumina for all the colonies with:

sed -E 's/ID:ColA/ID: ColA SN: ColA PL: Illumina/' headers.txt | sed -E 's/ID:ColB/ID: ColB SN: ColB PL: Illumina/' | other samples > newheader.txt
#then I used the samtools reheader tools

samtools reheader newheader.txt merged2_sorted.bam > merged2_newheader_sorted.bam


#!/bin/bash
#SBATCH --partition=compute         # Queue selection
#SBATCH --job-name=gatkhaplo       # Job name
#SBATCH --mail-type=END             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hrivera@whoi.edu  # Where to send mail
#SBATCH --ntasks=1                  # Run X CPU
#SBATCH -N 1
#SBATCH --mem=56gb                   # Job memory request
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=gatkhaplo_%j.log  # Standard output/error

module load bio
module load gatk/4.0.4.0

gatk HaplotypeCaller
     -R /vortexfs1/scratch/hrivera/BU/Sarah_scripts/mapping/host_transcriptome_for_mapping/host_newheader_transcriptome.fasta
     -I  merged2_newheader_sorted.bam
     -O  merged2_newheader_sorted.snps.indels.vcf


module load ddocent

# leave only SNPs where all samples have info
vcftools --vcf merged2_newheader_sorted.snps.indels.vcf --max-missing 1 --recode --recode-INFO-all --out no_missing_snps
# kept 34091 out of 107387 sites

# of loci with minimum quality of 30
vcftools --vcf no_missing_snps.recode.vcf --minQ 30 --recode --recode-INFO-all --out no_missing_snps_Q30
# kept 32327 out of a possible 34091 Sites

# loci with min depth of 10, didn't remove any data so didn't do it

# Convert variant calls to snps
vcfallelicprimitives no_missing_snps_Q30.recode.vcf --keep-info --keep-geno > no_missing_snps_Q30_genos.vcf

# Remove indels
vcftools --vcf no_missing_snps_Q30_genos.vcf --remove-indels --recode --recode-INFO-all --out Filtered_snps
# kept 30173 out of a possible 32327 Sites

# Only biallelic SNPs
vcftools --vcf Filtered_snps.recode.vcf --min-alleles 2 --max-alleles 2 --out Oculina_snps --recode --recode-INFO-all
# After filtering, kept 29663 out of a possible 30173 Sites
