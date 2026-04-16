import subprocess
import os

# File paths
raw_dir = '../data/raw'
processed_dir = '../data/processed'
chrom_sizes = f'{raw_dir}/hg19.chrom.sizes'

# The five histone marks
marks = ['H3K4me3', 'H3K4me1', 'H3K36me3', 'H3K9me3', 'H3K27me3']

for mark in marks:
    print(f"\nProcessing {mark}...")
    
    input_file  = f'{raw_dir}/E047-{mark}.tagAlign.gz'
    bam_file    = f'{processed_dir}/E047-{mark}.bam'
    sorted_bam  = f'{processed_dir}/E047-{mark}.sorted.bam'
    
    # Step 1: Convert tagAlign to BAM
    print(f"  Converting to BAM...")
    cmd1 = f"bedtools bedtobam -i {input_file} -g {chrom_sizes} > {bam_file}"
    subprocess.run(cmd1, shell=True, check=True)
    
    # Step 2: Sort the BAM file
    print(f"  Sorting BAM...")
    cmd2 = f"samtools sort {bam_file} -o {sorted_bam}"
    subprocess.run(cmd2, shell=True, check=True)
    
    # Step 3: Index the sorted BAM
    print(f"  Indexing BAM...")
    cmd3 = f"samtools index {sorted_bam}"
    subprocess.run(cmd3, shell=True, check=True)
    
    # Step 4: Remove unsorted BAM to save space
    print(f"  Cleaning up...")
    os.remove(bam_file)
    
    print(f"  Done! {mark} → {sorted_bam}")

print("\n All 5 marks converted successfully!")
