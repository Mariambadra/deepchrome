import subprocess
import os

# File paths
processed_dir = '../data/processed'
bins_bed = f'{processed_dir}/tss_bins.bed'

# The five histone marks in the exact order the paper uses
marks = ['H3K4me3', 'H3K4me1', 'H3K36me3', 'H3K9me3', 'H3K27me3']

print("Counting reads per bin for each histone mark...")
print("This will take a while — please be patient!\n")

for mark in marks:
    print(f"Processing {mark}...")
    
    bam_file    = f'{processed_dir}/E047-{mark}.sorted.bam'
    output_file = f'{processed_dir}/E047-{mark}_counts.txt'
    
    # bedtools multicov counts reads from BAM file
    # that overlap with each bin in our BED file
    cmd = f"bedtools multicov -bams {bam_file} -bed {bins_bed} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Quick sanity check
    with open(output_file) as f:
        line_count = sum(1 for line in f)
    
    print(f"  Done! {line_count:,} bins counted")
    print(f"  Saved to: {output_file}\n")

print("All 5 marks counted successfully!")
