import subprocess
import os

# Input and output files
tss_bed = '../data/processed/tss_coordinates.bed'
bins_bed = '../data/processed/tss_bins.bed'

print("Creating 100bp bins around TSS regions...")

# Use bedtools makewindows to split each TSS window into 100bp bins
# -b: input BED file
# -w: window size (100bp)
# -i: use the gene ID from column 4 as bin name
cmd = f"bedtools makewindows -b {tss_bed} -w 100 -i srcwinnum > {bins_bed}"

subprocess.run(cmd, shell=True, check=True)

# Count the result
with open(bins_bed) as f:
    total_bins = sum(1 for line in f)

print(f"Done! Created {total_bins} bins total")
print(f"Expected: 20,345 genes × 100 bins = {20345 * 100:,} bins")
print(f"Saved to: {bins_bed}")
