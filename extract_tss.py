import gzip
import csv

# Input and output files
gtf_file = '../data/raw/gencode.v19.annotation.gtf.gz'
output_file = '../data/processed/tss_coordinates.bed'

# Store TSS coordinates
tss_list = []

print("Extracting TSS coordinates from GTF file...")

with gzip.open(gtf_file, 'rt') as f:
    for line in f:
        
        # Skip header lines
        if line.startswith('#'):
            continue
        
        # Split line into fields
        fields = line.strip().split('\t')
        
        # We only want gene-level entries
        if fields[2] != 'gene':
            continue
        
        # Extract attributes
        attributes = fields[8]
        
        # We only want protein-coding genes
        if 'gene_type "protein_coding"' not in attributes:
            continue
        
        # Extract basic coordinates
        chrom  = fields[0]
        start  = int(fields[3])
        end    = int(fields[4])
        strand = fields[6]
        
        # Extract gene_id
        for attr in attributes.split(';'):
            if 'gene_id' in attr:
                gene_id = attr.strip().split('"')[1]
                break
        
        # TSS depends on strand direction
        # For + strand: TSS is at the start
        # For - strand: TSS is at the end
        if strand == '+':
            tss = start
        else:
            tss = end
        
        # Store as BED format
        # (chromosome, tss-5000, tss+5000, gene_id, score, strand)
        # Calculate window coordinates
        window_start = tss - 5000
        window_end = tss + 5000

        # Skip genes too close to chromosome start
        if window_start < 0:
            continue

        # Store as BED format
        tss_list.append([
            chrom,
            window_start,
            window_end,
            gene_id,
            '.',
            strand
        ])
# Write to BED file
with open(output_file, 'w') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerows(tss_list)

print(f"Done! Extracted {len(tss_list)} protein-coding genes")
print(f"Saved to: {output_file}")
