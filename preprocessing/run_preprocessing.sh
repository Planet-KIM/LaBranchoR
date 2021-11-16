# Input files:
# Gencode v19 annotatations:     ../anno/gencode.v19.annotation.gtf
# Mercer et. al., 2015 bp anno:  ../anno/mercer.suppdataS2.bed
# Taggart et. al., 2017 bp anno: ../anno/taggart.supptable1.csv  (Exported from excel)
# hg19 Phastcons 100-way scores: ../anno/hg19.100way.phastcons.bw


## Join BP annotations to 3'ss.

# Create a bed file with an entry for each unique 3'ss present in gencode v19.
# Entries consist of:
# CHROM (THREE-100) THREE GENE_ID FIVE STRAND
# The coordinates are transposed if the 3'ss aligns to the minus strand.
#cat ../anno/gencode.v19.annotation.gtf | python bp_intron_bed.py | sort -u -k1,1 -k2,3n > gencode.v19.introns.bpregion.bed

# Put Tagart et al annotations into bed format
#cat ../anno/taggart.supptable1.csv | python taggart_to_bed.py | sort -k1,1 -k2,3n > ../anno/taggart.supptable1.bed

# Put Signal et al predictions into bed format
#tail -n +2 ../anno/gencode_v19_branchpoints.csv | awk -F , '{print $1"\t"$2"\t"$2+1"\t"$5"\t.\t"$4}' > ../anno/gencode_v19_branchpoints.bed

# Join bp annotations to introns.
#bedtools intersect -loj -s -a gencode.v19.introns.bpregion.bed -b ../anno/mercer.suppdataS2.bed > introns_to_mercer.tsv
#bedtools intersect -loj -s -a gencode.v19.introns.bpregion.bed -b ../anno/taggart.supptable1.bed > introns_to_taggart.tsv
#bedtools intersect -loj -s -a gencode.v19.introns.bpregion.bed -b ../anno/gencode_v19_branchpoints.bed > introns_to_signal.tsv


## Prepare conservation scores.
## Note that this will take longer than the above.
## Comment out if previously run

# Convert conservation scores to bedgraph.
./bigWigToBedGraph ../anno/hg19.100way.phastcons.bw ../anno/hg19.100way.phastcons.bedGraph
./bigWigToBedGraph ../anno/hg19.100way.phyloP100way.bw ../anno/hg19.100way.phyloP100way.bedGraph

# Extract only entries we are interested in.
bedtools intersect -u -a ../anno/hg19.100way.phyloP100way.bedGraph -b gencode.v19.introns.bpregion.bed > ../anno/hg19.100way.phyloP100way.bp.bed
bedtools intersect -u -a ../anno/hg19.100way.phastcons.bedGraph -b gencode.v19.introns.bpregion.bed > ../anno/hg19.100way.phastCons.bp.bed
#bedtools intersect -u -a ../../primates/hg19.phyloP46way.primate.bedGraph -b gencode.v19.introns.bpregion.bed > ../anno/hg19.46way.phylop.primates.bp.bed