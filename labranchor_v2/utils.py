from labranchor_v2.genome import Genome


def load_datas(fasta='/home/dwkim/gencode_v19_raw_data/labranchor.gencode_v19.fa',
            intron_data='/home/dwkim/anno/introns_to_mercer.tsv'):
    genome = Genome(fasta)
    introns = {}
    with open(intron_data) as fp:
        for line in fp:
            chrom, start, end, _, pos, strand, _, bp = line.split('\t')[:8]
            bp, start, end = int(bp), int(start), int(end)
                                    
            three = end if strand == '+' else start  
            key = (chrom, three, strand)
                                                            
            if not 5 < abs(bp - three) < 60:
                bp = -1
                                                                                            
            if key not in introns: introns[key] = []
            assert bp not in introns[key], bp
            if bp != -1: introns[key] += [bp]
    
    known   = {key: value for key, value in introns.items() if value}
    missing = {key: value for key, value in introns.items() if not value}
    return known, missing

