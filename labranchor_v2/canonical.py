import json
import sys
import re
import traceback
from jklib.genome import locus
from jklib.bioDB.canonical_transcript import CanonicalTranscriptSelector


def canonical_variant(variant):
    try:
        variant = variant.split(' ')
        changed = variant[1].split('>')
        ref = changed[0]
        alt = changed[1]

        if type(variant[0])==str:
            variant = locus(variant[0])
       
        if variant.twoBitFrag().upper() != ref:
            print(f'{variant.twoBitFrag().upper()} -> {ref}')
            raise TypeError({'msg' : 'Posision Error'})

        regions = variant.regionType()
        check = None
        for region in regions:
            for item in region[3]:
                if check == 'canonical':
                    break
                if 'intron' in item:
                    sta, end = item[1]
                    if variant.strand == '+':
                        if end == -1:
                            check = 'canonical'
                    else:
                        if sta == -1:
                            check = 'canonical'
                else:
                    continue
        if check != 'canonical':
            raise Exception('This variant is not canonical')
        else:
            print(sta,end)
            if variant.strand == "+":
                intron_70base = f'{variant.chrom}:{variant.chrSta+end-69}-{variant.chrSta+end}{variant.strand}'
            else:
                intron_70base = f'{variant.chrom}:{variant.chrSta-sta}-{variant.chrSta-sta+69}{variant.strand}'
            print(intron_70base, locus(intron_70base).twoBitFrag().upper())
            return intron_70base, locus(intron_70base).twoBitFrag().upper()

    except TypeError as te:
        print(traceback.format_exc())
        return {"Error" : te }
    
    except Exception as e:
        print(traceback.format_exc())
        return {"Error" : e } 


def write2fasta(intron_70base):
    variant, fasta = intron_70base
    variant = locus(variant)
    with open('./variant.fa', 'w') as fasta_file:
        fasta_file.write(f'>{variant.chrom}:{variant.chrSta}:{variant.strand}\n{fasta}')
    print('success to write fasta file')


if __name__ == '__main__':
    try:
        variant = sys.argv[1:]
    except ValueError as ve:
        print(ve)
        exit()
    intron_70base = canonical_variant(variant[0])
    variant2fasta(intron_70base)
    print('success')

#variant="chr11:108389037-108389037 G>A"
#variant="chrX:15581390-15581390- G>T"
#variant="chr11:108365082-108365082+ G>T"
#intron_70base = canonical_variant(variant)
#variant2fasta(intron_70base)



