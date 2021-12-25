import sys
import traceback

class Genome:
    def __init__(self, fasta_fn):
        self.genome = {}
        with open(fasta_fn) as fp:
            chrom = fp.readline().strip()[1:].split(" ")[0]
            seq = ''
            for line in fp:
                if line[0] == '>':
                    self.genome[chrom] = seq
                    chrom = line.strip()[1:].split(" ")[0]
                    seq = ''
                else:
                    seq += line.strip().upper()


    def get_seq(self, chrom, start, end, strand = '+'):
        seq = self.genome[chrom][start:end]
        print(seq)
        if strand == '-':
            seq = ''.join(map(self.revcomp, seq[::-1]))
        return seq



    def get_snv_seq(self, chrom, pos, cassette, length, ref = False):
        assert type(pos) == int and type(length) == int
        center = pos + len(cassette) / 2
        start = center - length / 2
        seq = self.genome[chrom][start: start + length]
        if len(seq) != length:
            return ''

        cassette_start = length / 2 - len(cassette) / 2
        if ref or not cassette:
            assert seq[cassette_start:cassette_start+len(cassette)] == cassette, \
            "Cassette is {}, but genome is {}".format(cassette,
                                                      seq[cassette_start:cassette_start+len(cassette)])
            return seq
        if cassette == '?':
            cassette = self.revcomp(seq[cassette_start])

        return (seq[:cassette_start]
                + cassette
                + seq[cassette_start+len(cassette):])


    def coords2fasta(self, coords, output_fasta):
        with open(output_fasta, 'w') as out:
            for chrom, three, strand in coords:
                if strand == '+':
                    seq = self.get_seq(chrom, three-70, three, strand)
                else:
                    seq = self.get_seq(chrom, three+1, three+70+1, strand)
                out.write(">{}:{}:{}\n".format(chrom, three, strand))
                out.write(seq + '\n')


    @staticmethod
    def revcomp(char):
        char = char.upper()
        if   char == 'A': return 'T'
        elif char == 'C': return 'G'
        elif char == 'G': return 'C'
        elif char == 'T': return 'A'
        assert char == 'N', "{} is not a valid base".format(char)
        return 'N'

    @staticmethod
    def swapNs(seq):
        if 'N' in seq:
            return seq.replace('N', 'A')
        return seq

class GTF:
    def __init__(gtf):
        self.gtf = gtf

    def read_gtf(self):
        coords = set()
        with open(self.gtf) as fp:
            for line in fp:
                if line[0] == '#': continue
                chrom, feature, start, end, strand, info = self._parse2gtflines(line)
                if 'gene_type' not in info or info['gene_type'] != 'protein_coding': continue
                if 'transcript_type' not in info or info['transcript_type'] != 'protein_coding': continue

                if feature == 'transcript':
                    cur_transcript = info['transcript_id']
                    last_exon = None

                elif feature == 'exon':
                    assert info['transcript_id'] == cur_transcript
                    if last_exon:
                        if strand == '+':
                            pos, three = last_exon, start
                        else:
                            three, pos = end-1, last_exon
                        coords.add((chrom, three, strand))
                    if strand == '+':
                        last_exon = end
                    else:
                        last_exon = start
        return list(coords)


    def _parse2gtflines(self, line):
        chrom, source, feature, start, end, _, strand, _, info = line.strip().split('\t')
        info = {el.split(' ')[0].strip('"'): el.split(' ')[1].strip('"')
                for el in info.strip(';').split('; ')}
        start, end = int(start)-1, int(end) # bed coordinates
        return chrom, feature, start, end, strand, info

"""
def coords2fasta(coords, genome, oname):
    with open(oname, 'w') as out:
        for chrom, three, strand in coords:
            if strand == '+':
                seq = genome.get_seq(chrom, three-70, three, strand)
            else:
                seq = genome.get_seq(chrom, three+1, three+70+1, strand)
            out.write(">{}:{}:{}\n".format(chrom, three, strand))
            out.write(seq + '\n')
"""

def gtf2fasta(gtf, genome, oname):
    """
    gtf file to matching fasta

    Parameters
    ----------
    gtf: string
        gtf file path
    genome: string
        genome file path
        default value is hg19.fa
    oname: string
        fasta file

    Results
    -------
    log: dictionary
        exception or program result parameters

    Examples
    --------
    log = gtf2fasta(gtf, genome, oname)
    print(log)
    """
    try:
        print("Parsing gtf file at {}".format(gtf))
        coords = GTF(gtf).read_gtf()
        print(f'coods : {coords}')
        print(coords[0])
        print("Found {} 3'ss.".format(len(coords)))
        print("Loading genome fasta file from {}".format(genome))
        genome = Genome(genome)
        print(f"genome : {genome}")
        print("Writing output to {}".format(oname))
        genome.coords2fasta(coords, oname)
        print('Success!')

    except ValueError as ve:
        print(traceback.format_exc())
        return {"conditions" : "error", "genome" : genome, "coords" : coords,  "log" : ve }

    except Exception as e:
        print(traceback.format_exc())
        return {"conditions" : "error", "genome" : genome, "coords" : coords, "log" : e }
    return { "conditions" : "success", "genome" : genome, "coords" : coords}


if __name__ == '__main__':
    try:
        genome, gtf, oname = sys.argv[1:]
    except ValueError as ve:
        print(ve)
        exit()
    gtf2fasta(gtf, genome, oname)

