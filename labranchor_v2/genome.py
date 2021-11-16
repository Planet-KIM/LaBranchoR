
class Genome:
    def __init__(self, fasta_fn):
        self.genome = {}
        with open(fasta_fn) as fp:
            chrom = fp.readline().strip()[1:]
            seq = ''
            for line in fp:
                if line[0] == '>':
                    self.genome[chrom] = seq
                    chrom = line.strip()[1:]
                    seq = ''
                else:
                    seq += line.strip().upper()

    def get_seq(self, chrom, start, end, strand = '+'):
        #print(chrom, start, end, strand)
        key_hash = [ item for item in range(int(start)-1, int(end)+1)]
        seq =''
        for k, v in self.genome.items():
            variant = k.split(":")
            if int(variant[1]) in key_hash:
                if variant[0] == chrom:
                    if variant[2] == strand:
                        seq = v
        #seq = self.genome[chrom][start:end]
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
