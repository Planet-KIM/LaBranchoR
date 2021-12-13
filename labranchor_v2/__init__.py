__version__ = 'v0.0.2' 

import sys
import numpy as np
from tensorflow.keras.models import load_model

class labranchor:
    def __init__(self, weights, top, fasta, output):
        try:
            self.L = 70
            if (weights == None or top == None or fasta == None or output or None):
                raise Exception("Parameter Error")
            if top not in ['all', 'top', 'top-bed']:
                raise Exception("2nd argument must be 'all', 'top', or 'top-bed'!")
            self.weights = weights
            self.top = top
            self.fasta = fasta
            self.output = output
        except Exception as e:
            self.usage()


    def usage(self):
        print("Usage: python labrachor.py <weights> <'top-bed'/'top'/'all'> <fasta_file> <output>")
        print()
        print("Weights file should likely be path to '2layer.h5'")
        print("'top' returns the top scoring branchpoint per 3'ss, while 'all' gives 70 probabilities for positions -70 through -1.")
        print("The fasta file should contain length 70 sequences right-aligned to a 3'ss. I.e. most should end in 'AG'.")
        print("The output must be a file name, pipes are not supported.")
        exit()


    def onehot(self, seq):
        bases = ['A', 'C', 'G', 'T']
        X = np.zeros((len(seq), len(bases)))
        for i, char in enumerate(seq):
            X[i, bases.index(char)] = 1
        return X


    def encode(self, seqs):
        return np.vstack(self.onehot(seq).reshape(1, self.L, 4) for seq in seqs)


    def read_fasta(self):
        names, seqs = [], []
        with open(self.fasta) as fp:
            first = True
            for line in fp:
                line = line.strip()
                if line[0] == '>':
                    names += [line[1:]]
                    if not first and len(seqs[-1]) != self.L:
                        print('All sequences must have length 70.')
                        usage()
                    first = False
                    seqs += ['']
                else:
                    line = line.upper()
                    if 'N' in line: line = line.replace('N', 'A')
                    seqs[-1] += line
        assert len(names) == len(seqs)
        return names, seqs


    def write_output(self, names, preds):
        with open(self.output, 'w') as out:
            if self.top == 'top-bed':
                for name, pos, score in zip(names, np.argmax(preds, axis = 1), np.max(preds, axis = 1)):
                    chrom, three, strand = name.split(':')
                    three = int(three)
                    bp = three + (pos-self.L) if strand == '+' else three - (pos-self.L)
                    out.write('\t'.join(map(str, [chrom, bp, bp+1, three, score, strand])) + '\n')
            elif self.top == 'top':
                for name, pos, score in zip(names, np.argmax(preds, axis = 1), np.max(preds, axis = 1)):
                    out.write("{}\t{}\t{}\n".format(name, pos-self.L, score))
            else:
                for name, pred in zip(names, preds):
                    out.write("{}\t{}\n".format(name, ','.join(map(str, pred))))
        print("Success to write the prediction file")


    def predict(self):
        print(f"Loading sequences from {self.fasta}")
        names, seqs = self.read_fasta(self.fasta)
        print(f"Read {len(names)} sequences")
        
        # one hot encoding
        print('Encoding sequences.')
        X = self.encode(seqs)
        
        # load the tensorflow model
        # predict the branch points according to LaBranchoR model
        print(f"Loading model from {self.weights}")
        model = load_model(self.weights)
        print('Making predictions.')
        preds = model.predict(X).reshape(-1, self.L)
        
        print(f"Writing predictions to {self.output}")
        self.write_output(names, preds)
        
        print('Success to predict!')


if __name__ == '__main__':
    try:
        weights, top, fasta, output = sys.argv[1:]
    except ValueError as ve:
        print(ve)
        exit()
    
    print(weights, top, fasta, output)
    
    labranchor_run = labranchor(weights=weights, top=top, fasta=fasta, output=output)
    labranchor_run.predict()
