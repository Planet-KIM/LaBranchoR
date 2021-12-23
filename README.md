# LaBrachoR (LSTM Branchpoint Retriever)
**LaBranchoR uses a LSTM network built with keras to predict the
position of RNA splicing branchpoints relative to a three prime
splice site. Precisely evaluating LaBranchoR was challenging due
to pervasive noise in the experimental data, but as we show in our
paper, we estimate that LaBranchoR correcty predicts a branchpoint
for over 90% of 3'ss.**


## 1. Download existing branchpoint annotations
**See our website linked above to download branchpoint predictions
for introns in gencode v19 (hg19) or view LaBranchoR predicted
branchpoints in the UCSC genome browser.**


## 2. Running LaBranchoR
**All of the code and model weights needed to run LaBranchoR are available in
the 'labranchor' directory. Running LaBranchoR requires keras and numpy to be installed.**


### 2-1). Predicting branchpoints
#### ■ Terminal
The script __init__.py makes predictions for a fasta file of sequences upstream of 3'ss.
```sh
python labranchor_v2/__init__.py weights 'top-bed'/'top'/'all' fasta_file output
```

### 2-2). Description of Parameters (labrachor_v2/__init__.py)
#### weights
- The path to the h5 weights files
  - labranchor/2layer.h5 (original model)
  - labranchor/2layer2.h5 (batch model)

#### 'top-bed'/'top'/'all'
- top-bed:
  - produces a bed file of predicted branchpoints. 
  - Assumes fasta names are chrom:3'ss_coord:strand (ex. chr1:1000:+)
- top:
  - reports the shift of the top scoring branchpoint from the associated 3'ssfor each fasta entry
- all:
  - reports a comma seperated list of branchpoint probabilities corresponding to positions -70 to -1 from each 3'ss

#### fasta_file
- Path to a fasta file of sequences upstream of 3'ss. Input sequences are required to be 70 base pairs and should not contain characters other than 'A', 'C', 'G', 'T', or 'N'. Any Ns will be considered A's during prediction.

#### output
- Path to the output file. See the above options for formatting.


## 3. Creating 3'ss sequence fasta files
**The script genome.py can be used to create fasta files suitable for branchpoint prediction for all introns in given gtf file.**

### 3-1). Running genome
#### ■ Terminal
```sh
python labrachor_v2/genome.py <genome> <gtf> <output>
```

### 3-2). Description of Parameters(genome.py)
#### genome
- A path to a genome fasta file consistent with the gtf file.

#### gtf
- The path to the gtf file you wish to predict branchpoints in.

#### output
- The path to the output fasta file.


## 4. Analysis Included in Paper

Model training: notebooks/train_model.ipynb

Model performance: notebooks/performance_*

Cases where LaBranchoR disagrees with experimental data: notebooks/disagreement_*

Genome-wide properties and overlap with pathogenic variants: notebooks/landscape_*

Properties of C and no -2 U branchpoints: notebooks/landscape_C_and_noT.ipynb

Enrichments of ExAC variants: notebooks/ExAC_variant_enrichments.ipynb

Generation of ISM supplmentary data: notebooks/supp_data.ipynb

## 5. Analysis not included in paper

Exploration of nucleotide importances: notebooks/importance.ipynb

Analysis of secondary structure near branchpoints: notebooks/secondary_*

## 6. References
Paggi J.M., Bejerano, G. A sequence-based, deep learning model accurately predicts RNA splicing branchpoints. bioRxiv 185868 (2017). DOI:[10.1101/185868](http://www.biorxiv.org/content/early/2017/09/07/185868)

## 7. Contacts
If having to run the model yourself would stop you from using LaBranchoR,
please open an issue requesting the desired predictions or contact the
authors via email.
