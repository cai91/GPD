The Gut Phage Database (GPD)
============================

Scripts used for characterizing human gut bacteriophages in the following manuscript:

Camarillo-Guerrero LF, Almeida A, Rangel-Pineros G, Finn RD, Lawley TD (2020) [Massive expansion of human gut bacteriophage diversity]

Associated data can also be found in our [FTP server](http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/)

<b>Classifier/</b>

## Neural network that distinguishes phages from integrative and conjugative elements (ICEs)

<b>Requirements:</b>

* Python (tested v3.6.7)
* TensorFlow (tested v1.10)
* Keras (tested v.2.2.4)

<b>Usage:</b> 
```
classifier.py <input_features_file.txt>
```

<b>Notes:</b>

input_features_file.txt: It contains a feature vector of 1026 dimensions: <fraction of hypothetical proteins> (1) <gene density> (1) <5-kmer signature> (1024) that represents a phage or an ICE (1 vector per line). 

## Other analysis and plotting scripts

<b>Figures/</b>
* 'Figure 1.py': Distribution of MIUViG scores from CheckV analysis
* 'Figure 2.py': Viral diversity patterns across gut bacteria genera and broad host range VCs
* 'Figure 3.py': Gut phageome profiling across human populations and correlation with gut bacteria enterotypes
* 'Figure 4.py': Crass-like family global distribution and host-phage network of globally distributed VCs
* 'Figure 5.py': Phylogenetic structure of the pX phage and global distribution
* 'Figure S1.py': Quality control assessment of GPD
* 'Figure S2.py': Viral diversity patterns across gut bacteria phyla and host range analysis of gut phages
* 'Figure S3.py': Correlation between sequencing depth and number of phages detected in a sample
* 'Figure S4.py': Host range analysis of globally distributed phages

