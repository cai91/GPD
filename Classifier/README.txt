Classifier from "Massive expansion of human gut bacteriophage diversity "
Author: Luis Fernando Camarillo Guerrero

Dependencies: TensorFlow v1.10, Keras v.2.2.4
Input: A vector of 1026 dimensions: <fraction of hypothetical proteins> (1) <gene density> (1) <5-kmer signature> (1024)
Output: Probability of a prediction being a phage (Only predictions with scores > 0.8 were considered for this manuscript)

Demo files:

KGDU_phages_features.txt: It contains the features of 50 phages.
KGDU_ICEs_features.txt: It contains the features of 50 ICEs.
