# AFKOST
**A**lignment-**f**ree **K*mer **o**mposition **S**earch **T**ool


## Kmer
### Initialise
Initialise a `Sequence` instance with the sequence, sequence type, and entropy filtering.
```
sequence = Sequence("MDRGGRGGGSGGGRGSGAGGRGRGGAKMLARIAVISERY", kmer_max_length=2)
print(sequence)
```
By default, `Sequence` assumes a protein sequence. This can be set using `sequence_type` to `"protein"`, `"dna"` or `"rna"`. DNA and RNA are largely untested.
By default, kmers up to length `3` are analysed. This can be set using `kmer_max_length`.
Shannon entropy is calculated with a rolling window and the sequence may be trimmed to only residues above this threshold. By default, no thresholding is used as it removes some local sequence structure information. `entropy_threshold` and `entropy_window` can be set.

### Kmer-based sequence composition
The sequence is analysed by preserving local structure and discarding long range structure by analysing in short kmers, ie. subsequences of length `k`.
For `k = 1` (and a protein sequence), this is equivalent to a count of number of each amino acid.
For `k = 2`, it is all amino acid pairs, ie. GG, GA, GL, GI, ... AG, AA, AL, AI, ... etc.
And so on, for higher values of `k`.
Counts of occurrence of each kmer can be viewed in the `Sequence` `kmers` attribute:
```
print(sequence.kmers)
```

To handle different sequence lengths, kmer counts are normalised to the number of kmers of length `k` within the sequence. Number of kmers is `len(sequence) - k + 1`. The sum of normalised kmer frequency is therefore `1`, for each `k`. For sequences of length `0`, normalised kmer frequencies are all set to `0`.
This can be viewed in the `Sequence` `kmers_normalised` attribute.
```
print(sequence.kmers_normalised)
```

### Kmer comparison
Statistical tests for composition can be run:
```
print(sequence.kmer_outlier_stats())
print(sequence.composition_dissimilarity_stats())
```
The former tests for kmers with significantly higher or lower frequency than expected, the latter tests overall composition.
By default, it assumes a naive random sequence where all kmers for each kmer length are expected to be equally frequent.

Of course, protein (and DNA/RNA) sequences aren't naive random sequences, and have preferred amino acid/bases and local sequence structure.
Average kmer composition can be loaded from a fasta file, which can use built-in sequence retrieval tools.
```
from afkost import KmerMatrix, TriTrypDB
tritrypdb = TriTrypDB()
tritrypdb.fetch_fasta("TbruceiTREU927")
matrix = KmerMarix()
matrix = matrix.composition_from_fasta("_tritrypdb/TbruceiTREU927")
```

Protein composition can then be compared to this matrix:
```
print(sequence.kmer_outlier_stats(matrix))
```

This example uses a proteome-wide average composition. More intelligently generated input fasta files can give more useful comparison compositions.