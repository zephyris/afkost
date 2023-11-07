# AFKOST
**A**lignment-**f**ree **K*mer **o**mposition **S**earch **T**ool

## Getting started
Install the `afkost` package using pip from github. This requires `git` to be installed.
```shell
pip install git+https://github.com/zephyris/afkost
```

`afkost` has python modules dependencies: `numpy`, `matplot`. These are automatically installed when using `pip`.

To reinstall and upgrade use `pip`:
```shell
pip install --upgrade --force-reinstall git+https://github.com/zephyris/afkost
```

To uninstall also use `pip`:

```shell
pip uninstall afkost
```

## Kmer
### Initialise
Initialise a `Sequence` instance with the sequence, sequence type, and entropy filtering.
```python
sequence = Sequence("MDRGGRGGGSGGGRGSGAGGRGRGGAKMLARIAVISERY", kmer_max_length=2)
print(sequence)
```
By default, `Sequence` assumes a protein sequence. This can be set using `sequence_type` to `"protein"`, `"dna"` or `"rna"`, or three reduced complexity protein alphabets: `protein_reduced1`, `protein_reduced2` and `protein_reduced3`. DNA and RNA are largely untested.

Length of kmers to analyse are selected based on the number of possible amino acids/bases in the `sequence_type`, by default up to `3` for proteins. This can be set using `kmer_max_length`. Similarly maximum kgap length can be set with `kgap_max_length`.
Shannon entropy is calculated with a rolling window and the sequence may be trimmed to only residues above this threshold. By default, no thresholding is used as it removes some local sequence structure information. `entropy_threshold` and `entropy_window` can be set.

### Kmer-based sequence composition
The sequence is analysed by preserving local structure and discarding long range structure through two similar approaches.

First, by analysing in short kmers, ie. subsequences of length `k`.
For `k = 1` (and a protein sequence), this is equivalent to a count of number of each amino acid.
For `k = 2`, it is all amino acid pairs, ie. GG, GA, GL, GI, ... AG, AA, AL, AI, ... etc.
And so on, for higher values of `k`.
Counts of occurrence of each kmer can be viewed in the `Sequence` instance `kmers` attribute:
```python
print(sequence.kmers)
```

Second, by analysing `k`-spaced amino acid pairs. These are handled as the first and last position of a kmer.
For `k = 2`, this is equivalent to a kmer length two.
For `k = 3`, (and a protein sequence), this is all amino acid pairs with any amino acid between, ie. GXG, GXA, GXL, ... etc.
And so on, for higher values of `k`.
Counts of occurrence of each kgap can be viewed in the `Sequence` instance `kgaps` attribute:
```python
print(sequence.kgaps)
```

To handle different sequence lengths, kmer counts are normalised to the number of kmers of length `k` within the sequence. Number of kmers is `len(sequence) - k + 1`. The sum of normalised kmer frequency is therefore `1`, for each `k`. For sequences of length `0`, normalised kmer frequencies are all set to `0`.
This can be viewed in the `Sequence` `kmers_normalised` and `kgaps_normalised` attribute.
```python
print(sequence.kmers_normalised)
print(sequence.kgaps_normalised)
```

### Kmer comparison
Statistical tests for composition can be run:
```python
print(sequence.kmer_outlier_stats())
print(sequence.composition_dissimilarity_stats())
```
The former tests for kmers with significantly higher or lower frequency than expected, the latter tests overall composition.
By default, it assumes a naive random sequence where all kmers for each kmer length are expected to be equally frequent.

Of course, protein (and DNA/RNA) sequences aren't naive random sequences, and have preferred amino acid/bases and local sequence structure.
Average kmer composition can be loaded from a fasta file, which can use built-in sequence retrieval tools.
```python
from afkost import KmerMatrix
from afkost.databases import TriTrypDB
tritrypdb = TriTrypDB()
tritrypdb.fetch_fasta("TbruceiTREU927")
matrix = KmerMatrix()
matrix = matrix.composition_from_fasta(tritrypdb.fasta_path("TbruceiTREU927"))
```

Protein composition can then be compared to this matrix:
```python
print(sequence.kmer_outlier_stats(matrix))
```

This example uses a proteome-wide average composition. More intelligently generated input fasta files can give more useful comparison compositions.
