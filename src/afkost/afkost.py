# base modules
import os
import math
#import subprocess # currently unused
import time # currently unused
#import random # currently unused
import warnings
import hashlib
from functools import cached_property

# modules from pip
import numpy
from scipy import stats
#import requests # currently unused
#from tqdm.auto import tqdm # currently unused
#import matplotlib.pyplot as plt # slow import, so import only in `_Plot`

class Alphabet:
    def __init__(self, sequence_type: str = "protein"):
        """Initialise object containing alphabet information

        Keyword arguments:
        sequence_type -- "protein", "dna" or "rna" sequence type, defining expected alphabet (default `"protein"`)
        """    
        # setup sequence_type
        self.sequence_type = "".join(sequence_type.lower().split())
        self.alphabets = {
            "protein": { # 20 amino acids
                "alphabet": "RKHDESTNQAVILMFYWCGP",
                "alphabet_replace": {"U": "C", "O": "K"},
                "alphabet_any": "X",
                "typical_length": 1000
            }, "protein_reduced": { # simplified to 9 amino acid classes
                "alphabet": "RHDSAFCGP",
                "alphabet_replace": {
                    "K": "R", # KR -> R, basic, pKa ~11
                    "E": "D", # DE -> D, acidic pKa ~4
                    "T": "S", "N": "S", "Q": "S", # STNQ -> S, polar uncharged
                    "V": "A", "I": "A", "L": "A", "M": "A", # AVILM -> A, small hydrophobic
                    "Y": "F", "W": "F" # FYW -> F, hydrophobic aromatic
                },
                "alphabet_any": "X",
                "typical_length": 1000
            }, "protein_reduced2": { # simplified to 7 amino acid classes
                "alphabet": "RDSAFGP",
                "alphabet_replace": {
                    "K": "R", "H": "R", # KRH -> R, basic
                    "E": "D", # DE -> D, acidic
                    "T": "S", "N": "S", "Q": "S", # STNQ -> S, polar uncharged
                    "V": "A", "I": "A", "L": "A", "M": "A", "C": "A", # AVILMC -> A, small hydrophobic
                    "Y": "F", "W": "F" # FYW -> F, hydrophobic aromatic
                },
                "alphabet_any": "X",
                "typical_length": 1000
            }, "protein_reduced3": { # simplified to 4 amino acid classes
                "alphabet": "RDSA",
                "alphabet_replace": {
                    "K": "R", "H": "R", # KRH -> R, basic
                    "E": "D", # DE -> D, acidic
                    "T": "S", "N": "S", "Q": "S", "G": "S", # STNQG -> S, hydrophilic
                    "V": "A", "I": "A", "L": "A", "M": "A", "C": "A", "F": "A", "Y": "A", "W": "A", "P": "A" # AVILMCFYWP -> A, hydrophobic
                },
                "alphabet_any": "X",
                "typical_length": 1000
            }, "dna": {
                "alphabet": "ATCG",
                "alphabet_replace": {},
                "alphabet_any": "N",
                "typical_length": 3000
            }, "rna": {
                "alphabet": "AUCG".upper(),
                "alphabet_replace": {},
                "alphabet_any": "N",
                "typical_length": 3000
            }
        }
        if self.sequence_type not in self.alphabets:
            raise ValueError("`sequence_type` must be `\"protein\"` `\"dna\"` or `\"rna\"`")
        self.alphabet_any = self.alphabets[self.sequence_type]["alphabet_any"]
        self.alphabet = self.alphabets[self.sequence_type]["alphabet"] + self.alphabet_any
        self.alphabet_replace = self.alphabets[self.sequence_type]["alphabet_replace"]
        self.typical_length = self.alphabets[self.sequence_type]["typical_length"]

    def __repr__(self):
        """
        Represent alphabet as the alphabets entry
        """
        return str(self.alphabets[self.sequence_types])
    
    def __str__(self):
        """
        Human readable alphabet description
        """
        return self.sequence_type + " alphabet: " + self.alphabet + ", any/unknown: " + self.alphabet_any

    def all_kmers(self, k: int):
        """
        List of all possible kmers of length k

        Required arguments:
        k -- kmer length

        Returns:
        List of all kmers of length k as list of strings
        """
        all_kmers = []
        for i in range(0, len(self.alphabet) ** k):
                #Construct current kmer
                kmer = ""
                for j in range(0, k):
                    kmer = self.alphabet[math.floor(i / len(self.alphabet) ** j) % len(self.alphabet)] + kmer
                all_kmers.append(kmer)
        return all_kmers

    def num_kmers(self, k: int):
        """
        Number of unique possible kmers of length k
        
        Required arguments:
        k -- kmer length

        Returns
        Count of the number of possible kmers
        """
        return self.alphabet_length ** k
        
    def default_k(self, typical_length: int = None):
        """
        Automatic selection of `kmer_max_length` and `kgap_max_length` based on sequence typical length and alphabet complexity.

        Named arguments:
        typical_length -- typical sequence length for this dataset, default 1000 for protein, 3000 for nucleotide

        Returns tuple of default kmer and kgap max lengths.
        """
        if typical_length is None:
            typical_length = self.typical_length
        kmer_max_length = 0
        while typical_length > self.num_kmers(kmer_max_length):
            # increase kmer_max_length until number of possible kmers just exceeds typical sequence length
            # for 20 aas, typical protein length 1000aas, kmer_max_length = 3
            # for simplified to 10aa classes, typical protein length 1000aas, 
            # for 4 bases, typical transcript length 3000bps, kmer_Max_lenght = 6
            kmer_max_length += 1
        kgap_max_length = kmer_max_length * 2
        return kmer_max_length, kgap_max_length

    @cached_property
    def alphabet_lookup(self):
        """
        Lookup to convert sequence characters (kmers of length 1) to numerical indices
        """
        alphabet_lookup = {}
        for index, character in enumerate([x for x in self.alphabet]):
            alphabet_lookup[character] = index
        return alphabet_lookup
    
    @cached_property
    def alphabet_length(self):
        """
        Length, ie. number of unique characters (kmers of length 1) in the alphabet
        """
        return len(self.alphabet)

class Fasta:
    def __init__(self, path):
        """Initialise class containing sequence information from a fasta file path
        
        Required arguments:
        path -- the location of the fasta file to load
        """

        # record settings
        self.path = path
    
    @cached_property
    def sequences(self):
        """
        Dict of all sequences in the fasta file, in the form {"sequence_name": "sequence", ...}
        """
        with open(self.path, "r") as fasta_file:
            data = "\n"+fasta_file.read() # ensure prefixed with a newline
            data = data.replace('\r\n', '\n').replace('\r', '\n').split("\n>")[1:] # tidy newline, then split by `"\n>"``
            names = [x.splitlines()[0].split(" ")[0] for x in data] # name from up to the first space in the first line, per sequence
            seqs = ["".join(x.splitlines()[1:]) for x in data] # sequence from concatenated remaining lines, per sequence
            sequences = {}
            for i in range(len(names)):
                sequences[names[i]] = seqs[i]
            return sequences
    
    @cached_property
    def first_sequence(self):
        """
        Sequence of the first entry in the fasta file, as a string
        """
        #return self.sequences[next(iter(self.sequences))]
        return self.sequences[list(self.sequences.keys())[0]]

    @cached_property
    def first_name(self):
        """
        Name of the first entry in the fasta file, as a sting
        """
        return list(self.sequences.keys())[0]

    @cached_property
    def count(self):
        """
        Number of sequences in the fasta file
        """
        return len(self.sequences)

class Kmer:
    def __init__(self, kmer_max_length, kgap_max_length, sequence_type):
        """
        Initialise class for internal handling of kmers
            
        Required arguments:
        sequence_type -- "protein", "dna" or "rna" sequence type, defining expected alphabet (default `"protein"`)
        kmer_max_length -- analyse composition in kmers from 1 up to and equal to this length
        kgap_max_length -- analyse composition in kgaps from 3 up to and equal to this length
        """
        self.kmer_max_length = kmer_max_length
        self.kgap_max_length = kgap_max_length
        self.sequence_type = sequence_type
        self.alphabet = Alphabet(self.sequence_type)
        self.sequence = None
        # setup empty kmer object of correct size, using numpy to store the lists
        # assumes that alphabet.all_kmers gives a stable lookup for kmer sequence lookup
        self.kmer_arrays = {}
        self.kmer_normalised_arrays = {}
        for k in range(1, kmer_max_length + 1):
            self.kmer_arrays[k] = numpy.zeros(self.alphabet.num_kmers(k), dtype=float)
            self.kmer_normalised_arrays[k] = numpy.zeros(self.alphabet.num_kmers(k), dtype=float)
        # setup empty kgap object of correct size, using numpy to store the lists
        self.kgap_arrays = {}
        self.kgap_normalised_arrays = {}
        for k in range(3, kgap_max_length + 1):
            self.kgap_arrays[k] = numpy.zeros(self.alphabet.num_kmers(2), dtype=float)
            self.kgap_normalised_arrays[k] = numpy.zeros(self.alphabet.num_kmers(2), dtype=float)
    
    def __str__(self):
        """
        Human readable summary of most common kmers
        """
        # append common kmers
        indicence_min = 0.005
        string = "KMERs with a frequency of at least " + str(indicence_min * 100) + "%:\n"
        for k in range(1, self.kmer_max_length + 1):
            all_kmers = self.alphabet.all_kmers(k)
            string += "  KMER length " + str(k) + ":\n    "
            # list abundant kmers
            count = 0
            for i in range(len(all_kmers)):
                if self.kmer_normalised_arrays[k][i] > indicence_min:
                    string += all_kmers[i] + ":" + str(self.kmer_arrays[k][i]) + " "
                    count += 1
            if count == 0:
                string += "[none]"
            string += "\n"
        string += "KGAPs with a frequency of at least " + str(indicence_min * 100) + "%:\n"
        for k in range(3, self.kgap_max_length + 1):
            all_kmers = self.alphabet.all_kmers(2)
            string += "  KGAP length " + str(k) + ":\n    "
            # list abundant kgaps
            count = 0
            for i in range(len(all_kmers)):
                if self.kgap_normalised_arrays[k][i] > indicence_min:
                    string += all_kmers[i][0] + "-" + all_kmers[i][1] + ":" + str(self.kgap_arrays[k][i]) + " "
                    count += 1
            if count == 0:
                string += "[none]"
            string += "\n"
        return string[:-1]
    
    def _type_match_check(self, kmer):
        """
        Check that `kmer` is a `Kmer` instance with matching `sequence_type` and `max_kmer_length` to self
        """
        if not isinstance(kmer, Kmer): raise TypeError("`kmer` must be a `Kmer` instance")
        if self.kmer_max_length != kmer.kmer_max_length: raise ValueError("`kmer` must have the same `kmer_max_length`")
        if self.kgap_max_length != kmer.kgap_max_length: raise ValueError("`kmer` must have the same `kgap_max_length`")
        if self.sequence_type != kmer.sequence_type: raise ValueError("`kmer` must have the same `sequence_type`")

    def __add__(self, kmer):
        """
        Add kmer counts and normalised kmer frequencies in self to kmer
        """
        self._type_match_check(kmer)
        for k in range(1, self.kmer_max_length + 1):
            self.kmer_arrays[k] = self.kmer_arrays[k] + kmer.kmer_arrays[k]
            self.kmer_normalised_arrays[k] = self.kmer_normalised_arrays[k] + kmer.kmer_normalised_arrays[k]
        for k in range(3, self.kgap_max_length + 1):
            self.kgap_arrays[k] = self.kgap_arrays[k] + kmer.kgap_arrays[k]
            self.kgap_normalised_arrays[k] = self.kgap_normalised_arrays[k] + kmer.kgap_normalised_arrays[k]
        return self
    
    def __sub__(self, kmer):
        """
        Subtracts kmer counts and normalised kmer frequencies in kmer from self
        """
        self._type_match_check(kmer)
        for k in range(1, self.kmer_max_length + 1):
            self.kmer_arrays[k] = self.kmer_arrays[k] - kmer.kmer_arrays[k]
            self.kmer_normalised_arrays[k] = kmer.kmer_normalised_arrays[k] - self.kmer_normalised_arrays[k]
        for k in range(3, self.kgap_max_length + 1):
            self.kgap_arrays[k] = self.kgap_arrays[k] - kmer.kgap_arrays[k]
            self.kgap_normalised_arrays[k] = self.kgap_normalised_arrays[k] - kmer.kgap_normalised_arrays[k]
        return self

    def __mul__(self, value):
        """
        Multiply kmer counts and normalised kmer frequencies by a number
        """
        if not isinstance(value, int) and not isinstance(value, float): raise TypeError("`value` must be an `int` or `float`")
        for k in range(1, self.kmer_max_length + 1):
            self.kmer_arrays[k] = self.kmer_arrays[k] * value
            self.kmer_normalised_arrays[k] = self.kmer_normalised_arrays[k] * value
        for k in range(3, self.kgap_max_length + 1):
            self.kgap_arrays[k] = self.kgap_arrays[k] + value
            self.kgap_normalised_arrays[k] = self.kgap_normalised_arrays[k] * value
        return self

    def __truediv__(self, value):
        """
        Divide kmer counts and normalised kmer frequencies by a number
        """
        if not isinstance(value, int) and not isinstance(value, float): raise TypeError("`value` must be an `int` or `float`")
        for k in range(1, self.kmer_max_length + 1):
            self.kmer_arrays[k] = self.kmer_arrays[k] / value
            self.kmer_normalised_arrays[k] = self.kmer_normalised_arrays[k] / value
        for k in range(3, self.kgap_max_length + 1):
            self.kgap_arrays[k] = self.kgap_arrays[k] / value
            self.kgap_normalised_arrays[k] = self.kgap_normalised_arrays[k] / value
        return self

    def abs(self):
        """
        Return the absolute value of kmer and kgap array entries
        """
        for k in range(1, self.kmer_max_length + 1):
            self.kmer_arrays[k] = numpy.abs(self.kmer_arrays[k])
            self.kmer_normalised_arrays[k] = numpy.abs(self.kmer_normalised_arrays[k])
        for k in range(3, self.kgap_max_length + 1):
            self.kgap_arrays[k] = numpy.abs(self.kgap_arrays[k])
            self.kgap_normalised_arrays[k] = numpy.abs(self.kgap_normalised_arrays[k])
        return self

    def _kmer_array_to_dict(self, kmer_arrays):
        """
        Dict formatted version of kmer_arrays or kmer_normalised_arrays
        """
        kmers = {}
        for k in range(1, self.kmer_max_length + 1):
            all_kmers = self.alphabet.all_kmers(k)
            kmers[k] = {}
            for i in range(len(all_kmers)):
                kmers[k][all_kmers[i]] = kmer_arrays[k][i]
        return kmers

    def kmers(self):
        """
        Dict formatted version of kmer_arrays
        """
        return self._kmer_array_to_dict(self.kmer_arrays)
    
    def kmers_normalised(self):
        """
        Dict formatted version of kmer_normalised_arrays
        """
        return self._kmer_array_to_dict(self.kmer_normalised_arrays)

    def _kgap_array_to_dict(self, kgap_arrays):
        """
        Dict formatted version of kgap_arrays or kgap_normalised_arrays
        """
        kgaps = {}
        all_kmers = self.alphabet.all_kmers(2)
        for k in range(3, self.kgap_max_length + 1):
            kgaps[k] = {}
            for i in range(len(all_kmers)):
                kgaps[k][all_kmers[i]] = kgap_arrays[k][i]
        return kgaps

    def kgaps(self):
        """
        Dict formatted version of kgap_arrays
        """
        return self._kgap_array_to_dict(self.kgap_arrays)
    
    def kgaps_normalised(self):
        """
        Dict formatted version of kgap_normalised_arrays
        """
        return self._kgap_array_to_dict(self.kgap_normalised_arrays)

    def kmer_difference(self, query):
        """
        Difference between self and query, using kmer frequency normalised to length (kmer_normalised).

        Positional arguments:
        query -- Kmer instance to which to compare

        Returns:
        Returns `Kmer` instance with change in kmer proportion (query-self) and change normalised frequence in self ((query-self)/self)
        """
        change = self - query
        for k in range(1, self.kmer_max_length + 1):
            for i in range(len(change.kmer_arrays)):
                # change_normalised: if self is non-zero, (query - self) / self
                if self.kmer_arrays[k][i] == 0:
                    change.kmer_normalised_arrays[k][i] = 0
                else:
                    change.kmer_normalised_arrays[k][i] = change.kmer_arrays[k][i] / self.kmer_arrays[k][i]
        for k in range(3, self.kgap_max_length + 1):
            for i in range(len(change.kgap_arrays)):
                # change_normalised: if self is non-zero, (query - self) / self
                if self.kgap_arrays[k][i] == 0:
                    change.kgap_normalised_arrays[k][i] = 0
                else:
                    change.kgap_normalised_arrays[k][i] = change.kgap_arrays[k][i] / self.kgap_arrays[k][i]
        return change

    def from_sequence(self, sequence: str):
        """
        Determines kmers from a sequence, populating self.kmer_arrays, self.kmer_normalised_arrays, self.kgap_arrays and self.kgap_normalised_arrays

        Required arguments:
        sequence -- sequence as a string
        """

        self.sequence = sequence
        # clear arrays
        for k in range(1, self.kmer_max_length + 1):
            kmers = self.alphabet.alphabet_length ** k
            self.kmer_arrays[k] = numpy.zeros(kmers, dtype=float)
            self.kmer_normalised_arrays[k] = numpy.zeros(kmers, dtype=float)
        # for each position in sequence
        for i in range(len(self.sequence) - 1):
            # for each kmer length
            # A, AA, AAA, etc.
            for k in range(1, self.kmer_max_length + 1):
                # if kmer is fully within the sequence length
                if i + k < len(self.sequence):
                    # get kmer
                    characters = self.sequence[i:i + k]
                    # determine array index
                    array_index = 0
                    for j in range(k):
                        array_index += self.alphabet.alphabet_lookup[characters[j]] * (self.alphabet.alphabet_length ** (k - 1 - j))
                    # record simple count
                    self.kmer_arrays[k][array_index] += 1
            # for each kgap length
            # AXA, AXXA, AXXXA, etc.
            for k in range(3, self.kgap_max_length + 1):
                # if kgap is fully within the sequence length
                if i + k < len(self.sequence):
                    # get kgap
                    characters = self.sequence[i] + self.sequence[i + k]
                    # determine array index
                    array_index = 0
                    for j in range(2):
                        array_index += self.alphabet.alphabet_lookup[characters[j]] * (self.alphabet.alphabet_length ** (2 - 1 - j))
                    # record simple count
                    self.kgap_arrays[k][array_index] += 1
        for k in range(1, self.kmer_max_length + 1):
            # record normalised count, kmers per number of kmers of that length within a sequence of that length
            self.kmer_normalised_arrays[k] = self.kmer_arrays[k] / (len(sequence) - k)
        for k in range(3, self.kgap_max_length + 1):
            # kgaps per number of kgaps of that length within a sequence of that length
            self.kgap_normalised_arrays[k] = self.kgap_arrays[k] / (len(sequence) - k)

class Sequence:
    def __init__(self, sequence: str, kmer_max_length: int = None, kgap_max_length: int = None, sequence_type: str = "protein", entropy_threshold: float = None, entropy_window: int = None):
        """
        Initialise class containing compositional information for kmers of different lengths
        
        Required arguments:
        sequence -- the sequence, can be protein, DNA or RNA

        Keyword arguments:
        sequence_type -- "protein", "dna" or "rna" sequence type, defining expected alphabet (default `"protein"`)
        max_kmer_length -- analyse composition in kmers up to and equal to this length (default `3`)
        """

        # record settings
        self.entropy_threshold = entropy_threshold
        self.entropy_window = entropy_window

        # default entropy window of 9
        if self.entropy_threshold is not None:
            if self.entropy_window is None:
                self.entropy_window = 9
            if self.entropy_window % 2 == 0:
                raise ValueError("Window for Shannon entropy calculation must be odd")
        
        # setup alphabet
        self.sequence_type = "".join(sequence_type.lower().split())
        self.alphabet = Alphabet(self.sequence_type)

        # setup kmer and kgap max lengths
        if kmer_max_length is None or kgap_max_length is None:
            default_kmer_max_length, default_kgap_max_length = self.alphabet.default_k()
            if kmer_max_length is None: kmer_max_length = default_kmer_max_length
            if kgap_max_length is None: kgap_max_length = default_kgap_max_length
        self.kmer_max_length = kmer_max_length
        self.kgap_max_length = kgap_max_length
        
        # check and tidy sequence input
        self.sequence = ""
        self.sequence_any_count = 0
        for char in sequence.upper():
            if char in self.alphabet.alphabet_replace:
                #print(char, self.alphabet.alphabet_replace)
                self.sequence = self.sequence + self.alphabet.alphabet_replace[char]
            elif char not in self.alphabet.alphabet:
                self.sequence = self.sequence + self.alphabet.alphabet_any
                self.sequence_any_count += 1
            else:
                self.sequence = self.sequence + char
        if len(self.sequence) == 0:
            warnings.warn("`sequence` length is zero, is the input formatted correctly or the entropy threshold too low?")
        elif self.sequence_any_count / len(self.sequence) > 0.1:
            warnings.warn("`sequence` contains more than 10% entries not in the " + self.sequence_type + " alphabet, is it the correct sequence type?")

    def __repr__(self):
        """
        Represent Sequence as the string serialisation of self.kmers
        """
        return str(self.kmers)

    def __str__(self):
        """
        Human-readable Sequence description including sequence, sequence type and common kmers
        """
        # simple description
        string = "KMERS within the " + self.sequence_type + " sequence:\n"
        string += self.sequence + "\n"
        # show entropy thresholding
        if self.entropy_threshold is not None:
            string += "Sequence filtered to Shannon entropy under " + str(self.entropy_threshold) + "\n"
            string += self.entropy_filtered_sequence + "\n"
        # append common kmers
        string += str(self.kmers_instance)
        return string[:-1]

    @cached_property
    def kmers_instance(self):
        """
        Instance of Kmer containing the actual kmer data
        """
        kmers_instance = Kmer(self.kmer_max_length, self.kgap_max_length, self.sequence_type)
        kmers_instance.from_sequence(self.entropy_filtered_sequence)
        return kmers_instance

    @cached_property
    def kmers(self):
        """
        Analysis of frequency of kmers of different lengths in the sequence.
        Only retains parts of the sequence below the entropy threshold, if set.

        Returns:
        object in the form {kmer_length: {kmer_sequence: kmer_count, ...}, ...}
        """
        return self.kmers_instance.kmers()
    
    @cached_property
    def kmers_normalised(self):
        """
        Analysis of kmers, but normalised to total number of kmers of that length in the sequence.

        Returns:
        object in the form {kmer_length: {kmer_sequence: kmer_proportion, ...}, ...}
        """
        return self.kmers_instance.kmers_normalised()

    @cached_property
    def kgaps(self):
        """
        Analysis of frequency of kgaps of different lengths in the sequence.
        Only retains parts of the sequence below the entropy threshold, if set.

        Returns:
        object in the form {kgap_length: {kgap_sequence: kgap_count, ...}, ...}
        """
        return self.kmers_instance.kgaps()
    
    @cached_property
    def kgaps_normalised(self):
        """
        Analysis of kgaps, but normalised to total number of kgaps of that length in the sequence.

        Returns:
        object in the form {kgap_length: {kgap_sequence: kgap_count, ...}, ...}
        """
        return self.kmers_instance.kgaps_normalised()
    
    @cached_property
    def entropy_filtered_sequence(self):
        """
        Sequence filtered to only include positions under the threshold Shannon entropy

        Returns:
        Concatenated sequence fragments under the threshold Shannon entropy
        """
        if self.entropy_threshold is not None:
            filtered_sequence = ""
            for i in range(len(self.sequence)):
                if self.entropy[i] < self.entropy_threshold: filtered_sequence += self.sequence[i]
            return filtered_sequence
        else:
            return self.sequence

    @cached_property
    def entropy(self):
        """
        Rolling window Shannon entropy of the protein/DNA/RNA sequence

        Returns:
        List of entropy per position, based on window centre, and extrapolated to the full sequence length
        """
        entropy = []
        for i in range(len(self.sequence) - self.entropy_window + 1):
            # calculate entropy per window position
            # https://saturncloud.io/blog/calculating-shannon-entropy-of-an-array-using-pythons-numpy/
            array = numpy.array([x for x in self.sequence[i:i + self.entropy_window]])
            unique, counts = numpy.unique(array, return_counts = True)
            probabilities = counts / self.entropy_window
            entropy.append(-numpy.sum(probabilities * numpy.log2(probabilities)))
        # extrapolate start and end to full sequence length
        output = [entropy[0]] * int(((self.entropy_window - 1) / 2)) + entropy + [entropy[-1]] * int(((self.entropy_window - 1) / 2))
        return output

    def _prepare_matrix_chisq(self, matrix = None):
        """
        Prepare a matrix for chi squared tests.
        """
        # TODO: Handle zeroes in the matrix - any zeroes in the matrix break Chi squared
        # if not set, default to a unity composition matrix (ie. comparison to uniform kmer distribution)
        if matrix is None:
            matrix_instance = KmerMatrix(sequence_type = self.sequence_type, kmer_max_length = self.kmer_max_length)
            matrix = matrix_instance.composition_unity
        # check if the kmers match type
        self.kmers_instance._type_match_check(matrix)
        # ensure all kmers frequencies are non-zero, by changing zero to 1/100 minimum
        for k in range(1, self.kmer_max_length + 1):
            min_frequency = numpy.min(matrix.kmer_normalised_arrays[k])
            matrix.kmer_normalised_arrays[k][matrix.kmer_normalised_arrays == 0] = min_frequency / 100
        for k in range(3, self.kgap_max_length + 1):
            min_frequency = numpy.min(matrix.kgap_normalised_arrays[k])
            matrix.kgap_normalised_arrays[k][matrix.kgap_normalised_arrays == 0] = min_frequency / 100
        return matrix

    def _kmer_outlier_stats(self, matrix = None):
        """
        Statistical comparison of each individual kmer with a reference matrix.

        Named arguments:
        matrix -- Kmer instance encoding composition matrix (defauly, unity composition matrix)

        Returns:
        Per kmer, [fold change, p value] where:
        fold change is observed kmer count / expected kmer
        p value from chi squared test of count of that kmer vs. all others of the same length, in comparison to `matrix`.
        """
        matrix = self._prepare_matrix_chisq()
        # do the stats
        enrichment = {"kmer": {}, "kgap": {}}
        for k in range(1, self.kmer_max_length + 1):
            enrichment["kmer"][k] = {}
            # for every kmer length
            for i in range(self.alphabet.num_kmers(k)):
                # for each individual kmer
                observed = [self.kmers_instance.kmer_arrays[k][i]]
                observed.append(sum(self.kmers_instance.kmer_arrays[k] - observed[0]))
                observed = numpy.array(observed)
                expected = numpy.array([matrix.kmer_normalised_arrays[k][i], 1 - matrix.kmer_normalised_arrays[k][i]]) * sum(observed)
                statistic, pvalue = stats.chisquare(observed, expected)
                enrichment["kmer"][k][i] = [observed[0] / expected[0], pvalue]
        for k in range(3, self.kgap_max_length + 1):
            enrichment["kgap"][k] = {}
            # for every kmer length
            for i in range(self.alphabet.num_kmers(2)):
                # for each individual kmer
                observed = [self.kmers_instance.kgap_arrays[k][i]]
                observed.append(sum(self.kmers_instance.kgap_arrays[k] - observed[0]))
                observed = numpy.array(observed)
                expected = numpy.array([matrix.kgap_normalised_arrays[k][i], 1 - matrix.kgap_normalised_arrays[k][i]]) * sum(observed)
                statistic, pvalue = stats.chisquare(observed, expected)
                enrichment["kgap"][k][i] = [observed[0] / expected[0], pvalue]
        return enrichment

    def kmer_outlier_stats(self, matrix = None, threshold = 0.05):
        """
        Statistical comparison of each individual kmer with a reference matrix.

        Named arguments:
        matrix -- Kmer instance encoding composition matrix (defauly, unity composition matrix)
        threshold -- minimum p value for reporting

        Returns:
        Human-readable summary of kmer enrichment
        """
        enrichment = self._kmer_outlier_stats(matrix)
        ranges = {"kmer": {"min": 1, "max": self.kmer_max_length + 1}, "kgap": {"min": 3, "max": self.kgap_max_length + 1}}
        for type in ["kmer", "kgap"]:
            for k in range(ranges[type]["min"], ranges[type]["max"]):
                print(type + " length " + str(k))
                print("    [fold change]\t[p value]")
                if type == "kmer":
                    all_kmers = self.alphabet.all_kmers(k)
                elif type == "kgap":
                    all_kmers = self.alphabet.all_kmers(2)
                count = 0
                for i in range(len(all_kmers)):
                    if enrichment[type][k][i][1] <= 0.05:
                        print("    " + all_kmers[i] + "\t" +str(enrichment[type][k][i][0]) + "\t" + str(enrichment[type][k][i][1]))
                        count += 1
                if count == 0:
                    print("    [none]")

    def composition_dissimilarity_stats(self, matrix = None):
        """
        Statistical comparison of sequence kmers with a reference matrix.

        Named arguments:
        matrix -- Kmer instance encoding composition matrix (default, unity composition matrix)
        
        Returns:
        Chi squared test results for each kmer length `k`, in the form `{k: p_value, ...}`
        """
        matrix = self._prepare_matrix_chisq()
        # do the stats
        enrichment = {}
        for k in range(1, self.kmer_max_length + 1):
            statistic, enrichment[k] = stats.chisquare(self.kmers_instance.kmer_arrays[k], matrix.kmer_normalised_arrays[k] * sum(self.kmers_instance.kmer_arrays[k]))
        return enrichment

    def kmer_similarity(self, kmers, matrix):
        """
        Similarity score of self.kmer with kmers, using a composition change matrix
        """
        a=1
    
    def kmer_difference(self, query):
        """
        Difference between self.kmers_instance and query, using kmer frequency normalised to length (kmer_normalised).

        Positional arguments:
        query -- Kmer instance to which to compare

        Returns:
        Returns `Kmer` instance with change in kmer proportion (query-self) and change normalised frequence in self ((query-self)/self)
        """
        return self.kmers_instance.kmer_difference(query)

    def sequence_difference(self, query):
        """
        Difference between self.kmers_instance and query.kmers_instance.
        
        Positional arguments:
        query -- Sequence instance to which to compare

        Returns:
        Returns `Kmer` instance with change in kmer proportion (query-self) and change normalised frequence in self ((query-self)/self)
        """
        return self.kmers_instance.kmer_difference(query.kmers_instance)

    def plot(self, save_path: str = None, show_plot: bool = False):
        """
        Plot frequency of kmers and kgaps of different lengths in the sequence

        Named arguments:
        save_path -- save plots to this base file name (default is `None`, do not save)
        show_plot -- show plot, when using a GUI backend (default `False`)
        """
        plot = _Plot(self)
        plot.plot(save_path = save_path, show_plot = show_plot)

class KmerMatrix:
    def __init__(self, sequence_type: str = "protein", kmer_max_length: int = 3, kgap_max_length: int = 5):
        """Initialise object containing kmer composition and kmer substitution rate information.

        Keyword arguments:
        sequence_type -- "protein", "dna" or "rna" sequence type, defining expected alphabet (default `"protein"`)
        max_kmer_length -- analyse composition in kmers up to and equal to this length (default `3`)

        All substitution matrices are in the form:
        {kmer_length: {original_kmer: {substituted_kmer: substitution_fraction, ...}, ...}, ...}
        """

        # record settings
        self.kmer_max_length = kmer_max_length
        self.kgap_max_length = kgap_max_length
        
        # setup alphabet
        self.sequence_type = "".join(sequence_type.lower().split())
        self.alphabet = Alphabet(self.sequence_type)

    def _unity(self):
        """
        Generates a normalised unity matrix.
        Dummy data (all ones) in Kmer.kmer_arrays and Kmer.kmer_normalised_arrays
        """
        matrix = Kmer(kmer_max_length = self.kmer_max_length, kgap_max_length = self.kgap_max_length, sequence_type = self.sequence_type)
        matrix.kmer_arrays = {}
        matrix.kmer_normalised_arrays = {}
        for k in range(1, self.kmer_max_length + 1):
            all_kmers = self.alphabet.all_kmers(k)
            matrix.kmer_arrays[k] = numpy.ones(len(all_kmers), dtype=float) * len(self.alphabet.alphabet) ** k
            matrix.kmer_normalised_arrays[k] = numpy.ones(len(all_kmers), dtype=float) / len(all_kmers)
        for k in range(3, self.kgap_max_length + 1):
            all_kmers = self.alphabet.all_kmers(2)
            matrix.kgap_arrays[k] = numpy.ones(len(all_kmers), dtype=float) * len(self.alphabet.alphabet) ** 2
            matrix.kgap_normalised_arrays[k] = numpy.ones(len(all_kmers), dtype=float) / len(all_kmers)
        return matrix

    @cached_property
    def composition_unity(self):
        """
        Unity composition matrix, all kmers are equally prevalent.
        Normalised, such that prevalence of all kmers of a particular length sum to 1.
        """
        return self._unity()
    
    @cached_property
    def change_unity(self):
        """
        Unity change matrix, all changes to kmers are equally prevalent.
        Normalised, such that changes of all kmers of a particular length sum to 1.
        """
        return self._unity()
    
    def _composition_from_fasta(self, fasta):
        """
        Average composition from all sequences in a fasta instance.
        """
        average = None
        for name in fasta.sequences:
            kmer = Sequence(fasta.sequences[name], kmer_max_length=self.kmer_max_length, sequence_type=self.sequence_type)
            if average is None:
                average = kmer.kmers_instance
            else:
                average = average + kmer.kmers_instance
        return average / fasta.count

    def composition_from_fasta(self, path: str):
        """
        Average composition matrix from all sequences in a fasta file.
        """
        fasta = Fasta(path)
        return self._composition_from_fasta(fasta)

    def composition_from_fastas(self, path_list: list):
        """
        Average composition from all sequences in a list of fasta files.
        """
        average = None
        for path in path_list:
            fasta = Fasta(path)
            current_average = self._composition_from_fasta(fasta)
            if average is None:
                average = current_average / len(path_list)
            else:
                average = average + current_average / len(path_list)
        return average
    
    def _change_from_fasta(self, fasta):
        """
        Average composition change from pairwise comparison of all sequences in a `Fasta` instance.
        """
        change = None
        count = 0
        for name1 in fasta.sequences:
            kmer1 = Sequence(fasta.sequence[name1], kmer_max_length=self.kmer_max_length, sequence_type=self.sequence_type)
            for name2 in fasta.sequences:
                kmer2 = Sequence(fasta.sequence[name2], kmer_max_length=self.kmer_max_length, sequence_type=self.sequence_type)
                if name1 != name2:
                    kmer_difference = kmer1 - kmer2
                    change = change + kmer_difference.abs()
                    count += 1
        return change / count
    
    # Additional matrices... standard pre-calculated matrices, tools for calculating matrices

class _Plot:
    def __init__(self, kmers: Kmer, plot_font_size = 8, plot_dpi = 100):
        """Initialise object for plotting of kmer-like data

        Required arguments:
        kmer -- Sequence instance to plot
        """
        self.kmers = kmers
        self.plot_font_size = plot_font_size
        self.plot_dpi = plot_dpi

    def _plot_1d(self, list, aas, labels, xlabel: str, ylabel: str, save_name: str = None, show_plot: bool = False, title: str = None):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(3 * len(labels) / len(aas), 2), dpi=self.plot_dpi)
        ax=plt.subplot(1, 1, 1)
        plt.bar(aas, list)
        plt.xticks(rotation=90, fontsize=self.plot_font_size)
        plt.yticks(fontsize=self.plot_font_size)
        if title != None:
            plt.title = title
        ax.set_xlabel(xlabel, fontsize=self.plot_font_size)
        ax.set_ylabel(ylabel, fontsize=self.plot_font_size)
        if show_plot: plt.show()
        if save_name != None: plt.savefig("%s.png" % save_name)
        plt.close()

    def _plot_2d(self, matrix, aas, labelsx, labelsy, xlabel: str, ylabel: str, zlabel, save_name: str = None, show_plot: bool = False, title: str = None):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4 * len(labelsx) / len(aas), 3.5), dpi=self.plot_dpi)
        ax=plt.subplot(1, 1, 1)
        plt.imshow(matrix, cmap="plasma")
        plt.xticks(range(len(labelsx)), labelsx, rotation=90, fontsize=self.plot_font_size)
        plt.yticks(range(len(labelsy)), labelsy, fontsize=self.plot_font_size)
        ax.set_xlabel(xlabel, fontsize=self.plot_font_size)
        ax.set_ylabel(ylabel, fontsize=self.plot_font_size)
        if title != None:
            plt.title = title
        clb=plt.colorbar()
        clb.ax.tick_params(labelsize=self.plot_font_size)
        clb.ax.set_title(zlabel, fontsize=self.plot_font_size)
        if show_plot: plt.show()
        if save_name != None: plt.savefig("%s.png" % save_name)
        plt.close()

    def plot(self, save_path: str = None, show_plot: bool = False):
        """
        Plot frequency of kmers of different lengths in the sequence

        Named arguments:
        save_path -- save plots to this base file name (default is `None`, do not save)
        show_plot -- show plot, when using a GUI backend (default `False`)
        """

        aa = self.kmers.alphabet.all_kmers(1)
        for k in range(1, self.kmers.kmer_max_length + 1):
            if k == 1:
                # length 1, ie. single amino acids
                # bar plot
                matrix = numpy.zeros(len(aa))
                index = 0
                for entry in self.kmers.kmers[k]:
                    matrix[index] = self.kmers.kmers[k][entry]
                    index += 1
                current_path = None
                if save_path != None: current_path = save_path + ".k1"
                self._plot_1d(matrix, self.kmers.alphabet.all_kmers(1), self.kmers.alphabet.all_kmers(1), "amino acid", "frequency", save_name = current_path, show_plot = show_plot, title = "Query kmers, k = 1")
            elif k == 2:
                # length 2
                # matrix heat map
                matrix = numpy.zeros((len(aa), len(aa)))
                index = 0
                for entry in self.kmers.kmers[k]:
                    matrix[index % len(aa)][math.floor(index/len(aa))] = self.kmers.kmers[k][entry]
                    index += 1
                current_path = None
                if save_path != None: current_path = save_path + ".k2"
                self._plot_2d(matrix, self.kmers.alphabet.all_kmers(1), self.kmers.alphabet.all_kmers(1), self.kmers.alphabet.all_kmers(1), "amino acid 1", "amino acid 2", "frequency", save_name = current_path, show_plot = show_plot, title = "Query kmers k = 2")
        for k in range(3, self.kmers.kgap_max_length + 1):
            # always length 2
            # matrix heat map
            matrix = numpy.zeros((len(aa), len(aa)))
            index = 0
            for entry in self.kmers.kgaps[k]:
                matrix[index % len(aa)][math.floor(index/len(aa))] = self.kmers.kgaps[k][entry]
                index += 1
            current_path = None
            if save_path != None: current_path = save_path + ".g" + str(k)
            self._plot_2d(matrix, self.kmers.alphabet.all_kmers(1), self.kmers.alphabet.all_kmers(1), self.kmers.alphabet.all_kmers(1), "amino acid 1", "amino acid 2", "frequency", save_name = current_path, show_plot = show_plot, title = "Query kmers k = 2")

class Search:
    def __init__(self, sequence_type: str = "protein", kmer_max_length: int = 3, kgap_max_length: int = 5, entropy_threshold: float = None, entropy_window: int = 9):
        """Initialise object for carrying out composition-based searches

        Keyword arguments:
        sequence_type -- "protein", "dna" or "rna" sequence type, defining expected alphabet (default `"protein"`)
        kmer_max_length -- analyse composition in kmers from 1 up to and equal to this length (default `3`)
        kgap_max_length -- analyse composition in kgaps from 3 up to and equal to this length (default `5`)
        """

        # record settings
        self.kmer_max_length = kmer_max_length
        self.kgap_max_length = kgap_max_length
        self.entropy_threshold = entropy_threshold
        self.entropy_window = entropy_window
        
        # setup alphabet
        self.sequence_type = "".join(sequence_type.lower().split())
        self.alphabet = Alphabet(self.sequence_type)
    
        self.query = None
        self.database = None

    def new_database(self, sequences):
        """
        Creates the database of composition from every provided sequence

        Required arguments:
        sequences -- dict of name/sequences, in the form {"sequence_name": "sequence", ...}
        """
        self.database = []
        for name in sequences:
            self.database.append({"name": name, "kmers": Sequence(sequences[name], kmer_max_length = self.kmer_max_length, sequence_type = self.sequence_type).kmers})

    def new_database_from_fasta(self, path: str):
        """
        Creates the database of composition from a fasta file

        Required arguments:
        path -- path of the fasta file to load
        """
        fasta = Fasta(path)
        self.new_database(fasta.sequences)

    _string_columnindividual_kmer_cutoff = 3

    def _format_kmers_string(self, kmers):
        """
        Format a kmers object as a machine-readable string

        Required arguments:
        kmers -- kmers instance

        Returns a string representation. Columns of counts for each kgap then each kmer, in k then alphabet order for k < 4, kmer:count pairs for k >= 4
        """
        #TODO: refactor to write from a Kmer instance
        string = ""
        for k in range(3, self.kgap_max_length + 1):
            all_kmers = self.alphabet.all_kmers(2)
            for i in range(self.alphabet.num_kmers(2)):
                string += str(self.kgaps[k][all_kmers[i]]) + "\t"
        for k in range(1, self.kmer_max_length + 1):
            all_kmers = self.alphabet.all_kmers(k)
            for i in range(self.alphabet.num_kmers(k)):
                if k < self._string_columnindividual_kmer_cutoff:
                    # for short kmers, record all counts as columns
                    string += str(kmers[k][all_kmers[i]]) + "\t"
                elif kmers[k][all_kmers[i]] != 0:
                    # for long kmers, record kmers with non-zero counts as kmer:count pairs
                    string += all_kmers[i] + ":" + str(kmers[k][all_kmers[i]]) + "\t"
        return string

    def _parse_kmers_string(self, string: str):
        """
        Load a kmers object from a machine-readable string

        Required arguments:
        string -- string encoding kmers object

        Returns a kmers object
        """
        #TODO: Load from string into a `Kmer` instance

    def save_database(self, path: str, md5: str = None):
        """
        Saves the database of composition in a text format to a file

        Required arguments:
        path -- path to save the database

        Named arguments:
        md5 -- md5 of the input file
        """
        with open(path, "w") as output_file:
            # header
            output_file.write("cki_database_version\t1\n")
            output_file.write("sequence_type\t" + self.sequence_type + "\n")
            output_file.write("kmer_max_length\t" + str(self.kmer_max_length) + "\n")
            output_file.write("kgap_max_length\t" + str(self.kmer_max_length) + "\n")
            if self.entropy_threshold is not None:
                output_file.write("entropy_threshold\t" + str(self.entropy_threshold) + "\n")
                output_file.write("entropy_window\t" + str(self.entropy_window) + "\n")
            if md5 is not None:
                output_file.write("source_md5\t" + md5 + "\n")
            # end of header indicator
            output_file.write("###\n")
            # data
            for entry in self.database:
                output_file.write(entry["name"] + "\t" + self._format_kmers_string(entry["kmers"]) + "\n")
            # end of file indicator
            output_file.write("###\n")

    def _read_index_header(self, path: str):
        """
        Reads index properties from the index header

        Required arguments:
        path -- path of the .cki index file to read

        Returns
        cki_database_version, sequence_type, kmer_max_length and source_md5
        """
        with open(path, "r") as index_file:
            cki_database_version, sequence_type, kmer_max_length, source_md5, entropy_window, entropy_threshold = None, None, None, None, None, None
            while True:
                line = index_file.readline()
                if not line or line == "###":
                    break
                line = line.split("\t")
                if line[0] == "cki_database_version": cki_database_version = int(line[1])
                if line[0] == "sequence_type": sequence_type = line[1]
                if line[0] == "kmer_max_length": kmer_max_length = int(line[1])
                if line[0] == "kgap_max_length": kgap_max_length = int(line[1])
                if line[0] == "entropy_window": entropy_window = int(line[1])
                if line[0] == "entropy_threshold": entropy_threshold = int(line[1])
                if line[0] == "source_md5": source_md5 = line[1]
        return cki_database_version, sequence_type, kmer_max_length, kgap_max_length, source_md5, entropy_window, entropy_threshold

    def file_md5(self, path: str):
        def file_as_bytes(file):
            with file:
                return file.read()
        return hashlib.md5(file_as_bytes(open(path, 'rb'))).hexdigest()

    def index_fasta_file(self, path: str):
        """
        Creates a database from a fasta file and saves the result.
        Checks the md5 of the input fasta file and does not update the output if the md5 matches an existing saved database (index).

        Required arguments:
        path -- path of the fasta file to load
        """
        do_indexing = True
        md5 = None
        # check the header of an existing .cki file and overwrite if necessary
        if os.path.isfile(path + ".cki"):
            do_indexing = False
            cki_database_version, sequence_type, kmer_max_length, kgap_max_length, source_md5, entropy_window, entropy_threshold = self._read_index_header(path)
            # check for matching versions and types
            if cki_database_version != 1: do_indexing = True
            if sequence_type != self.sequence_type: do_indexing = True
            if kmer_max_length != self.kmer_max_length: do_indexing = True
            if kgap_max_length != self.kgap_max_length: do_indexing = True
            if entropy_window != self.entropy_window: do_indexing = True
            if entropy_threshold != self.entropy_threshold: do_indexing = True
            # if necessary, check md5
            if not do_indexing:
                md5 = self.file_md5(path)
                if source_md5 != md5: do_indexing = True
        # do indexing
        if do_indexing:
            if md5 is None: md5 = self.file_md5(path)
            self.new_database_from_fasta(path)
            self.save_database(path + ".cki", md5 = md5) # composition kmer index

    def load_database(self, path: str):
        """
        Load a database from an index file.

        Required arguments:
        path -- path of the .cki index file to load
        """
        cki_database_version, sequence_type, kmer_max_length, source_md5, entropy_window, entropy_threshold = self._read_index_header(path)
        if cki_database_version != 1: raise ValueError("Unknown database version " + str(cki_database_version) + " encountered loading .cki file")
        if kmer_max_length != self.kmer_max_length: raise ValueError("Database constructed for kmers of length " + str(kmer_max_length) + " and we want " + str(self.kmer_max_length))
        if kmer_max_length != self.kmer_max_length: raise ValueError("Database constructed for sequences of type " + sequence_type + " and we want " + self.sequence_type)
        if entropy_threshold is not None:
            if entropy_threshold != self.entropy_threshold: raise ValueError("Database constructed with entropy filtering threshold of " + str(entropy_threshold) + " and we want " + str(self.entropy_threshold))
            if entropy_window != self.entropy_window: raise ValueError("Database constructed with entropy filtering window of " + str(entropy_window) + " and we want " + str(self.entropy_window))
        database = []
        with open(path, "r") as index_file:
            lines = index_file.readlines()
            # skip header, then pass valid lines
            in_header = True
            for line in lines:
                if line == "###": in_header = False
                if not in_header and line != "###":
                    line = line.split("\t")
                    database.append({"name": line[0], "kmers": self._parse_kmers_string("\t".join(line[1:]))})
            