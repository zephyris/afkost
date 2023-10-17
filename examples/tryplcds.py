from tryptag import TrypTag
from afkost import Sequence
from afkost.databases import TriTrypDB
import numpy

tryptag = TrypTag(verbose=False)
tritrypdb = TriTrypDB()
fasta = tritrypdb.sequences("TbruceiTREU927")

localisations = ["glycosome", "cortical cytoskeleton", "nucleoplasm", "nucleolus", "axoneme", "paraflagellar rod", "mitochondrion", "endoplasmic reticulum"]
queries = {}
for localisation in localisations:
    queries[localisation] = tryptag.localisation_search(localisation, match_subterms=False)

queries["cytoplasm points"] = tryptag.localisation_search("cytoplasm", required_modifiers=["points", "strong"], match_subterms=False)
queries["cytoplasm not points"] = tryptag.localisation_search("cytoplasm", exclude_modifiers=["<10%", "weak", "points"], match_subterms=False)


for query in queries:
    result = []
    for hit in queries[query]:
        name = hit.gene_id + ":mRNA-p1"
        if name in fasta.sequences:
            #sequence = Sequence(fasta.sequences[name], kmer_max_length=1, entropy_threshold=2.2, entropy_window=9)
            sequence = Sequence(fasta.sequences[name], kmer_max_length=1, entropy_threshold=2.0, entropy_window=7)
            result.append(len(sequence.entropy_filtered_sequence) / len(sequence.sequence))
    result = numpy.array(result)
    print(query, len(result), numpy.mean(result), numpy.percentile(result, 50), numpy.percentile(result, 25), numpy.percentile(result, 75), numpy.percentile(result, 5), numpy.percentile(result, 95))
