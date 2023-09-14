import subprocess
import os

class Blast:
    def __init__(self, search_tool: str = "diamond", sequence_type: str = "protein"):
        """
        Initialise a sequence search instance

        Named arguments:
        search_tool -- the search tool to use, either `"blast"` or `"diamond"`
        sequence_type -- the sequence type, either `"nucleotide"` or `"protein"`
        """
        # check search_tool value is valid
        search_tools = ["diamond", "blast"]
        self.search_tool = search_tool
        if self.search_tool not in search_tools:
            raise ValueError("`search_tool` must be " + " or ".join(search_tools))
        # check sequence_type value is valid
        sequence_types = ["protein", "nucleotide"]
        self.sequence_type = sequence_type
        if self.sequence_type not in sequence_types:
            raise ValueError("`search_tool` must be " + " or ".join(sequence_types))
    
    def _check_install(self, verbose: bool = False):
        """
        Check if the necessary programs are installed, and thow an error if not.
        ncbi-blast+ is required for `"blast"` search_type
        diamond is required for `"diamond"` search_type

        Named arguments:
        verbose -- print a verbose output, boolean, default `False`
        """
        command = {
            "diamond": ["diamond", "--help"], 
            "blast": ["makeblastdb", "-help"]
        }
        proc = subprocess.Popen(command[self.search_tool], stderr = subprocess.PIPE, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
        result = proc.communicate()[0].decode("utf-8").splitlines()
        if proc.returncode != 0:
            print(result)
            assert AssertionError(command[self.search_tool] + " did not run correctly, is it installed and available in the system path?")
        if verbose:
            print("'" + " ".join(command[self.search_tool]) + "' ran successfully, " + self.search_tool + " appears to be installed correctly")

    def search(self, query_sequence, fasta_path, query_name: str = None):
        self._check_install()
        if query_name is None:
            query_name = "none"
        search_program = {
            "protein": "blastp",
            "nucleotide": "blastn"
        }
        if self.search_tool == "diamond":
            # make the database
            if not os.path.isfile(fasta_path + ".dmnd"):
                command = ["diamond", "makedb", fasta_path]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                result = proc.communicate()[0].decode("utf-8")
                #if proc.returncode != 0:
                #    print(result)
                #    assert AssertionError("diamond database generation failed")
            # do the search
            search_command = ["-d", fasta_path, "-e", "0.00001", "-q", "-", "--outfmt", "6", "sseqid", "evalue", "sseq"]
            command = ["diamond", search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            proc.stdin.write((">%s\n%s" % (query_name, query_sequence)).encode('utf-8'))
            result = proc.communicate()[0].decode("utf-8")
            #if proc.returncode != 0:
            #    print(result)
            #    assert AssertionError("diamond search failed")
        if self.search_tool == "blast":
            # make the database
            if not os.path.isfile(fasta_path + ".pdb"):
                dbtype = {
                    "protein": "prot",
                    "nucleotide": "nucl"
                }
                command = ["makeblastdb", "-dbtype", dbtype[self.sequence_type], "-in", fasta_path]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                result = proc.communicate()[0].decode("utf-8")
                #if proc.returncode != 0:
                #    print(result)
                #    assert AssertionError("blast database generation failed")
            # do the search
            search_command = ["-db", fasta_path, "-evalue", "0.00001", "-query", "-", "-outfmt", "6 sseqid evalue sseq"]
            command = [search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            proc.stdin.write((">%s\n%s" % (query_name, query_sequence)).encode('utf-8'))
            result = proc.communicate()[0].decode("utf-8")
            #if proc.returncode != 0:
            #    print(result)
            #    assert AssertionError("diamond search failed")
            print(result)
        result = result.splitlines()
        result = [x.split("\t") for x in result]
        return result