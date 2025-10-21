"""
UE - Algorithm in Bioinformatics
CCA3 
Name: Harshvardhan Salunke
TY B.TECH CSE Ai & DS
PRN - 1032230341

"""

import string
import random
from collections import Counter

# --- Part A ---

class DNA:
    """
    Represents a DNA sequence, providing validation and basic stats.
    
    Uses a string to store the sequence and a dictionary for statistics.
    """
    
    VALID_NUCLEOTIDES = 'ATGC'
    
    def __init__(self, sequence: str):
        """
        Initializes and validates the DNA sequence.
        
        Args:
            sequence (str): The DNA sequence string.
        
        Raises:
            ValueError: If the sequence contains invalid nucleotides.
        """
        self.sequence = sequence.upper()
        self.validate()
        
    def validate(self):
        """Internal method to validate the sequence."""
        for nucleotide in self.sequence:
            if nucleotide not in self.VALID_NUCLEOTIDES:
                raise ValueError(f"Invalid nucleotide detected: '{nucleotide}'")
                
    def __str__(self) -> str:
        """Returns the sequence as a string."""
        return self.sequence
        
    def __len__(self) -> int:
        """Returns the length of the sequence."""
        return len(self.sequence)
        
    def nucleotide_count(self) -> dict:
        """
        Counts the occurrences of each valid nucleotide (A, T, G, C).
        
        Returns:
            dict: A dictionary with counts for each nucleotide.
        """
        # Using collections.Counter is efficient and clean
        counts = Counter(self.sequence)
        
        # Ensure all four bases are in the dict, even if count is 0
        for nuc in self.VALID_NUCLEOTIDES:
            if nuc not in counts:
                counts[nuc] = 0
                
        return counts

    def basic_statistics(self) -> dict:
        """
        Provides a basic statistical report of nucleotide frequencies.
        
        Returns:
            dict: A dictionary with nucleotide frequencies (percentages).
        """
        return calculate_frequencies(self.sequence)


# Question 2

def count_nucleotides(seq: str) -> dict:
    """
    Counts individual nucleotides in a DNA sequence.
    
    Args:
        seq (str): The DNA sequence.
        
    Returns:
        dict: Counts for 'A', 'T', 'G', 'C'.
    """
    # Simple dictionary-based counting
    seq = seq.upper()
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for nucleotide in seq:
        if nucleotide in counts:
            counts[nucleotide] += 1
    return counts

def calculate_frequencies(seq: str) -> dict:
    """
    Calculates nucleotide frequencies as percentages.
    
    Args:
        seq (str): The DNA sequence.
        
    Returns:
        dict: Frequencies (percentages) for each nucleotide.
    """
    seq = seq.upper()
    counts = count_nucleotides(seq)
    total_len = len(seq)
    
    if total_len == 0:
        return {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}
        
    frequencies = {}
    for nucleotide, count in counts.items():
        frequencies[nucleotide] = (count / total_len) * 100
        
    return frequencies

def generate_nucleotide_report(seq: str) -> str:
    """
    Generates a comprehensive nucleotide analysis report.
    
    Args:
        seq (str): The DNA sequence.
        
    Returns:
        str: A formatted report string.
    """
    seq = seq.upper()
    counts = count_nucleotides(seq)
    frequencies = calculate_frequencies(seq)
    total_len = len(seq)
    
    report = f"--- Nucleotide Analysis Report ---\n"
    report += f"Sequence Length: {total_len}\n"
    report += "\nCounts:\n"
    report += f"  A: {counts['A']}\n"
    report += f"  T: {counts['T']}\n"
    report += f"  G: {counts['G']}\n"
    report += f"  C: {counts['C']}\n"
    report += "\nFrequencies:\n"
    report += f"  A: {frequencies['A']:.2f}%\n"
    report += f"  T: {frequencies['T']:.2f}%\n"
    report += f"  G: {frequencies['G']:.2f}%\n"
    report += f"  C: {frequencies['C']:.2f}%\n"
    report += "----------------------------------\n"
    return report

def compare_composition(seq1: str, seq2: str):
    """
    Prints a comparison of the nucleotide composition of two sequences.
    
    Args:
        seq1 (str): First DNA sequence.
        seq2 (str): Second DNA sequence.
    """
    print("--- Composition Comparison ---")
    print("\nSequence 1 Report:")
    print(generate_nucleotide_report(seq1))
    print("\nSequence 2 Report:")
    print(generate_nucleotide_report(seq2))
    print("------------------------------")


# Question 3

def to_uppercase(seq: str) -> str:
    """Converts a DNA sequence to uppercase."""
    return seq.upper()

def to_lowercase(seq: str) -> str:
    """Converts a DNA sequence to lowercase."""
    return seq.lower()

def remove_non_nucleotides(seq: str, valid_chars: str = 'ATGC') -> str:
    """
    Removes non-nucleotide characters from a sequence.
    
    Args:
        seq (str): The input sequence.
        valid_chars (str, optional): A string of allowed characters. 
                                     Defaults to 'ATGC'.
    
    Returns:
        str: The cleaned sequence.
    """
    seq = seq.upper()
    # Use a list comprehension for a clean, readable filter
    cleaned = [nuc for nuc in seq if nuc in valid_chars]
    return "".join(cleaned)

def split_into_codons(seq: str) -> list:
    """
    Splits a long sequence into codons (groups of 3).
    
    Handles sequences that are not a multiple of 3.
    
    Args:
        seq (str): The DNA sequence.
        
    Returns:
        list: A list of codon strings.
    """
    codons = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        codons.append(codon)
        
    # Optional: check for incomplete final codon
    if len(codons) > 0 and len(codons[-1]) < 3:
        print(f"Warning: Final codon '{codons[-1]}' is incomplete.")
        
    return codons

def merge_fragments(fragments: list) -> str:
    """
    Merges multiple DNA fragments (strings) into a single sequence.
    
    Args:
        fragments (list): A list of DNA sequence strings.
        
    Returns:
        str: The combined, single sequence.
    """
    return "".join(fragments)

# --- Part B: Essential DNA Algorithms ---

# Question 4: DNA Transcription

def transcribe(seq: str, strand_type: str = 'coding') -> str:
    """
    Converts a DNA sequence to its corresponding RNA sequence.
    
    Args:
        seq (str): The DNA sequence.
        strand_type (str, optional): The strand to transcribe. 
                                     'coding' (T->U) or 'template' (A->U, T->A, C->G, G->C).
                                     Defaults to 'coding'.
                                     
    Returns:
        str: The RNA sequence.
        
    Raises:
        ValueError: If strand_type is not 'coding' or 'template'.
    """
    seq = seq.upper()
    
    if strand_type == 'coding':
        # Easiest way: just replace T with U
        return seq.replace('T', 'U')
        
    elif strand_type == 'template':
        # Map template DNA base to RNA base
        template_map = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
        rna_list = []
        for nuc in seq:
            if nuc not in template_map:
                raise ValueError(f"Invalid nucleotide for template: {nuc}")
            rna_list.append(template_map[nuc])
        return "".join(rna_list)
        
    else:
        raise ValueError("strand_type must be 'coding' or 'template'")

def batch_transcribe(sequences: list, strand_type: str = 'coding') -> list:
    """
    Runs transcription on a list of DNA sequences.
    
    Args:
        sequences (list): A list of DNA sequence strings.
        strand_type (str, optional): 'coding' or 'template'. Defaults to 'coding'.
        
    Returns:
        list: A list of the resulting RNA sequences.
    """
    return [transcribe(seq, strand_type) for seq in sequences]


# Question 5: Reverse Complement Generation

# Define the full complement map, including degenerate nucleotides
COMPLEMENT_MAP = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', # Purine <-> Pyrimidine
    'S': 'S', 'W': 'W', # Strong <-> Weak (self-complementary)
    'K': 'M', 'M': 'K', # Keto <-> Amino
    'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D',
    'N': 'N', # Any
    # Add lowercase just in case, though we primarily uppercase
    'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
}

# Add all other degenerate complements
COMPLEMENT_MAP.update({k.lower(): v.lower() for k, v in COMPLEMENT_MAP.items()})


def reverse_complement(seq: str) -> str:
    """
    Generates the reverse complement of a DNA sequence.
    
    Handles standard (ATGC) and degenerate (R, Y, S, etc.) nucleotides.
    Assumes input is 5'-3' and returns 5'-3' reverse complement.
    
    Args:
        seq (str): The 5'-3' DNA sequence.
        
    Returns:
        str: The 5'-3' reverse complement sequence.
    """
    # An efficient way using str.translate
    # 1. Create the translation table
    # We must use all 256 byte values for maketrans
    
    # Simple, human-readable loop approach (more student-like)
    seq = seq.upper()
    complement_list = []
    
    for nuc in seq:
        if nuc in COMPLEMENT_MAP:
            complement_list.append(COMPLEMENT_MAP[nuc])
        else:
            # Handle any characters not in our map (like '-')
            complement_list.append(nuc) 
            
    # 1. Get complement
    complement_seq = "".join(complement_list)
    
    # 2. Reverse it
    reverse_comp_seq = complement_seq[::-1]
    
    return reverse_comp_seq

def reverse_complement_optimized(seq: str) -> str:
    """
    A more optimized version of reverse_complement using str.translate.
    (For discussion in Part C)
    
    Args:
        seq (str): The 5'-3' DNA sequence.
        
    Returns:
        str: The 5'-3' reverse complement sequence.
    """
    seq = seq.upper()
    
    # Define the characters to be replaced and their replacements
    in_chars = "ATGCYRSwKMBDHVN"
    out_chars = "TACGRYSWMKVHDBN"
    
    # Create the translation table
    translation_table = str.maketrans(in_chars, out_chars)
    
    # 1. Get complement using the translation table
    complement_seq = seq.translate(translation_table)
    
    # 2. Reverse the string
    return complement_seq[::-1]