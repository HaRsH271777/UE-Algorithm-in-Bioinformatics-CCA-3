"""
UE - Algorithm in Bioinformatics
CCA3 

Test File
Name: Harshvardhan Salunke
TY B.TECH CSE Ai & DS
PRN - 1032230341

"""

import unittest
import dna_toolkit as dnat # Import the main script

class TestDnaTools(unittest.TestCase):
    
    # --- Part A Tests ---
    
    def test_q1_dna_class_valid(self):
        """Test the DNA class with a valid sequence."""
        seq_str = "ATGCGGTA"
        dna = dnat.DNA(seq_str)
        self.assertEqual(str(dna), seq_str)
        self.assertEqual(len(dna), 8)

    def test_q1_dna_class_invalid(self):
        """Test the DNA class validation by expecting a failure."""
        with self.assertRaises(ValueError):
            dnat.DNA("ATGCX") # X is not a valid nucleotide

    def test_q1_dna_class_counting(self):
        """Test the nucleotide_count method of the DNA class."""
        dna = dnat.DNA("AATTCG")
        counts = dna.nucleotide_count()
        self.assertEqual(counts['A'], 2)
        self.assertEqual(counts['T'], 2)
        self.assertEqual(counts['C'], 1)
        self.assertEqual(counts['G'], 1)

    def test_q2_frequencies(self):
        """Test the standalone calculate_frequencies function."""
        freq = dnat.calculate_frequencies("AATT")
        self.assertAlmostEqual(freq['A'], 50.0)
        self.assertAlmostEqual(freq['T'], 50.0)
        self.assertAlmostEqual(freq['G'], 0.0)
        self.assertAlmostEqual(freq['C'], 0.0)

    def test_q2_frequencies_empty(self):
        """Test frequency calculation on an empty sequence."""
        freq = dnat.calculate_frequencies("")
        self.assertAlmostEqual(freq['A'], 0.0)

    def test_q3_remove_non_nucleotides(self):
        """Test the cleaning function."""
        dirty = " A-T\nG-C-X-Y-R"
        # Test with standard ATGC
        cleaned = dnat.remove_non_nucleotides(dirty, valid_chars='ATGC')
        self.assertEqual(cleaned, "ATGC")
        # Test with degenerate chars allowed
        cleaned_degen = dnat.remove_non_nucleotides(dirty, valid_chars='ATGCXYR')
        self.assertEqual(cleaned_degen, "ATGCXYR")

    def test_q3_split_to_codons(self):
        """Test the codon splitting function."""
        seq_complete = "ATGCGT"
        self.assertEqual(dnat.split_into_codons(seq_complete), ['ATG', 'CGT'])
        
        seq_incomplete = "ATGCGTA"
        self.assertEqual(dnat.split_into_codons(seq_incomplete), ['ATG', 'CGT', 'A'])
        
        seq_empty = ""
        self.assertEqual(dnat.split_into_codons(seq_empty), [])
    
    # --- Part B Tests ---

    def test_q4_transcription_coding(self):
        """Test transcription of the coding strand."""
        seq = "ATGCTTGATAG"
        rna = dnat.transcribe(seq, strand_type='coding')
        self.assertEqual(rna, "AUGCUUGAUAG")

    def test_q4_transcription_template(self):
        """Test transcription of the template strand."""
        # Coding:   5'-ATGC-3'
        # Template: 3'-TACG-5'
        seq = "TACG" 
        rna = dnat.transcribe(seq, strand_type='template')
        # RNA:      5'-AUGC-3'
        self.assertEqual(rna, "AUGC")

    def test_q4_transcription_invalid(self):
        """Test transcription with an invalid strand type."""
        with self.assertRaises(ValueError):
            dnat.transcribe("ATGC", strand_type='invalid')

    def test_q5_reverse_complement_simple(self):
        """Test reverse complement with a simple sequence."""
        seq = "ATGC"
        # Complement: TACG
        # Reverse: GCAT
        self.assertEqual(dnat.reverse_complement(seq), "GCAT")

    def test_q5_reverse_complement_degenerate(self):
        """Test reverse complement with degenerate nucleotides."""
        # Seq:      5'-AART-Y-SW-3'
        # Comp:     3'-TTYR-A-SW-5' (S/W are self)
        # RevComp:  5'-WS-A-RYTT-3'
        seq = "AARTYSW"
        self.assertEqual(dnat.reverse_complement(seq), "WSARYTT")
        
    def test_q5_reverse_complement_optimized(self):
        """Test the optimized reverse complement function."""
        seq = "AARTYSW"
        # Compare the optimized and simple versions
        self.assertEqual(
            dnat.reverse_complement(seq), 
            dnat.reverse_complement_optimized(seq)
        )
        self.assertEqual(dnat.reverse_complement_optimized("ATGC"), "GCAT")

    def test_q5_reverse_complement_empty(self):
        """Test reverse complement on an empty sequence."""
        self.assertEqual(dnat.reverse_complement(""), "")
        


if __name__ == '__main__':
    unittest.main()