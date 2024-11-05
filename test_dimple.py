import unittest
import os
from Bio.Seq import Seq

from DIMPLE.DIMPLE import (
    print_all,
    post_qc,
    addgene,
    DIMPLE,
    generate_DMS_fragments,
)
from DIMPLE.utilities import parse_custom_mutations


class TestDIMPLE(unittest.TestCase):
    def setUp(self) -> None:
        # Set up test.
        return super().setUp()

    def tearDown(self) -> None:
        # Clean up test.
        # Delete files created during test.
        os.remove("tests/test_dimple/Kir_DMS_Gene_Primers.fasta")
        os.remove("tests/test_dimple/Kir_DMS_Oligo_Primers.fasta")
        os.remove("tests/test_dimple/Kir_DMS_Oligos.fasta")
        os.remove("tests/test_dimple/Kir_mutations.csv")

        return super().tearDown()

    def test_generate_oligos(self):
        # Test the generate_oligos function
        # Set parameters for test

        wDir = "tests/test_dimple"
        geneFile = "Kir.fa"

        DIMPLE.handle = ""
        dis = False
        matchSequences = "nomatch"
        DIMPLE.overlap = 3
        DIMPLE.synth_len = 230
        DIMPLE.maxfrag = DIMPLE.synth_len - 62 - DIMPLE.overlap
        DIMPLE.primerBuffer += DIMPLE.overlap

        DIMPLE.DMS = True
        custom_mutations = None
        # Human usage table
        DIMPLE.usage = {
            "TTT": 0.45,
            "TTC": 0.55,
            "TTA": 0.07,
            "TTG": 0.13,
            "TAT": 0.43,
            "TAC": 0.57,
            "TAA": 0.28,
            "TAG": 0.2,
            "CTT": 0.13,
            "CTC": 0.2,
            "CTA": 0.07,
            "CTG": 0.41,
            "CAT": 0.41,
            "CAC": 0.59,
            "CAA": 0.25,
            "CAG": 0.75,
            "ATT": 0.36,
            "ATC": 0.48,
            "ATA": 0.16,
            "ATG": 1,
            "AAT": 0.46,
            "AAC": 0.54,
            "AAA": 0.42,
            "AAG": 0.58,
            "GTT": 0.18,
            "GTC": 0.24,
            "GTA": 0.11,
            "GTG": 0.47,
            "GAT": 0.46,
            "GAC": 0.54,
            "GAA": 0.42,
            "GAG": 0.58,
            "TCT": 0.18,
            "TCC": 0.22,
            "TCA": 0.15,
            "TCG": 0.06,
            "TGT": 0.45,
            "TGC": 0.55,
            "TGA": 0.52,
            "TGG": 1,
            "CCT": 0.28,
            "CCC": 0.33,
            "CCA": 0.27,
            "CCG": 0.11,
            "CGT": 0.08,
            "CGC": 0.19,
            "CGA": 0.11,
            "CGG": 0.21,
            "ACT": 0.24,
            "ACC": 0.36,
            "ACA": 0.28,
            "ACG": 0.12,
            "AGT": 0.15,
            "AGC": 0.24,
            "AGA": 0.2,
            "AGG": 0.2,
            "GCT": 0.26,
            "GCC": 0.4,
            "GCA": 0.23,
            "GCG": 0.11,
            "GGT": 0.16,
            "GGC": 0.34,
            "GGA": 0.25,
            "GGG": 0.25,
        }

        insertions = False
        deletions = False
        # include_substitutions = True
        DIMPLE.stop_codon = True
        include_synonymous = True
        DIMPLE.make_double = False
        DIMPLE.maximize_nucleotide_change = False

        barcode_start = 0
        DIMPLE.barcodeF = DIMPLE.barcodeF[int(barcode_start) :]
        DIMPLE.barcodeR = DIMPLE.barcodeR[int(barcode_start) :]

        restriction_sequence = "CGTCTC(G)1/5"
        tmp_cutsite = restriction_sequence.split("(")
        DIMPLE.cutsite = Seq(tmp_cutsite[0])
        DIMPLE.cutsite_buffer = Seq(tmp_cutsite[1].split(")")[0])
        tmp_overhang = tmp_cutsite[1].split(")")[1].split("/")
        DIMPLE.cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])

        avoid_sequence = ["CGTCTC", "GGTCTC"]
        DIMPLE.avoid_sequence = [Seq(x) for x in avoid_sequence]

        print(os.path.join(wDir, geneFile).strip())

        OLS = addgene(os.path.join(wDir, geneFile).strip())

        # Call function. Writes outputs to file.

        generate_DMS_fragments(
            OLS,
            DIMPLE.overlap,
            DIMPLE.overlap,
            include_synonymous,
            custom_mutations,
            DIMPLE.DMS,
            insertions,
            deletions,
            dis,
            wDir,
        )

        # Check the output. Expected output is in tests/test_dimple/expected and is compared to the output in tests/test_dimple
        # Files:
        # Kir_DMS_Gene_Primers.fasta
        # Kir_DMS_Oligo_Primers.fasta
        # Kir_DMS_Oligos.fasta
        # Kir_mutations.csv

        with open("tests/test_dimple/expected/Kir_DMS_Gene_Primers.fasta", "r") as f:
            expected_output_gene_primers = f.read()

        with open("tests/test_dimple/Kir_DMS_Gene_Primers.fasta", "r") as f:
            output_gene_primers = f.read()

        self.assertEqual(expected_output_gene_primers, output_gene_primers)

        with open("tests/test_dimple/expected/Kir_DMS_Oligo_Primers.fasta", "r") as f:
            expected_output_oligo_primers = f.read()

        with open("tests/test_dimple/Kir_DMS_Oligo_Primers.fasta", "r") as f:
            output_oligo_primers = f.read()

        self.assertEqual(expected_output_oligo_primers, output_oligo_primers)


if __name__ == "__main__":
    unittest.main()
