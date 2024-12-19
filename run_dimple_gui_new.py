import sys
import ast
import re

from DIMPLE.DIMPLE import (
    print_all,
    post_qc,
    addgene,
    DIMPLE,
    generate_DMS_fragments,
)
from DIMPLE.utilities import parse_custom_mutations

from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QWidget,
    QLabel,
    QLineEdit,
    QComboBox,
    QPushButton,
    QRadioButton,
    QCheckBox,
    QTextEdit,
    QFileDialog,
)
from PyQt5.QtGui import QIcon
from Bio import SeqIO
from Bio.Seq import Seq

import logging
from datetime import datetime

# Set up logging
logger = logging.getLogger(__name__)
log_file = "logs/DIMPLE-{:%Y-%m-%d-%s}.log".format(datetime.now())
logger.basicConfig = logging.basicConfig(filename=log_file, level=logging.INFO)

logger.info("Started")

# Set up lists and dicts
amino_acids = [
    "Cys",
    "Asp",
    "Ser",
    "Gln",
    "Met",
    "Asn",
    "Pro",
    "Lys",
    "Thr",
    "Phe",
    "Ala",
    "Gly",
    "Ile",
    "Leu",
    "His",
    "Arg",
    "Trp",
    "Val",
    "Glu",
    "Tyr",
]

e_coli_usage = {
    "TTT": 0.58,
    "TTC": 0.42,
    "TTA": 0.14,
    "TTG": 0.13,
    "TAT": 0.59,
    "TAC": 0.41,
    "TAA": 0.61,
    "TAG": 0.09,
    "CTT": 0.12,
    "CTC": 0.1,
    "CTA": 0.04,
    "CTG": 0.47,
    "CAT": 0.57,
    "CAC": 0.43,
    "CAA": 0.34,
    "CAG": 0.66,
    "ATT": 0.49,
    "ATC": 0.39,
    "ATA": 0.11,
    "ATG": 1,
    "AAT": 0.49,
    "AAC": 0.51,
    "AAA": 0.74,
    "AAG": 0.26,
    "GTT": 0.28,
    "GTC": 0.2,
    "GTA": 0.17,
    "GTG": 0.35,
    "GAT": 0.63,
    "GAC": 0.37,
    "GAA": 0.68,
    "GAG": 0.32,
    "TCT": 0.17,
    "TCC": 0.15,
    "TCA": 0.14,
    "TCG": 0.14,
    "TGT": 0.46,
    "TGC": 0.54,
    "TGA": 0.3,
    "TGG": 1,
    "CCT": 0.18,
    "CCC": 0.13,
    "CCA": 0.2,
    "CCG": 0.49,
    "CGT": 0.36,
    "CGC": 0.36,
    "CGA": 0.07,
    "CGG": 0.11,
    "ACT": 0.19,
    "ACC": 0.4,
    "ACA": 0.17,
    "ACG": 0.25,
    "AGT": 0.16,
    "AGC": 0.25,
    "AGA": 0.07,
    "AGG": 0.04,
    "GCT": 0.18,
    "GCC": 0.26,
    "GCA": 0.23,
    "GCG": 0.33,
    "GGT": 0.35,
    "GGC": 0.37,
    "GGA": 0.13,
    "GGG": 0.15,
}

human_usage = {
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


class DimpleApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DIMPLE Mutational Library Design")
        self.setGeometry(100, 100, 800, 600)

        # Central widget
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout()
        self.central_widget.setLayout(self.layout)

        # Input fields
        self.create_inputs()

        # Output text area
        self.create_output_area()

    def create_inputs(self):
        # General Inputs
        general_layout = QVBoxLayout()

        # File and directory selections

        # Working directory / output directory
        output_layout = QHBoxLayout()
        self.output_dir_label = QLabel("Output directory: not selected")
        self.output_dir_button = QPushButton("Select output directory")
        self.output_dir_button.clicked.connect(self.select_output_directory)
        output_layout.addWidget(self.output_dir_label)
        output_layout.addWidget(self.output_dir_button)
        general_layout.addLayout(output_layout)

        # Target gene file
        input_layout = QHBoxLayout()
        self.input_file_label = QLabel("Target gene file: not selected")
        self.input_file_button = QPushButton("Select target gene file")
        self.input_file_button.clicked.connect(self.select_input_file)
        input_layout.addWidget(self.input_file_label)
        input_layout.addWidget(self.input_file_button)
        general_layout.addLayout(input_layout)

        # Oligo and generation options
        oligo_layout = QGridLayout()
        self.oligo_label = QLabel("Oligo length:")
        self.oligo_input = QLineEdit("250")

        self.fragment_label = QLabel("Fragment length:")
        self.fragment_input = QLineEdit("auto")

        self.barcode_label = QLabel("Barcode start:")
        self.barcode_input = QLineEdit("0")

        self.fragment_label = QLabel("Fragment length:")
        self.fragment_input = QLineEdit("auto")

        self.overlap_label = QLabel("Fragment overlap:")
        self.overlap_input = QLineEdit("4")

        # Melting temperature options
        self.melting_label = QLabel("Gene primer Tm:")
        self.melting_high_label = QLabel("High:")
        self.melting_low_input = QLineEdit("58")
        self.melting_low_label = QLabel("Low:")
        self.melting_high_input = QLineEdit("62")

        melting_temp_layout = QHBoxLayout()
        melting_temp_layout.addWidget(self.melting_low_label)
        melting_temp_layout.addWidget(self.melting_low_input)
        melting_temp_layout.addWidget(self.melting_high_label)
        melting_temp_layout.addWidget(self.melting_high_input)

        oligo_layout.addWidget(self.oligo_label, 0, 0)
        oligo_layout.addWidget(self.oligo_input, 0, 1)

        oligo_layout.addWidget(self.barcode_label, 1, 0)
        oligo_layout.addWidget(self.barcode_input, 1, 1)

        oligo_layout.addWidget(self.overlap_label, 0, 2)
        oligo_layout.addWidget(self.overlap_input, 0, 3)

        oligo_layout.addWidget(self.melting_label, 1, 2)
        oligo_layout.addLayout(melting_temp_layout, 1, 3)

        general_layout.addLayout(oligo_layout)
        self.layout.addLayout(general_layout)

        # Restriction site options
        restriction_layout = QGridLayout()

        self.restriction_label = QLabel("Type IIS restriction sequence:")
        self.restriction_input = QLineEdit("GGTCTC(G)1/5")
        self.restriction_label.setToolTip(
            "Specify the restriction site and overhang in the format 'GGTCTC(G)1/5'"
        )
        restriction_layout.addWidget(self.restriction_label, 0, 0)
        restriction_layout.addWidget(self.restriction_input, 1, 0)

        self.enzyme_label = QLabel("Enzyme:")
        self.enzyme_selection_box = QComboBox()
        enzyme_list = ["---", "BsaI", "BsmBI", "BbsI", "PaqCI"]
        self.enzyme_selection_box.addItems(enzyme_list)

        restriction_layout.addWidget(self.enzyme_label, 0, 1)
        restriction_layout.addWidget(self.enzyme_selection_box, 1, 1)

        self.avoid_label = QLabel("Avoid sites:")
        self.avoid_input = QLineEdit("GGTCTC")
        self.avoid_label.setToolTip("Specify sites to avoid in the library design")
        restriction_layout.addWidget(self.avoid_label, 0, 2)
        restriction_layout.addWidget(self.avoid_input, 1, 2)

        # Detect if enzyme is selected and then update the restriction site field.
        # Also add the selected enzyme sequence to the avoid list.
        enzyme_sequence_dict = {
            "---": ("", ""),
            "BsaI": ("GGTCTC(G)1/5", "GGTCTC"),
            "BsmBI": ("CGTCTC(G)1/5", "CGTCTC"),
            "BbsI": ("GAAGAC(GG)2/6", "GAAGAC"),
            "PaqCI": ("CACCTGC(GGGG)4/8", "CACCTGC"),
        }

        # Update the sequence field.

        self.enzyme_selection_box.currentIndexChanged.connect(
            lambda: self.restriction_input.setText(
                enzyme_sequence_dict[self.enzyme_selection_box.currentText()][0]
            )
        )

        # If the sequence field is changed (except if it is due to an enzyme being selected),
        # reset the enzyme selection box to default.
        self.restriction_input.textChanged.connect(
            lambda: (
                self.enzyme_selection_box.setCurrentIndex(0)
                if self.restriction_input.text()
                not in [
                    enzyme_sequence_dict[enzyme][0]
                    for enzyme in enzyme_sequence_dict.keys()
                ]
                else None
            )
        )

        # Add to the the avoid field. If the avoid sequence is already present, do not add it again.

        self.enzyme_selection_box.currentIndexChanged.connect(
            lambda: self.avoid_input.setText(
                self.avoid_input.text()
                + ","
                + enzyme_sequence_dict[self.enzyme_selection_box.currentText()][1]
                if enzyme_sequence_dict[self.enzyme_selection_box.currentText()][1]
                not in self.avoid_input.text()
                else self.avoid_input.text()
            )
        )

        self.layout.addLayout(restriction_layout)

        # Mutation Options
        mutation_layout = QGridLayout()
        self.design_label = QLabel("Mutation options:")
        self.include_sub = QCheckBox("Include substitutions: Deep mutational scan")
        self.include_ins = QCheckBox("Include insertions (in nucleotides)")
        self.include_del = QCheckBox("Include deletions (in nucleotides)")
        self.include_dis = QCheckBox("Include domain insertions")

        self.include_stop = QCheckBox("Include stop codons")
        self.include_synonymous = QCheckBox("Include synonymous mutations")

        # Mutation field expansions
        self.sub_field = QLineEdit(
            "Cys,Asp,Ser,Gln,Met,Asn,Pro,Lys,Thr,Phe,Ala,Gly,Ile,Leu,His,Arg,Trp,Val,Glu,Tyr"
        )
        self.ins_field = QLineEdit("GGC,GGCTCT,GGCTCTGGA")
        self.del_field = QLineEdit("3,6,9")
        self.dis_field = QLineEdit("AGCGGGAGACCGGGGTCTCTGAGC")

        # Additional mutation options
        self.double_fragments = QCheckBox("Double fragments per oligo")
        self.max_nucleotide_diff = QCheckBox("Maximize nucleotide differences")
        self.double_mutations = QCheckBox("Make double mutations")

        self.sub_field.hide()
        self.ins_field.hide()
        self.del_field.hide()
        self.dis_field.hide()

        self.include_sub.stateChanged.connect(
            lambda: self.toggle_field(self.include_sub, self.sub_field)
        )
        self.include_ins.stateChanged.connect(
            lambda: self.toggle_field(self.include_ins, self.ins_field)
        )
        self.include_del.stateChanged.connect(
            lambda: self.toggle_field(self.include_del, self.del_field)
        )
        self.include_dis.stateChanged.connect(
            lambda: self.toggle_field(self.include_dis, self.dis_field)
        )

        mutation_layout.addWidget(self.design_label, 0, 0)
        mutation_layout.addWidget(self.include_sub, 1, 0)
        mutation_layout.addWidget(self.sub_field, 1, 1)
        mutation_layout.addWidget(self.include_ins, 2, 0)
        mutation_layout.addWidget(self.ins_field, 2, 1)
        mutation_layout.addWidget(self.include_del, 3, 0)
        mutation_layout.addWidget(self.del_field, 3, 1)

        mutation_layout.addWidget(self.include_stop, 4, 0)
        mutation_layout.addWidget(self.include_synonymous, 5, 0)
        mutation_layout.addWidget(self.double_fragments, 6, 0)
        mutation_layout.addWidget(self.max_nucleotide_diff, 7, 0)

        self.layout.addLayout(mutation_layout)

        # Select codon usage
        codon_layout = QGridLayout()
        self.codon_label = QLabel("Select codon usage:")
        # Radio box for codon usage
        self.e_coli_usage = QRadioButton("E. coli")
        self.human_usage = QRadioButton("Human")
        self.human_usage.setChecked(True)
        self.custom_usage = QRadioButton("Custom codon usage")

        # Custom codon file selection
        self.custom_codon_button = QPushButton("Select custom codon table")
        self.custom_codon_button.clicked.connect(self.select_custom_codon_usage)

        # Avoid breaksites in custom mutations
        # Reveal avoid field if checkbox is selected
        self.avoid_checkbox = QCheckBox("Avoid breaksites in custom mutations")
        self.avoid_field = QLineEdit("Comma separated AA positions")

        self.avoid_field.hide()

        self.avoid_checkbox.stateChanged.connect(
            lambda: self.toggle_field(self.avoid_checkbox, self.avoid_field)
        )

        codon_layout.addWidget(self.codon_label, 0, 0)
        codon_layout.addWidget(self.e_coli_usage, 1, 0)
        codon_layout.addWidget(self.human_usage, 2, 0)
        codon_layout.addWidget(self.custom_usage, 3, 0)
        codon_layout.addWidget(self.custom_codon_button, 3, 1)
        codon_layout.addWidget(self.avoid_checkbox, 4, 0)
        codon_layout.addWidget(self.avoid_field, 4, 1)

        self.layout.addLayout(codon_layout)

        # Other options
        other_options_layout = QHBoxLayout()
        self.random_seed_label = QLabel("Random seed")
        self.random_seed = QLineEdit("None")
        self.random_seed.setToolTip("Set a random seed for reproducibility")

        other_options_layout.addWidget(self.random_seed_label)
        other_options_layout.addWidget(self.random_seed)

        self.layout.addLayout(other_options_layout)

        # Run Button
        self.run_button = QPushButton("Run DIMPLE")
        self.run_button.clicked.connect(self.run_dimple)
        self.layout.addWidget(self.run_button)

    def select_input_file(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self, "Select Input File", "", "All Files (*.*)"
        )
        if file_name:
            self.input_file = file_name
            self.input_file_label.setText(f"Input File: {file_name}")

            self.gene_list = self.check_fasta(file_name)

    def select_output_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_dir = directory
            self.output_dir_label.setText(f"Output directory: {directory}")
            self.output_area.append(f"Output directory selected: {directory}")

    def select_custom_codon_usage(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Custom codon usage")
        if file_name:
            self.custom_codon_file = file_name
            self.custom_codon_label.setText(f"Custom codon table: {file_name}")
            self.output_area.append(f"Custom codon table selected: {file_name}")

    def create_output_area(self):
        # Output Area
        self.output_area = QTextEdit()
        self.output_area.setReadOnly(True)
        self.layout.addWidget(self.output_area)

    def toggle_field(self, checkbox, field):
        if checkbox.isChecked():
            field.show()
        else:
            field.hide()

    def validate_inputs(self):
        errors = []

        try:
            oligo_len = int(self.oligo_input.text())
            if oligo_len <= 0:
                errors.append("Oligo Length must be a positive integer.")
        except ValueError:
            errors.append("Oligo Length must be a valid integer.")

        try:
            barcode_start = int(self.barcode_input.text())
            if barcode_start < 0:
                errors.append("Barcode start must be a non-negative integer.")
            if barcode_start >= 3000:
                errors.append(
                    "Barcode start must be less than total number of barcodes (3000)."
                )
        except ValueError:
            errors.append("Barcode Start must be a valid integer.")

        if self.include_sub.isChecked():
            if not self.sub_field.text().strip():
                errors.append("Substitutions field cannot be empty when selected.")
            if any(
                aa not in amino_acids for aa in self.sub_field.text().strip().split(",")
            ):
                errors.append("Invalid amino acid substitutions specified.")

        if self.include_ins.isChecked():
            if not self.ins_field.text().strip():
                errors.append("Insertions must be specified.")
            if any(
                len(ins) % 3 != 0 for ins in self.ins_field.text().strip().split(",")
            ):
                self.output_area.append(
                    "Warning: selected insertions will not maintain reading frame!"
                )
            if any(
                base not in ["A", "C", "G", "T"]
                for ins in self.ins_field.text().strip().split(",")
                for base in ins
            ):
                errors.append("Non-nucleotide characters found in insertions.")

        if self.include_del.isChecked():
            try:
                deletions = [int(x) for x in self.del_field.text().strip().split(",")]
                if any(d % 3 != 0 for d in deletions):
                    self.output_area.append(
                        "Warning: selected deletions will not maintain reading frame!"
                    )
                    self.output_area.append(
                        "Deletion lengths are in nucleotides, not codons!"
                    )
                if any(d < 0 for d in deletions):
                    errors.append("Deletions must be non-negative integers.")
                # Warn if deletions are not codons
                if any(d > oligo_len for d in deletions):
                    errors.append("Deletions must be less than the oligo length.")
            except ValueError:
                errors.append("Deletions must be a comma-separated list of integers.")

        if self.include_dis.isChecked():
            if not self.dis_field.text().strip():
                errors.append("Domain insertions must be specified.")
            if any(
                base not in ["A", "C", "G", "T"]
                for ins in self.dis_field.text().strip().split(",")
                for base in ins
            ):
                errors.append("Non-nucleotide characters found in domain insertions.")

        if self.restriction_input.text().strip():
            try:
                self.parse_restriction_sequence(self.restriction_input.text())
            except ValueError:
                errors.append("Invalid restriction sequence specified.")

        if self.custom_usage.isChecked():
            if not hasattr(self, "custom_codon_file"):
                errors.append("Custom codon usage file not selected.")

        if self.avoid_checkbox.isChecked():
            if not self.avoid_field.text().strip():
                errors.append("Avoid field cannot be empty when selected.")
            if any(
                pos is not int for pos in self.avoid_field.text().strip().split(",")
            ):
                errors.append("Breaksites must be integers.")

        if self.random_seed.text().strip():
            if self.random_seed.text() != "None":
                try:
                    seed = int(self.random_seed.text())
                except ValueError:
                    errors.append("Random seed must be a valid integer.")

        return errors

    def check_fasta(self, file_name):
        # Read the selected fasta file. Check that it can be opened and parsed.
        # Print the number of genes found in the file and the name of each gene.
        # If start and end are specified, print them.
        # Check whether the start and end positions are a proper ORF.

        with open(file_name, "r") as file:
            try:
                gene_list = list(SeqIO.parse(file, "fasta"))
            except UnicodeDecodeError:
                self.output_area.append(f"Error reading file: {file_name}")
                return
            n_genes = len(gene_list)

            self.output_area.append(f"Found {n_genes} genes in {file_name}")

            for gene in gene_list:
                self.output_area.append(f"Target gene: {gene.id}")

                if "start:" in gene.description and "end:" in gene.description:
                    start = int(gene.description.split("start:")[1].split(" ")[0]) - 1
                    end = int(gene.description.split("end:")[1].split(" ")[0])
                    orf_length = end - start

                    self.output_area.append(f"Start: {start} and end: {end}")
                    if orf_length % 3 != 0:
                        self.output_area.append(
                            "Warning: ORF length is not a multiple of 3."
                        )
                else:
                    self.output_area.append("No start and end positions specified.")

        return gene_list

    def parse_restriction_sequence(self, restriction_sequence):
        # Parse the restriction sequence and return the cutsite, buffer, and return
        # in the format (cutsite, buffer, overhang)

        if re.match(r"[ACGT]+\([ACGT]\)\d+/\d+", restriction_sequence):

            tmp_cutsite = restriction_sequence.split("(")
            cutsite = Seq(tmp_cutsite[0])
            cutsite_buffer = Seq(tmp_cutsite[1].split(")")[0])
            tmp_overhang = tmp_cutsite[1].split(")")[1].split("/")
            cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])

            return cutsite, cutsite_buffer, cutsite_overhang
        else:
            return ValueError(
                f"Restriction sequence {restriction_sequence} not recognized. Please check input."
            )

    def run_dimple(self):
        # Validate inputs
        errors = self.validate_inputs()
        if not hasattr(self, "input_file"):
            errors.append("Input file not selected.")
        if not hasattr(self, "output_dir"):
            errors.append("Output directory not selected.")
        if errors:
            self.output_area.append("Validation errors:")
            for error in errors:
                self.output_area.append(f"- {error}")
            return

        # Run DIMPLE
        self.output_area.append("Running DIMPLE...")

        # Set up parameters
        wDir = self.output_dir
        geneFile = self.input_file

        oligoLen = int(self.oligo_input.text())
        overlap = int(self.overlap_input.text())
        barcode_start = int(self.barcode_input.text())
        fragmentLen = (
            int(self.fragment_input.text())
            if self.fragment_input.text() != "auto"
            else None
        )
        melting_low = int(self.melting_low_input.text())
        melting_high = int(self.melting_high_input.text())

        restriction_sequence = self.restriction_input.text()
        enzyme = (
            self.enzyme_selection_box.currentText()
            if self.enzyme_selection_box.currentText() != "---"
            else None
        )

        dms = self.include_sub.isChecked() if self.include_sub.isChecked() else None
        include_stop_codons = (
            self.include_stop.isChecked() if self.include_stop.isChecked() else None
        )
        include_synonymous = (
            self.include_synonymous.isChecked()
            if self.include_synonymous.isChecked()
            else None
        )
        insertions = (
            self.ins_field.text().split(",") if self.include_ins.isChecked() else None
        )
        deletions = [int(x) for x in self.del_field.text().strip().split(",")]
        handle = self.dis_field.test() if self.include_dis.isChecked() else None

        make_double = (
            self.double_fragments.isChecked()
            if self.double_fragments.isChecked()
            else None
        )
        maximize_nucleotide_change = (
            self.max_nucleotide_diff.isChecked()
            if self.max_nucleotide_diff.isChecked()
            else None
        )

        usage = (
            "human"
            if self.human_usage.isChecked()
            else "ecoli" if self.e_coli_usage.isChecked() else "custom"
        )

        avoid_sequence = (
            self.avoid_input.text().split(",")
            if self.avoid_input.text().strip()
            else None
        )
        custom_mutations = (
            self.custom_codon_file if hasattr(self, "custom_codon_file") else None
        )

        if self.random_seed.text().strip() and self.random_seed.text() != "None":
            seed = int(self.random_seed.text())
        else:
            seed = None

        # Set up DIMPLE parameters
        DIMPLE.handle = handle
        DIMPLE.synth_len = oligoLen
        if fragmentLen:
            DIMPLE.maxfrag = fragmentLen
        else:
            DIMPLE.maxfrag = (
                oligoLen - 64 - overlap
            )  # 64 allows for cutsites and barcodes

        DIMPLE.cutsite, DIMPLE.cutsite_buffer, DIMPLE.cutsite_overhang = (
            self.parse_restriction_sequence(restriction_sequence)
        )

        #  adjust primer primerBuffer
        DIMPLE.primerBuffer += overlap

        DIMPLE.gene_primerTm = (melting_low, melting_high)
        DIMPLE.enzyme = enzyme

        DIMPLE.avoid_sequence = avoid_sequence
        DIMPLE.barcodeF = DIMPLE.barcodeF[int(barcode_start) :]
        DIMPLE.barcodeR = DIMPLE.barcodeR[int(barcode_start) :]

        DIMPLE.avoid_sequence = [Seq(x) for x in avoid_sequence]

        # Check whether restriction sequence is included in the avoid list
        if DIMPLE.cutsite not in DIMPLE.avoid_sequence:
            DIMPLE.avoid_sequence.append(DIMPLE.cutsite)
            self.output_area.append(
                f"Restriction sequence {DIMPLE.cutsite} was not included in the avoid list. Adding before continuing."
            )
            logger.warning(
                f"Restriction sequence {DIMPLE.cutsite} was not included in the avoid list. Adding before continuing."
            )

        # Set up DMS parameters
        DIMPLE.dms = dms
        DIMPLE.stop_codon = include_stop_codons
        DIMPLE.make_double = make_double
        DIMPLE.maximize_nucleotide_change = maximize_nucleotide_change

        if custom_mutations:
            # load file with custom mutations
            with open(custom_mutations) as f:
                custom_mutations = f.readlines()
            # parse custom mutations
            custom_mutations = parse_custom_mutations(custom_mutations)
        else:
            custom_mutations = None

        # Set up random seed
        DIMPLE.random_seed = seed

        if usage == "ecoli":
            DIMPLE.usage = e_coli_usage
        elif usage == "human":
            DIMPLE.usage = human_usage
        else:
            with open(usage) as f:
                usage = f.readlines()
            DIMPLE.usage = ast.literal_eval(usage.strip("\n"))

        OLS = addgene(geneFile)

        # Not including match sequence option!

        # Run DIMPLE
        try:

            generate_DMS_fragments(
                OLS,
                overlap,
                overlap,
                include_synonymous,
                custom_mutations,
                dms,
                insertions,
                deletions,
                handle,
                wDir,
            )

            post_qc(OLS)
            print_all(OLS, wDir)

        except Exception as e:
            self.output_area.append(f"Error running DIMPLE: {e}")
            return

        self.output_area.append("DIMPLE run complete.")
        self.output_area.append(f"Output files saved to {wDir}")
        self.output_area.append(f"Log file saved to {log_file}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon("resources/dimple_icon.png"))

    window = DimpleApp()
    window.show()

    sys.exit(app.exec())
