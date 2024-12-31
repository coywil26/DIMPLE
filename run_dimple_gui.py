import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
import logging
import ast
from io import StringIO
import warnings

from pathlib import Path

from DIMPLE.DIMPLE import (
    align_genevariation,
    print_all,
    post_qc,
    addgene,
    DIMPLE,
    generate_DMS_fragments,
    switch_fragmentsize,
)
from DIMPLE.utilities import parse_custom_mutations, codon_usage

# Set up logging
logger = logging.getLogger(__name__)

# Sample data for dropdowns
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
enzyme_list = ["---", "BsaI", "BsmBI"]


def run_dimple():
    if not any(
        [
            dimple_gui.deletions.get(),
            dimple_gui.insertions.get(),
            dimple_gui.substitutions.get(),
            dimple_gui.dis.get(),
        ]
    ):
        messagebox.showerror("Python Error", "Error: You must select a mutation type.")
        raise ValueError("You must select a mutation type")
    input_file = str(dimple_gui.input_file.get())
    output_dir = str(dimple_gui.output_dir.get())

    if output_dir is None:
        output_dir = input_file.rsplit("/", 1)[0] + "/"

    if any(
        [x not in ["A", "C", "G", "T", "a", "c", "g", "t"] for x in dimple_gui.handle.get()]
    ):
        raise ValueError("Genetic handle contains non nucleic bases")

    DIMPLE.synth_len = int(dimple_gui.oligo_length.get())
    overlapL = int(dimple_gui.fragment_overlap.get())
    overlapR = int(dimple_gui.fragment_overlap.get())
    if dimple_gui.include_del.get():
        overlapR = max([int(x) for x in dimple_gui.deletions.get().split(",")]) + overlapR - 3
    if dimple_gui.fragment_length.get() != "auto":
        DIMPLE.maxfrag = int(dimple_gui.fragment_length.get())
    else:
        DIMPLE.maxfrag = (
            int(dimple_gui.oligo_length.get()) - 64 - overlapL - overlapR
        )  # 62 allows for cutsites and barcodes

    # adjust primer primerBuffer
    DIMPLE.primerBuffer += overlapL
    if dimple_gui.codon_usage.get() == "ecoli":
        DIMPLE.usage = codon_usage("ecoli")
    elif dimple_gui.codon_usage.get() == "human":
        DIMPLE.usage = codon_usage("human")
    elif dimple_gui.codon_usage.get() == "custom":
        DIMPLE.usage = dimple_gui.custom_codons

    DIMPLE.barcodeF = DIMPLE.barcodeF[int(dimple_gui.barcode_start.get()) :]
    DIMPLE.barcodeR = DIMPLE.barcodeR[int(dimple_gui.barcode_start.get()) :]
    tmp_cutsite = dimple_gui.restriction_site.get().split("(")
    DIMPLE.cutsite = Seq(tmp_cutsite[0])
    DIMPLE.cutsite_buffer = Seq(tmp_cutsite[1].split(")")[0])
    tmp_overhang = tmp_cutsite[1].split(")")[1].split("/")
    DIMPLE.cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])
    DIMPLE.avoid_sequence = [Seq(x) for x in dimple_gui.avoid_sites.get().split(",")]

    enzyme = dimple_gui.enzyme.get()
    if enzyme == "---":
        DIMPLE.enzyme = None
    else:
        DIMPLE.enzyme = enzyme

    DIMPLE.aminoacids = dimple_gui.substitutions.get().split(",")
    DIMPLE.stop_codon = dimple_gui.include_stop.get()
    DIMPLE.dms = dimple_gui.include_sub.get()
    DIMPLE.handle = dimple_gui.handle.get()
    DIMPLE.make_double = dimple_gui.make_double.get()
    DIMPLE.doublefrag = dimple_gui.double_fragments.get()
    DIMPLE.gene_primerTm = (
        int(dimple_gui.melting_low.get()),
        int(dimple_gui.melting_high.get()),
    )
    DIMPLE.maximize_nucleotide_change = dimple_gui.maximize_nucleotide_change.get()

    OLS = addgene(input_file)
    if dimple_gui.avoid_breaksites_custom_mutations.get():
        OLS[0].problemsites = set(int(x) for x in dimple_gui.custom_mutations.keys())
        # add extras
        if dimple_gui.avoid_others_list.get() != "":
            OLS[0].problemsites.update(
                [int(x) for x in dimple_gui.avoid_others_list.get().split(",")]
            )
        for i in range(len(OLS[0].breaksites)):
            switch_fragmentsize(OLS[0], 1, OLS)

    if dimple_gui.matchSequences.get() == "match":
        align_genevariation(OLS)
    if dimple_gui.include_del.get() == 0:
        deletions = False
    else:
        deletions = [int(x) for x in dimple_gui.deletions.get().split(",")]
    if dimple_gui.include_ins.get() == 0:
        insertions = False
    else:
        insertions = dimple_gui.insertions.get().split(",")

    DIMPLE.random_seed = int(dimple_gui.random_seed.get())

    generate_DMS_fragments(
        OLS,
        overlapL,
        overlapR,
        dimple_gui.include_synonymous.get(),
        dimple_gui.custom_mutations,
        dimple_gui.include_sub.get(),
        insertions,
        deletions,
        dimple_gui.dis.get(),
        output_dir,
    )

    post_qc(OLS)
    print_all(OLS, output_dir)


class RedirectedStdout:
    """Custom class to redirect stdout to a Tkinter Text widget."""
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, message):
        self.text_widget.config(state="normal")  # Temporarily enable
        self.text_widget.insert(tk.END, message)
        self.text_widget.config(state="disabled")  # Re-disable
        self.text_widget.see(tk.END)  # Auto-scroll
        # Also send to stdout.
        sys.__stdout__.write(message)

    def flush(self):
        """Flush method to maintain compatibility with stdout."""
        pass

class DimpleApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DIMPLE Mutational Library Design")
        self.geometry("800x800")
        self.configure(padx=10, pady=10)

        # Variables to store user inputs
        self.output_dir = tk.StringVar(value="Not selected")
        self.input_file = tk.StringVar(value="Not selected")
        self.oligo_length = tk.StringVar(value="250")
        self.fragment_length = tk.StringVar(value="auto")
        self.fragment_overlap = tk.StringVar(value="4")
        self.barcode_start = tk.StringVar(value="0")
        self.melting_low = tk.StringVar(value="58")
        self.melting_high = tk.StringVar(value="62")
        self.restriction_site = tk.StringVar(value="GGTCTC(G)1/5")
        self.enzyme = tk.StringVar(value=enzyme_list[0])
        self.avoid_sites = tk.StringVar(value="GGTCTC")

        # Mutation options variables
        self.include_sub = tk.BooleanVar(value=False)
        self.include_ins = tk.BooleanVar(value=False)
        self.include_del = tk.BooleanVar(value=False)
        self.include_stop = tk.BooleanVar(value=False)
        self.include_synonymous = tk.BooleanVar(value=False)
        self.substitutions = tk.StringVar(
            value="Cys,Asp,Ser,Gln,Met,Asn,Pro,Lys,Thr,Phe,Ala,Gly,Ile,Leu,His,Arg,Trp,Val,Glu,Tyr"
        )
        self.insertions = tk.StringVar(value="GGC,GGCTCT,GGCTCTGGA")
        self.deletions = tk.StringVar(value="3,6,9")
        self.make_double = tk.BooleanVar(value=False)

        # Domain insertions excluded from GUI
        self.dis = tk.BooleanVar(value=False)
        self.handle = tk.StringVar(value="AGCGGGAGACCGGGGTCTCTGAGC")

        self.double_fragments = tk.BooleanVar(value=False)
        self.maximize_nucleotide_change = tk.BooleanVar(value=False)
        self.custom_mutations = {}
        self.avoid_breaksites_custom_mutations = tk.BooleanVar(value=False)

        self.avoid_others_list = tk.StringVar(value="")

        self.matchSequences = tk.IntVar()

        self.random_seed = tk.IntVar(value=1848)

        # Codon usage settings
        self.codon_usage = tk.StringVar(value="human")

        self.custom_codons = {
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

        # Flag to avoid circular resetting
        self.updating = False

        self.create_scrollable_frame()

        self.setup_ui()

    def create_scrollable_frame(self):
        """Create a scrollable frame."""
        # Create a canvas and a vertical scrollbar
        self.canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)

        # Configure canvas to work with the scrollbar
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")),
        )

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)

        # Pack the canvas and scrollbar
        self.canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def setup_ui(self):
        # Top-level input fields
        self.setup_general_options()

        # Enzyme options
        self.setup_enzyme_options()

        self.restriction_site.trace_add("write", self.reset_enzyme_field)

        # Mutation options
        self.setup_mutation_options()

        # Codon usage settings
        self.setup_codon_usage_settings()

        # Run button
        tk.Button(
            self.scrollable_frame, text="Run DIMPLE", command=self.start_dimple
        ).pack(pady=10)

        # Output area
        self.output_text = tk.Text(
            self.scrollable_frame, wrap="word", height=15, state="disabled"
        )
        self.output_text.pack(fill="both", expand=True, pady=10)

        # Redirect stdout to the Text widget
        sys.stdout = RedirectedStdout(self.output_text)

        # Redirect warnings to the Text widget
        self.redirect_warnings()


    def setup_general_options(self):
        # Create a labeled frame for General Options
        frame = tk.LabelFrame(
            self.scrollable_frame, text="General Options", padx=10, pady=10
        )
        frame.pack(fill="x", pady=10)

        # Input file selection
        input_sub_frame = tk.Frame(frame)
        input_sub_frame.pack(fill="x", pady=5)

        tk.Label(input_sub_frame, text="Target Gene File:").grid(
            row=0, column=0, sticky="w", padx=5, pady=2
        )
        tk.Label(input_sub_frame, textvariable=self.input_file).grid(
            row=0, column=1, sticky="w", padx=5, pady=2
        )
        tk.Button(input_sub_frame, text="Select", command=self.select_input_file).grid(
            row=0, column=2, padx=5, pady=2
        )

        # Output directory selection
        output_dir_sub_frame = tk.Frame(frame)
        output_dir_sub_frame.pack(fill="x", pady=5)

        tk.Label(output_dir_sub_frame, text="Output Directory:").grid(
            row=0, column=0, sticky="w", padx=5, pady=2
        )
        tk.Label(output_dir_sub_frame, textvariable=self.output_dir).grid(
            row=0, column=1, sticky="w", padx=5, pady=2
        )
        tk.Button(
            output_dir_sub_frame, text="Select", command=self.select_output_directory
        ).grid(row=0, column=2, padx=5, pady=2)

        # Other top-level inputs (e.g., oligo length, fragment length)
        options_sub_frame = tk.Frame(frame)
        options_sub_frame.pack(fill="x", pady=5)

        # Oligo Length
        tk.Label(options_sub_frame, text="Oligo Length:").grid(
            row=0, column=0, sticky="w", padx=5, pady=2
        )
        tk.Entry(options_sub_frame, textvariable=self.oligo_length, width=10).grid(
            row=0, column=1, padx=5, pady=2
        )

        # Fragment Length
        tk.Label(options_sub_frame, text="Fragment Length:").grid(
            row=1, column=0, sticky="w", padx=5, pady=2
        )
        tk.Entry(options_sub_frame, textvariable=self.fragment_length, width=10).grid(
            row=1, column=1, padx=5, pady=2
        )

        # Fragment Overlap
        tk.Label(options_sub_frame, text="Fragment Overlap:").grid(
            row=0, column=2, sticky="w", padx=5, pady=2
        )
        tk.Entry(options_sub_frame, textvariable=self.fragment_overlap, width=10).grid(
            row=0, column=3, padx=5, pady=2
        )

        # Barcode Start
        tk.Label(options_sub_frame, text="Barcode Start:").grid(
            row=1, column=2, sticky="w", padx=5, pady=2
        )
        tk.Entry(options_sub_frame, textvariable=self.barcode_start, width=10).grid(
            row=1, column=3, padx=5, pady=2
        )

        # Primer Tm Range
        tk.Label(options_sub_frame, text="Gene primer melting temperature range:").grid(
            row=4, column=0, sticky="w", padx=5, pady=2
        )
        tk.Entry(options_sub_frame, textvariable=self.melting_low, width=10).grid(
            row=4, column=1, padx=5, pady=2
        )
        tk.Entry(options_sub_frame, textvariable=self.melting_high, width=10).grid(
            row=4, column=2, padx=5, pady=2
        )
        tk.Checkbutton(
            options_sub_frame, text="Double fragments per oligo", variable=self.double_fragments
        ).grid(row=5, column=0, sticky="w", padx=5, pady=2)


    def setup_enzyme_options(self):
        frame = tk.LabelFrame(
            self.scrollable_frame, text="Enzyme Options", padx=10, pady=10
        )
        frame.pack(fill="x", pady=10)

        # Enzyme selection
        enzyme_frame = tk.Frame(frame)
        enzyme_frame.pack(fill="x", pady=6)

        # Restriction Site
        tk.Label(enzyme_frame, text="Restriction Sequence:").grid(
            row=0, column=0, sticky="w", padx=4, pady=2
        )
        self.restriction_site_field = tk.Entry(
            enzyme_frame, textvariable=self.restriction_site
        ).grid(row=1, column=0, sticky="w", padx=4, pady=2)

        # Enzyme
        tk.Label(enzyme_frame, text="Enzyme:").grid(
            row=0, column=1, sticky="w", padx=4, pady=2
        )
        self.enzyme_combobox = ttk.Combobox(
            enzyme_frame, values=enzyme_list, textvariable=self.enzyme
        )
        self.enzyme_combobox.grid(row=1, column=1, sticky="w", padx=4, pady=2)

        self.enzyme_combobox.bind("<<ComboboxSelected>>", self.update_enzyme_field)

        # Avoid Sites
        tk.Label(enzyme_frame, text="Avoid Sites:").grid(
            row=0, column=2, sticky="w", padx=4, pady=2
        )
        tk.Entry(enzyme_frame, textvariable=self.avoid_sites).grid(
            row=1, column=2, sticky="w", padx=4, pady=2
        )

    def setup_mutation_options(self):
        frame = tk.LabelFrame(
            self.scrollable_frame, text="Mutation Options", padx=10, pady=10
        )
        frame.pack(fill="x", pady=10)

        # Substitutions
        sub_frame = tk.Frame(frame)
        sub_frame.pack(fill="x", pady=5)
        tk.Checkbutton(
            sub_frame,
            text="Include Substitutions (DMS)",
            variable=self.include_sub,
            command=lambda: self.toggle_field(self.include_sub, self.sub_field),
        ).pack(side="left")
        self.sub_field = tk.Entry(sub_frame, textvariable=self.substitutions)
        self.sub_field.pack(side="left", padx=5)
        self.sub_field.config(state="disabled")

        # Insertions
        ins_frame = tk.Frame(frame)
        ins_frame.pack(fill="x", pady=5)
        tk.Checkbutton(
            ins_frame,
            text="Include Insertions",
            variable=self.include_ins,
            command=lambda: self.toggle_field(self.include_ins, self.ins_field),
        ).pack(side="left")
        self.ins_field = tk.Entry(ins_frame, textvariable=self.insertions)
        self.ins_field.pack(side="left", padx=5)
        self.ins_field.config(state="disabled")

        # Deletions
        del_frame = tk.Frame(frame)
        del_frame.pack(fill="x", pady=5)
        tk.Checkbutton(
            del_frame,
            text="Include Deletions",
            variable=self.include_del,
            command=lambda: self.toggle_field(self.include_del, self.del_field),
        ).pack(side="left")
        self.del_field = tk.Entry(del_frame, textvariable=self.deletions)
        self.del_field.pack(side="left", padx=5)
        self.del_field.config(state="disabled")

        # Domain insertions not included in GUI

        # Other mutation options
        tk.Checkbutton(
            frame, text="Include Stop Codons", variable=self.include_stop
        ).pack(anchor="w")
        tk.Checkbutton(
            frame, text="Include Synonymous Mutations", variable=self.include_synonymous
        ).pack(anchor="w")
        tk.Checkbutton(
            frame, text="Make Double Mutations", variable=self.make_double
        ).pack(anchor="w")
        tk.Checkbutton(
            frame,
            text="Maximize Nucleotide Changes",
            variable=self.maximize_nucleotide_change,
        ).pack(anchor="w")
        self.custom_mutations_button = tk.Button(
            frame, text="Custom Mutations", command=self.custom_mutations_input
        )
        self.custom_mutations_button.pack(pady=10)
        tk.Checkbutton(
            frame,
            text="Avoid breaksites in custom mutations",
            variable=self.avoid_breaksites_custom_mutations,
        ).pack(anchor="w")
        tk.Entry(frame, textvariable=self.avoid_others_list).pack()
        tk.Label(frame, text="Random Seed:").pack(side="left")
        tk.Entry(frame, textvariable=self.random_seed).pack(side="left")


    def setup_codon_usage_settings(self):
        frame = tk.LabelFrame(
            self.scrollable_frame, text="Codon Usage Settings", padx=10, pady=10
        )
        frame.pack(fill="x", pady=10)

        # Codon Usage
        codon_frame = tk.Frame(frame)
        codon_frame.pack(fill="x", pady=5)
        tk.Label(codon_frame, text="Codon Usage:").pack(side="left", padx=5)
        tk.Radiobutton(
            codon_frame, text="Human", value="human", variable=self.codon_usage
        ).pack(side="left", padx=5)
        tk.Radiobutton(
            codon_frame, text="E. coli", value="ecoli", variable=self.codon_usage
        ).pack(side="left", padx=5)
        tk.Radiobutton(
            codon_frame, text="Custom", value="custom", variable=self.codon_usage
        ).pack(side="left", padx=5)
        tk.Button(
            codon_frame,
            text="Custom Codon Usage",
            command=self.select_custom_codon_usage,
        ).pack(side="left", padx=5)

    def toggle_field(self, variable, field):
        if variable.get():
            field.config(state="normal")
        else:
            field.config(state="disabled")

    def select_output_directory(self):
        directory = filedialog.askdirectory()
        if directory:
            self.output_dir.set(directory)

    def select_input_file(self):
        file_name = filedialog.askopenfilename(
            filetypes=[("Fasta Files", "*.fasta"), ("All Files", "*.*")]
        )
        if file_name:
            self.input_file.set(file_name)
        if not os.path.isdir(self.output_dir.get()):
            self.output_dir.set(os.path.dirname(file_name))
        # If output dir isn't set, set it to the directory of the input file
        if not self.output_dir.get():
            self.output_dir.set(os.path.dirname(file_name))

    def select_custom_codon_usage(self):
        """Open a custom codon usage editor in a new window."""
        # Create a new Toplevel window
        newWindow = tk.Toplevel(self)
        newWindow.title("Custom Codon Usage")
        newWindow.geometry("900x400")

        # Create a Text widget for editing the custom codon dictionary
        custom_codon = tk.Text(newWindow, width=90, height=20, wrap=tk.WORD)
        custom_codon.pack(fill="both", expand=True, padx=10, pady=10)

        # Populate the Text widget with the current custom codon usage
        if not isinstance(self.codon_usage.get(), dict):
            self.codon_usage.set(self.custom_codons)

        s = StringIO()
        print(self.codon_usage.get(), file=s)
        custom_codon.insert(tk.END, s.getvalue())

        def save_changes():
            """Save the custom codon usage and close the window."""
            try:
                # Save the custom codon usage from the Text widget
                codon_dict = ast.literal_eval(custom_codon.get("1.0", "end-1c").strip())
                if isinstance(codon_dict, dict):
                    self.codon_usage.set("custom")  # Set the radio button to "Custom"
                    self.custom_codons = codon_dict  # Save the custom codons
                    self.append_output("Custom codon usage updated.")
                    newWindow.destroy()
                else:
                    self.append_output(
                        "Invalid custom codon usage format. No changes made."
                    )
            except Exception as e:
                self.append_output(f"Error parsing custom codon usage: {e}")

        def cancel_changes():
            """Discard changes and close the window."""
            self.append_output("Custom codon usage edit canceled.")
            # Reset the codon usage to default (human)
            self.codon_usage.set("human")
            newWindow.destroy()

        # Add Save and Cancel buttons at the bottom of the window
        button_frame = tk.Frame(newWindow)
        button_frame.pack(fill="x", pady=10, padx=10)

        save_button = tk.Button(
            button_frame, text="Save", command=save_changes, width=10
        )
        save_button.pack(side="left", padx=5)

        cancel_button = tk.Button(
            button_frame, text="Cancel", command=cancel_changes, width=10
        )
        cancel_button.pack(side="right", padx=5)

        # Set protocol to handle closing the window
        newWindow.protocol("WM_DELETE_WINDOW", cancel_changes)

    # Create a function to input custom mutations
    def custom_mutations_input(self):
        # Create a new window
        newWindow = tk.Toplevel(self)
        newWindow.title("Custom Mutations")
        newWindow.geometry("700x1000")

        # Create a text box for custom mutations
        custom_mutations_window = tk.Text(newWindow, height=50, width=80)
        custom_mutations_window.pack(pady=10)

        # Populate the text box with example mutations or existing mutations
        if self.custom_mutations == {}:
            custom_mutations_window.insert(
                tk.END,
                "Positions:Mutations\n1-10:All\n11-20:A,C,D,E,F,G,H,I,K,L\n33:A,C",
            )
        else:
            custom_mutations_window.insert(tk.END, "Positions:Mutations\n")
            for key, value in self.custom_mutations.items():
                custom_mutations_window.insert(tk.END, f"{key}:{value}\n")

        def save_mutations():
            """Save the custom mutations and close the window."""
            try:
                mutation_text = (
                    custom_mutations_window.get("1.0", "end-1c").strip().split("\n")
                )
                self.custom_mutations = parse_custom_mutations(mutation_text[1:])
                self.custom_mutations_button.config(
                    bg="green", activebackground="green", relief=tk.SUNKEN
                )
                self.include_substitutions.set(1)
                self.append_output("Custom mutations saved.")
            except Exception as e:
                self.append_output(f"Error parsing custom mutations: {e}")
            newWindow.destroy()

        def cancel_mutations():
            """Discard changes and close the window."""
            self.append_output("Custom mutations input canceled.")
            newWindow.destroy()

        # Add Save and Cancel buttons
        button_frame = tk.Frame(newWindow)
        button_frame.pack(fill="x", pady=10)

        save_button = tk.Button(
            button_frame, text="Save", command=save_mutations, width=10
        )
        save_button.pack(side="left", padx=5)

        cancel_button = tk.Button(
            button_frame, text="Cancel", command=cancel_mutations, width=10
        )
        cancel_button.pack(side="right", padx=5)

        # Set the close protocol to cancel changes
        newWindow.protocol("WM_DELETE_WINDOW", cancel_mutations)

    def update_enzyme_field(self, event):
        self.updating = True
        selected = self.enzyme.get()
        if selected == "BsaI":
            self.restriction_site.set("GGTCTC(G)1/5")
        elif selected == "BsmBI":
            self.restriction_site.set("CGTCTC(G)1/5")
        else:
            self.restriction_site.set("")
        # Update avoid sequence based on selected enzyme
        avoid_sequence = "GGTCTC" if selected == "BsaI" else "CGTCTC"
        # Check if avoid sequence is already in the list
        if avoid_sequence not in self.avoid_sites.get():
            self.avoid_sites.set(self.avoid_sites.get() + "," + avoid_sequence)
        self.updating = False

    def reset_enzyme_field(self, *args):
        if not self.updating:
            self.enzyme.set(enzyme_list[0])

    def append_output(self, message):
        self.output_text.config(state="normal")
        self.output_text.insert("end", message + "\n")
        self.output_text.config(state="disabled")
        self.output_text.see("end")

    def validate_inputs(self):
        errors = []
        if not os.path.isdir(self.output_dir.get()):
            errors.append("Output directory not selected.")
        if not os.path.isfile(self.input_file.get()):
            errors.append("Input file not selected or invalid.")
        try:
            int(self.oligo_length.get())
        except ValueError:
            errors.append("Oligo length must be an integer.")
        return errors

    def redirect_warnings(self):
        """Redirect warnings to the Text widget."""
        def custom_showwarning(message, category, filename, lineno, file=None, line=None):
            warning_message = f"Warning: {message} (Category: {category.__name__}, File: {filename}, Line: {lineno})\n"
            self.output_text.config(state="normal")
            self.output_text.insert(tk.END, warning_message)
            self.output_text.config(state="disabled")
            self.output_text.see(tk.END)
            # Also print to stdout
            sys.__stderr__.write(warning_message)

        warnings.showwarning = custom_showwarning


    def start_dimple(self):

        log_file = (
            Path(str(self.output_dir)) / "logs" / f"DIMPLE-{datetime.now().strftime('%Y-%m-%d-%s')}.log"
        )
        if not Path(log_file.parent).exists():
            Path(log_file.parent).mkdir(parents=True)

        logger.basicConfig = logging.basicConfig(
            filename=Path(log_file), level=logging.INFO
        )

        logger.info("Started DIMPLE run.")

        errors = self.validate_inputs()
        if errors:
            messagebox.showerror("Input Errors", "\n".join(errors))
            return

        self.append_output("Running DIMPLE...")

        try:
            run_dimple()
            self.append_output("DIMPLE run completed successfully.")
        except Exception as e:
            self.append_output(f"Error: {str(e)}")

        self.append_output(f"Output Directory: {self.output_dir.get()}")
        self.append_output(f"Input File: {self.input_file.get()}")
        self.append_output(f"Log saved to: {self.output_dir.get()}/logs")
        self.append_output("DIMPLE run complete.")


if __name__ == "__main__":
    dimple_gui = DimpleApp()
    dimple_gui.mainloop()
