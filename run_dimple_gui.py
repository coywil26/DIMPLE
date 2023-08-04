# RUN DIMPLE
# script for GUI

from DIMPLE.DIMPLE import align_genevariation, print_all, post_qc, addgene, DIMPLE, generate_DMS_fragments
from DIMPLE.utilities import parse_custom_mutations
from Bio.Seq import Seq
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import ast
from io import StringIO

def run():
    if not any([app.delete.get(), app.insert.get(), app.include_substitutions.get(), app.dis.get()]):
        messagebox.showerror('Python Error', 'Error: You must select a mutation type.')
        raise ValueError('You must select a mutation type')
    if app.wDir is None:
        app.wDir = app.geneFile.rsplit('/', 1)[0]+'/'

    #if any([x not in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] for x in app.handle.get()]):
    #   raise ValueError('Genetic handle contains non nucleic bases')

    #DIMPLE.handle = app.handle.get()
    DIMPLE.synth_len = int(app.oligoLen.get())
    overlapL = int(app.overlap.get())
    overlapR = int(app.overlap.get())
    if app.delete.get():
        overlapR = max([int(x) for x in app.deletions.get().split(',')]) + overlapR - 3
    if app.fragmentLen.get() != 'auto':
        DIMPLE.maxfrag = int(app.fragmentLen.get())
    else:
        DIMPLE.maxfrag = int(app.oligoLen.get()) - 62 - overlapL - overlapR  # 62 allows for cutsites and barcodes

    #adjust primer primerBuffer
    DIMPLE.primerBuffer += overlapL
    if app.codon_usage == 'ecoli':
        DIMPLE.usage = {
            'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'TAT': 0.59, 'TAC': 0.41, 'TAA': 0.61, 'TAG': 0.09,
            'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47, 'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
            'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1, 'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
            'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
            'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'TGT': 0.46, 'TGC': 0.54, 'TGA': 0.3, 'TGG': 1,
            'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11,
            'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04,
            'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15
        }  # E.coli codon usage table
    elif app.codon_usage == 'human':
        DIMPLE.usage = {
            'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.2,
            'CTT': 0.13, 'CTC': 0.2, 'CTA': 0.07, 'CTG': 0.41, 'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75,
            'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1, 'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58,
            'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47, 'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
            'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.52, 'TGG': 1,
            'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11, 'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21,
            'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.2, 'AGG': 0.2,
            'GCT': 0.26, 'GCC': 0.4, 'GCA': 0.23, 'GCG': 0.11, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25
        }
    else:
        DIMPLE.usage = app.codon_usage
    DIMPLE.barcodeF = DIMPLE.barcodeF[int(app.barcode_start.get()):]
    DIMPLE.barcodeR = DIMPLE.barcodeR[int(app.barcode_start.get()):]
    DIMPLE.cutsite = Seq(app.restriction_sequence.get())
    DIMPLE.avoid_sequence = [Seq(x) for x in app.avoid_sequence.get().split(',')]
    DIMPLE.aminoacids = app.substitutions.get().split(',')
    DIMPLE.stop_codon = app.stop.get()
    DIMPLE.dms = app.include_substitutions.get()
    DIMPLE.make_double = app.make_double.get()
    DIMPLE.handle = app.handle
    print(DIMPLE.usage)

    OLS = addgene(app.geneFile)

    if app.matchSequences.get() == 'match':
        align_genevariation(OLS)
    if app.delete.get() == 0:
        deletions = False
    else:
        deletions = [int(x) for x in app.deletions.get().split(',')]
    if app.insert.get() == 0:
        insertions = False
    else:
        insertions = app.insertions.get().split(',')
    generate_DMS_fragments(OLS, overlapL, overlapR, app.synonymous.get(), app.custom_mutations, app.include_substitutions.get(), insertions, deletions, app.dis.get(), app.wDir)

    post_qc(OLS)
    print_all(OLS, app.wDir)


class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.winfo_toplevel().title("DIMPLE Deep Indel Missense Programmable Library Engineering")
        self.matchSequences = tk.IntVar()
        self.mutationType = tk.IntVar()
        self.usage = tk.IntVar()
        self.include_substitutions = tk.IntVar()
        self.delete = tk.IntVar()
        self.insert = tk.IntVar()
        self.stop = tk.IntVar()
        self.synonymous = tk.IntVar()
        self.dis = tk.IntVar()
        self.make_double = tk.IntVar()
        self.custom_mutations = {}

        self.wDir_file = tk.Button(self, text='Working Directory', command=self.browse_wDir)
        self.wDir_file.pack()
        self.wDir = None

        self.gene_file = tk.Button(self, text='Target Gene File', command=self.browse_gene)
        self.gene_file.pack()
        self.geneFile = None

        tk.Label(self, text='Oligo Length').pack()
        self.oligoLen = tk.Entry(self, textvariable=tk.StringVar(self, '230'))
        self.oligoLen.pack()

        tk.Label(self, text='Fragment Length').pack()
        self.fragmentLen = tk.Entry(self, textvariable=tk.StringVar(self, 'auto'))
        self.fragmentLen.pack()

        tk.Label(self, text='Fragment Overlap (This will change if deletions are selected)').pack()
        self.overlap = tk.Entry(self, textvariable=tk.StringVar(self, '4'))
        self.overlap.pack()

        tk.Label(self, text='Barcode Start position (3000 total available)').pack()
        self.barcode_start = tk.Entry(self, textvariable=tk.StringVar(self, '0'))
        self.barcode_start.pack()

        tk.Label(self, text='Type IIS restriction sequence').pack()
        self.restriction_sequence = tk.Entry(self, textvariable=tk.StringVar(self, 'CGTCTC'))
        self.restriction_sequence.pack()

        tk.Label(self, text='Sequences to avoid').pack()
        self.avoid_sequence = tk.Entry(self, width=50, textvariable=tk.StringVar(self, 'CGTCTC, GGTCTC'))
        self.avoid_sequence.pack()

        self.stop_codon = tk.Checkbutton(self, text="Include Stop Codons", variable=self.stop)
        self.stop_codon.pack()
        self.stop_codon.select()

        self.synonymous_check = tk.Checkbutton(self, text="Include Synonymous Mutations", variable=self.synonymous)
        self.synonymous_check.pack()
        self.synonymous_check.select()

        self.codon_usage = 'human'

        self.custom_codons =  {
            'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.2,
            'CTT': 0.13, 'CTC': 0.2, 'CTA': 0.07, 'CTG': 0.41, 'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75,
            'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1, 'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58,
            'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47, 'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
            'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.52, 'TGG': 1,
            'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11, 'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21,
            'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.2, 'AGG': 0.2,
            'GCT': 0.26, 'GCC': 0.4, 'GCA': 0.23, 'GCG': 0.11, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25
        }

        def openNewWindow():
            newWindow = tk.Toplevel(self)
            newWindow.title("Custom Codon Usage")
            newWindow.geometry("900x400")
            custom_codon = tk.Text(newWindow, width=750, height=400, wrap=tk.WORD)
            custom_codon.pack()
            # if self.codon_usage is not a dictionary then set it to self.custom_codons
            if type(self.codon_usage) is not dict:
                self.codon_usage = self.custom_codons
            s = StringIO()
            print(self.codon_usage, file=s)
            custom_codon.insert(tk.END, s.getvalue())
            def on_closing():
                # set codon usage
                self.codon_usage = ast.literal_eval(custom_codon.get("1.0", 'end-1c').strip('\n'))
                newWindow.destroy()
            newWindow.protocol("WM_DELETE_WINDOW", on_closing)

        def ecoli_ON():
            self.codon_usage = 'ecoli'

        def human_ON():
            self.codon_usage = 'human'

        tk.Label(self, text='Codon Usage', font="helvetica 12 underline").pack(pady=5)
        self.ecoli_check = tk.Radiobutton(self, text="E. coli", variable=self.usage, value=1, command=ecoli_ON)
        self.ecoli_check.pack()
        self.human_check = tk.Radiobutton(self, text="Human", variable=self.usage, value=0, command=human_ON)
        self.human_check.pack()
        self.custom_codon_button = tk.Button(self, text="Custom Codon Usage", command=openNewWindow)
        self.custom_codon_button.pack()

        tk.Label(self, text='Select Mutations', font="helvetica 12 underline").pack(pady=5)

        self.include_dis = tk.Checkbutton(self, text='Domain Insertion Scan', variable=self.dis)
        self.include_dis.pack()
        self.handle = tk.Entry(self, width=50, textvariable=tk.StringVar(self, 'AGCGGGAGACCGGGGTCTCTGAGC'))
        self.handle.pack()

        self.delete_check = tk.Checkbutton(self, text="List of Deletions", variable=self.delete)
        self.delete_check.pack()
        self.delete_check.deselect()
        self.deletions = tk.Entry(self, width=50, textvariable=tk.StringVar(self, '3,6'))
        self.deletions.pack()

        self.insert_check = tk.Checkbutton(self, text="List of Insertions", variable=self.insert)
        self.insert_check.pack()
        self.insert_check.deselect()
        self.insertions = tk.Entry(self, width=50, textvariable=tk.StringVar(self, 'GGG,GGGGGG'))
        self.insertions.pack()

        self.include_sub_check = tk.Checkbutton(self, text='Deep Mutational Scan', variable=self.include_substitutions)
        self.include_sub_check.pack()
        self.include_sub_check.deselect()
        self.double_it = tk.Checkbutton(self, text='Make Double Mutations', variable=self.make_double)

        self.substitutions = tk.Entry(self, width=80, textvariable=tk.StringVar(self, "Cys,Asp,Ser,Gln,Met,Asn,Pro,Lys,Thr,Phe,Ala,Gly,Ile,Leu,His,Arg,Trp,Val,Glu,Tyr"))
        self.substitutions.pack()

        self.double_it.pack()
        self.double_it.deselect()

        self.custom_mutations_button = tk.Button(self, text="Custom Mutations", command=self.custom_mutations_input)
        self.custom_mutations_button.pack(pady=10)

        #self.matchSequences_check = tk.Checkbutton(self, text='Match Sequences', variable=self.matchSequences)
        #self.matchSequences_check.pack()

        self.run = tk.Button(self, text='Run DIMPLE', command=run).pack(pady=10)

        self.output = tk.Text(self, height=10, width=60).pack()

    def browse_wDir(self):
        self.wDir = filedialog.askdirectory(title="Select a File")
        if self.wDir != '':
            self.wDir_file.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.wDir = None

    def browse_gene(self):
        self.geneFile = filedialog.askopenfilename(title="Select a File")
        if self.geneFile != '':
            self.gene_file.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.geneFile = None

    # create a function to input custom mutations
    def custom_mutations_input(self):
        # create a new window
        newWindow = tk.Toplevel(self)
        newWindow.title("Custom Mutations          (close window when finished)")
        newWindow.geometry('700x1000')
        # create a text box for custom mutations
        custom_mutations_window = tk.Text(newWindow, height=50, width=80)
        custom_mutations_window.pack()
        # example mutations
        if self.custom_mutations == {}:
            custom_mutations_window.insert(tk.END, "Positions:Mutations\n1-10:All\n11-20:A,C,D,E,F,G,H,I,K,L\n33:A,C")
        else:
            custom_mutations_window.insert(tk.END, "Positions:Mutations\n")
            for key, value in self.custom_mutations.items():
                custom_mutations_window.insert(tk.END, str(key)+':'+value+'\n')
        # create a button to set custom mutations
        def on_closing():
            # set custom mutations
            mutation_text = custom_mutations_window.get("1.0", 'end-1c').strip().split('\n')
            self.custom_mutations = parse_custom_mutations(mutation_text[1:])
            newWindow.destroy()
        newWindow.protocol("WM_DELETE_WINDOW", on_closing)
        self.custom_mutations_button.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        self.include_substitutions.set(1)


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry('700x1000')
    app = Application(master=root)
    app.mainloop()
