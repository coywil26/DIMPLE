# RUN DIMPLE
# script for GUI

from DIMPLE.DIMPLE import align_genevariation, print_all, post_qc, addgene, DIMPLE, generate_DMS_fragments
from Bio.Seq import Seq
import tkinter as tk
from tkinter import filedialog


def run():
    if not any([app.delete.get(), app.insert.get(), app.include_substitutions.get()]):
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

    DIMPLE.barcodeF = DIMPLE.barcodeF[int(app.barcode_start.get()):]
    DIMPLE.barcodeR = DIMPLE.barcodeR[int(app.barcode_start.get()):]
    DIMPLE.cutsite = Seq(app.restriction_sequence.get())
    DIMPLE.avoid_sequence = [Seq(x) for x in app.avoid_sequence.get().split(',')]
    DIMPLE.stop_codon = app.stop.get()
    if app.mutationType.get() == 1:
        DIMPLE.dms = True
    else:
        DIMPLE.dms = False

    if app.usage:
        DIMPLE.usage = 'ecoli'
    else:
        DIMPLE.usage = 'human'

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
    generate_DMS_fragments(OLS, overlapL, overlapR, app.include_substitutions.get(), insertions, deletions, app.wDir)

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
        #self.overlap.pack()

        tk.Label(self, text='Barcode Start position').pack()
        self.barcode_start = tk.Entry(self, textvariable=tk.StringVar(self, '0'))
        self.barcode_start.pack()

        tk.Label(self, text='Type IIS restriction sequence').pack()
        self.restriction_sequence = tk.Entry(self, textvariable=tk.StringVar(self, 'CGTCTC'))
        self.restriction_sequence.pack()

        tk.Label(self, text='Sequences to avoid').pack()
        self.avoid_sequence = tk.Entry(self, width=50, textvariable=tk.StringVar(self, 'CGTCTC, GGTCTC'))
        self.avoid_sequence.pack()

        self.handle = 'AGCGGGAGACCGGGGTCTCTGAGC'

        self.stop_codon = tk.Checkbutton(self, text="Include Stop Codons", variable=self.stop)
        self.stop_codon.pack()

        def sub_ON():
            self.include_substitutions.set(1)

        # tk.Label(self, text='Type of mutations to generate', font="helvetica 12 underline").pack()
        # self.DIS_check = tk.Radiobutton(self, text="Deep Insertional Scan (Genetic Handle)", variable=self.mutationType, value=0)
        # self.DIS_check.pack()
        # self.DMS_check = tk.Radiobutton(self, text="Deep Mutational Scan", variable=self.mutationType, value=1, command=sub_ON)
        # self.DMS_check.pack()
        # self.DIS_check.select()

        tk.Label(self, text='Codon Usage', font="helvetica 12 underline").pack()
        self.ecoli_check = tk.Radiobutton(self, text="E. coli", variable=self.usage, value=1)
        self.ecoli_check.pack()
        self.human_check = tk.Radiobutton(self, text="Human", variable=self.usage, value=0)
        self.human_check.pack()

        def DMS_ON():
            if self.delete.get() or self.insert.get():
                return
                #self.DMS_check.select()

        tk.Label(self, text='Select Mutations', font="helvetica 12 underline").pack()
        self.delete_check = tk.Checkbutton(self, text="List of Deletions", variable=self.delete, command=DMS_ON)
        self.delete_check.pack()
        self.delete_check.deselect()
        self.deletions = tk.Entry(self, width=50, textvariable=tk.StringVar(self, '3,6'))
        self.deletions.pack()

        self.insert_check = tk.Checkbutton(self, text="List of Insertions", variable=self.insert, command=DMS_ON)
        self.insert_check.pack()
        self.insert_check.deselect()
        self.insertions = tk.Entry(self, width=50, textvariable=tk.StringVar(self, 'GGG,GGGGGG'))
        self.insertions.pack()

        self.include_sub_check = tk.Checkbutton(self, text='Include Substitutions', variable=self.include_substitutions)
        self.include_sub_check.pack()
        self.include_sub_check.deselect()

        #self.matchSequences_check = tk.Checkbutton(self, text='Match Sequences', variable=self.matchSequences)
        #self.matchSequences_check.pack()

        self.run = tk.Button(self, text='Run DIMPLE', command=run).pack()

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


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry('700x900')
    app = Application(master=root)
    app.mainloop()
