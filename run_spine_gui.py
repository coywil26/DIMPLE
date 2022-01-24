# RUN SPINE
# script for usage with command line

import argparse
from SPINE.SPINE import align_genevariation, generate_DIS_fragments, print_all, post_qc, addgene, SPINEgene, generate_DMS_fragments
from Bio.Seq import Seq
import tkinter as tk
from tkinter import filedialog
from tkinter.tix import Balloon

def run():
    if app.wDir is None:
        app.wDir = app.geneFile.rsplit('/', 1)[0]+'/'

    if any([x not in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] for x in app.handle.get()]):
        raise ValueError('Genetic handle contains non nucleic bases')

    SPINEgene.handle = app.handle.get()
    SPINEgene.synth_len = int(app.oligoLen.get())
    overlapL = int(app.overlap.get())
    overlapR = int(app.overlap.get())
    if app.delete.get():
        overlapR = max([int(x) for x in app.deletions.get().split(',')])+overlapR
    if app.fragmentLen.get() != 'auto':
        SPINEgene.maxfrag = int(app.fragmentLen.get())
    else:
        SPINEgene.maxfrag = int(app.oligoLen.get()) - 62 - overlapL - overlapR  # 62 allows for cutsites and barcodes

    #adjust primer primerBuffer
    SPINEgene.primerBuffer += overlapL

    SPINEgene.barcodeF = SPINEgene.barcodeF[int(app.barcode_start.get()):]
    SPINEgene.barcodeR = SPINEgene.barcodeR[int(app.barcode_start.get()):]
    SPINEgene.cutsite = Seq(app.restriction_sequence.get())
    SPINEgene.avoid_sequence = [Seq(x) for x in app.avoid_sequence.get().split(',')]
    if app.usage:
        SPINEgene.usage = 'ecoli'
    else:
        SPINEgene.usage = 'human'

    OLS = addgene(app.geneFile)

    if app.matchSequences.get() == 'match':
        align_genevariation(OLS)
    if app.mutationType.get() == 0:
        generate_DIS_fragments(OLS, overlapL, app.wDir)
    elif app.mutationType.get() == 1:
        if app.delete.get() == 0:
            deletions = 0
        else:
            deletions = [int(x) for x in app.deletions.get().split(',')]
        if app.insert.get() == 0:
            insertions = 0
        else:
            insertions = app.insertions.get().split(',')
        generate_DMS_fragments(OLS, overlapL, overlapR, app.include_substitutions, insertions, deletions, app.wDir)
    else:
        raise AttributeError('Did not select type of mutation')
    post_qc(OLS)
    print_all(OLS, app.wDir)

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.winfo_toplevel().title("SPINE Saturated Programmable INsertion Engineering")
        self.matchSequences = tk.IntVar()
        self.mutationType = tk.IntVar()
        self.usage = tk.IntVar()
        self.include_substitutions = tk.IntVar()
        self.delete = tk.IntVar()
        self.insert = tk.IntVar()

        self.wDir_file = tk.Button(self, text='Working Directory', command=self.browse_wDir)
        self.wDir_file.pack()
        self.wDir = None

        self.gene_file = tk.Button(self, text='Target Gene File', command=self.browse_gene)
        self.gene_file.pack()
        self.geneFile = None
        #wDir_balloon = Balloon(root, bg="white", title="Help")
        #wDir_balloon.bind_widget(self.gene, balloonmsg="sequences including backbone in a fasta format. Place all in one fasta file. Name description can include start and end points (>gene1 start:1 end:2)")

        tk.Label(self, text='Oligo Length').pack()
        self.oligoLen = tk.Entry(self, textvariable=tk.StringVar(self, '230'))
        self.oligoLen.pack()

        tk.Label(self, text='Fragment Length').pack()
        self.fragmentLen = tk.Entry(self, textvariable=tk.StringVar(self, 'auto'))
        self.fragmentLen.pack()

        tk.Label(self, text='Fragment Overlap (This will change if deletions are selected)').pack()
        self.overlap = tk.Entry(self, textvariable=tk.StringVar(self, '3'))
        self.overlap.pack()

        tk.Label(self, text='Barcode Start position').pack()
        self.barcode_start = tk.Entry(self, textvariable=tk.StringVar(self, '0'))
        self.barcode_start.pack()

        tk.Label(self, text='Type II restriction sequence').pack()
        self.restriction_sequence = tk.Entry(self, textvariable=tk.StringVar(self, 'CGTCTC'))
        self.restriction_sequence.pack()

        tk.Label(self, text='Sequences to avoid').pack()
        self.avoid_sequence = tk.Entry(self, width=50, textvariable=tk.StringVar(self, 'CGTCTC, GGTCTC'))
        self.avoid_sequence.pack()

        tk.Label(self, text='Insertion Handle').pack()
        self.handle = tk.Entry(self, width=50, textvariable=tk.StringVar(self, 'AGCGGGAGACCGGGGTCTCTGAGC'))
        self.handle.pack()

        tk.Label(self, text='Type of mutations to generate', font="helvetica 12 underline").pack()
        self.DIS_check = tk.Radiobutton(self, text="Deep Insertional Scan (Genetic Handle)", variable=self.mutationType, value=0)
        self.DIS_check.pack()
        self.DMS_check = tk.Radiobutton(self, text="Deep Mutational Scan", variable=self.mutationType, value=1)
        self.DMS_check.pack()
        self.DIS_check.select()

        tk.Label(self, text='Codon Usage', font="helvetica 12 underline").pack()
        self.ecoli_check = tk.Radiobutton(self, text="E. coli", variable=self.usage, value=1)
        self.ecoli_check.pack()
        self.human_check = tk.Radiobutton(self, text="Human", variable=self.usage, value=0)
        self.human_check.pack()

        def DMS_ON():
            if self.delete.get() or self.insert.get():
                self.DMS_check.select()

        tk.Label(self, text='Settings for Indels', font="helvetica 12 underline").pack()
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

        #self.matchSequences_check = tk.Checkbutton(self, text='Match Sequences', variable=self.matchSequences)
        #self.matchSequences_check.pack()

        self.run = tk.Button(self, text='Run SPINE', command=run).pack()

        #self.output = tk.Text(self, height=10, width=60).pack()

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
    root.geometry('700x800')
    app = Application(master=root)
    app.mainloop()
