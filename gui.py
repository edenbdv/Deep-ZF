import tkinter as tk
from tkinter import scrolledtext

def create_gui(protein_seq, zf_list):
    root = tk.Tk()
    root.title("Protein Sequence Viewer")

    text_area = scrolledtext.ScrolledText(root, width=80, height=10, wrap=tk.WORD, font=("Courier", 12))
    text_area.pack(expand=True, fill="both")

    text_area.insert(tk.END, protein_seq)

    # Highlight ZF domains in different colors
    for zf_seq in zf_list:
        start = protein_seq.find(zf_seq)
        end = start + len(zf_seq)

        while start != -1:
            text_area.tag_add("zf", f"1.{start}", f"1.{end}")
            text_area.tag_config("zf", background="yellow", foreground="black")
            start = protein_seq.find(zf_seq, start + 1)
            end = start + len(zf_seq)

    root.mainloop()

# Example usage:
protein_seq = "MDAKSLTAWSRTLVTFKDVFVDFTREEWKLLDTAQQIVYRNVMLENYKNLVSLGYQLTKPDVILRLEKGEEPWLVEREIHQETHPDSETAFEIKSSVSSRSIFKDKQSCDIKMEGMARNDLWYLSLEEVWKCRDQLDKYQENPERHLRQVAFTQKKVLTQERVSESGKYGGNCLLPAQLVLREYFHKRDSHTKSLKHDLVLNGHQDSCASNSNECGQTFCQNIHLIQFARTHTGDKSYKCPDNDNSLTHGSSLGISKGIHREKPYECKECGKFFSWRSNLTRHQLIHTGEKPYECKECGKSFSRSSHLIGHQKTHTGEEPYECKECGKSFSWFSHLVTHQRTHTGDKLYTCNQCGKSFVHSSRLIRHQRTHTGEKPYECPECGKSFRQSTHLILHQRTHVRVRPYECNECGKSYSQRSHLVVHHRIHTGLKPFECKDCGKCFSRSSHLYSHQRTHTGEKPYECHDCGKSFSQSSALIVHQRIHTGEKPYECCQCGKAFIRKNDLIKHQRIHVGEETYKCNQCGIIFSQNSPFIVHQIAHTGEQFLTCNQCGTALVNTSNLIGYQTNHIRENAY"
extracted_zf_list = ['GKFFSWRSNLTR', 'GKSFSRSSHLIG', 'GKSFSWFSHLVT', 'GKSFVHSSRLIR', 'GKSFRQSTHLIL', 'GKSYSQRSHLVV', 'GKCFSRSSHLYS', 'GKSFSQSSALIV', 'GKAFIRKNDLIK', 'GIIFSQNSPFIV', 'GTALVNTSNLIG']

create_gui(protein_seq, extracted_zf_list)
