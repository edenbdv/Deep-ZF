from tkinter import *
from tkinter import scrolledtext
import os
from helper_functions_script import *


def show_gui(protein_seq, zf_list, predictions, threshold=0.5):
    # Function to create the Tkinter GUI for displaying protein sequence and highlighting domains
    root = Tk()
    root.geometry("900x600")
    root.title("Protein Sequence Viewer")
    root.resizable(0, 0)

    label = Label(root, text="Protein Sequence Viewer", font=('Times New Roman', 16))
    label.pack(pady=10)

    text_area = scrolledtext.ScrolledText(root, width=100, height=15, wrap=WORD, font=("Courier", 12))
    text_area.pack(expand=True, fill="both", padx=20, pady=10)

    # Split the protein sequence into lines and add a line space above each line
    lines = protein_seq.split('\n')
    spaced_protein_seq = '\n\n'.join(lines)

    text_area.insert(END, spaced_protein_seq)

    # Highlight all ZF domains and display their predictions
    for idx, (zf_seq, prediction) in enumerate(zip(zf_list, predictions), start=1):
        start = spaced_protein_seq.find(zf_seq)
        end = start + len(zf_seq)

        while start != -1:
            if prediction >= threshold:
                bg_color = "green"
                fg_color = "white"
            else:
                bg_color = "blue"
                fg_color = "white"

            text_area.tag_add(f"zf_{idx}", f"1.{start}", f"1.{end}")
            text_area.tag_config(f"zf_{idx}", background=bg_color, foreground=fg_color)

            # Display prediction next to the zinc finger domain
            prediction_pos = f"1.{end}"
            text_area.insert(prediction_pos, f"({prediction}) ", f"prediction_{idx}")
            text_area.tag_add(f"prediction_{idx}", prediction_pos, f"{prediction_pos}+{len(str(prediction)) + 3}c")
            text_area.tag_config(f"prediction_{idx}", foreground=bg_color)

            start = spaced_protein_seq.find(zf_seq, start + 1)
            end = start + len(zf_seq)
    # Add legend for color coding
    legend_frame = Frame(root)
    legend_frame.pack(pady=10)

    green_label = Label(legend_frame, text="Prediction >= Threshold", bg="green", fg="white", padx=10, pady=5)
    green_label.pack(side=LEFT, padx=10)

    blue_label = Label(legend_frame, text="Prediction < Threshold", bg="blue", fg="white", padx=10, pady=5)
    blue_label.pack(side=LEFT, padx=10)
    # Create buttons for each zinc finger domain that passed the threshold
    for idx, (zf_seq, prediction) in enumerate(zip(zf_list, predictions), start=1):
        if prediction >= threshold:
            pwm_button = Button(root, text=f"Display PWM {idx}",
                                command=lambda seq=idx: display_pwm_file(f"pred_protein{seq}.txt"))
            pwm_button.pack(pady=5)

    root.mainloop()


def display_pwm_file(file_path):
    # all the results are in pwm_per_zf dir


    # Function to display the contents of a PWM file
    pwm_root = Tk()
    pwm_root.geometry("700x500")
    pwm_root.title("PWM Viewer")

    text_area = scrolledtext.ScrolledText(pwm_root, width=100, height=30, wrap=WORD, font=("Courier", 12))
    text_area.pack(expand=True, fill="both", padx=20, pady=10)

    # Read the content of the file and display in the text area
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            pwm_content = file.read()
            text_area.insert(END, pwm_content)
    else:
        text_area.insert(END, f"Error: File '{file_path}' not found.")

    pwm_root.mainloop()


def get_input_and_display():
    # Function to get input from user and display the GUI
    protein_seq = input_seq.get("1.0", "end-1c")
    threshold = float(
        input_threshold.get()) if input_threshold.get().strip() else 0.5  # Default threshold is 0.5 if not provided

    extracted_zf_list = pre_bindzf(protein_seq,1)

    # Example list of extracted zinc finger domains and predictions
    # extracted_zf_list = [
    #     'GKFFSWRSNLTR', 'GKSFSRSSHLIG', 'GKSFSWFSHLVT', 'GKSFVHSSRLIR',
    #     'GKSFRQSTHLIL', 'GKSYSQRSHLVV', 'GKCFSRSSHLYS', 'GKSFSQSSALIV',
    #     'GKAFIRKNDLIK', 'GIIFSQNSPFIV', 'GTALVNTSNLIG'
    # ]


    predictions =  run_bindzf(1)

    # predictions = [
    #     0.8, 0.7, 0.6, 0.61,
    #     0.4, 0.49, 0.9, 1,
    #     0, 0.1, 0.7
    # ]

    # Clear the input fields
    input_seq.delete("1.0", "end")
    input_threshold.delete(0, "end")



    # Call function to create GUI with protein sequence, highlighted domains, and predictions
    show_gui(protein_seq, extracted_zf_list, predictions, threshold)


# Main Tkinter window setup
win = Tk()
win.geometry("700x400")
win.title("Protein Sequence Viewer Setup")
win.resizable(0, 0)

# Label and entry for protein sequence input
label_seq = Label(win, text="Enter protein sequence:")
label_seq.pack(pady=10)
input_seq = Text(win, height=10, width=70, wrap=WORD)
input_seq.pack(pady=10)

# Label for threshold input (optional)
label_threshold = Label(win, text="Enter threshold (optional):")
label_threshold.pack()

# Entry for threshold input
input_threshold = Entry(win, width=30)
input_threshold.pack(pady=10)

# Button to submit inputs and display GUI
submit_button = Button(win, text="Submit", command=get_input_and_display)
submit_button.pack()

win.mainloop()

# protein_seq = "MDAKSLTAWSRTLVTFKDVFVDFTREEWKLLDTAQQIVYRNVMLENYKNLVSLGYQLTKPDVILRLEKGEEPWLVEREIHQETHPDSETAFEIKSSVSSRSIFKDKQSCDIKMEGMARNDLWYLSLEEVWKCRDQLDKYQENPERHLRQVAFTQKKVLTQERVSESGKYGGNCLLPAQLVLREYFHKRDSHTKSLKHDLVLNGHQDSCASNSNECGQTFCQNIHLIQFARTHTGDKSYKCPDNDNSLTHGSSLGISKGIHREKPYECKECGKFFSWRSNLTRHQLIHTGEKPYECKECGKSFSRSSHLIGHQKTHTGEEPYECKECGKSFSWFSHLVTHQRTHTGDKLYTCNQCGKSFVHSSRLIRHQRTHTGEKPYECPECGKSFRQSTHLILHQRTHVRVRPYECNECGKSYSQRSHLVVHHRIHTGLKPFECKDCGKCFSRSSHLYSHQRTHTGEKPYECHDCGKSFSQSSALIVHQRIHTGEKPYECCQCGKAFIRKNDLIKHQRIHVGEETYKCNQCGIIFSQNSPFIVHQIAHTGEQFLTCNQCGTALVNTSNLIGYQTNHIRENAY"