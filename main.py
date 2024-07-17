

# Label and entry for protein sequence input
from tkinter import *
from tkinter import scrolledtext, filedialog
import os

from PIL import Image, ImageTk

from helper_functions_script import *


def show_logo():
    # Function to display the logo image in a new window

    global logo_image  # Ensure logo_image is not garbage collected
    logo_image = Image.open("./logo_output.png")
    logo_image = logo_image.resize((logo_image.width, logo_image.height),
                                   Image.LANCZOS)
    logo_image = ImageTk.PhotoImage(logo_image)

    logo_window = Toplevel()
    logo_window.title("Logo")

    # Set the size of the Toplevel window based on the image size
    logo_window.geometry(f"{logo_image.width()+100}x{logo_image.height()+100}")

    logo_label = Label(logo_window, image=logo_image)
    logo_label.pack(padx=10, pady=10)

    def save_logo():
        save_path = filedialog.asksaveasfilename(defaultextension=".png",
                                                 filetypes=[
                                                     ("PNG files", "*.png"),
                                                     ("All files", "*.*")])
        if save_path:
            # Convert PhotoImage back to PIL.Image to save
            pil_image = Image.open("./logo_output.png")
            pil_image.save(save_path)

    save_logo_button = Button(logo_window, text="Save Logo", command=save_logo)
    save_logo_button.pack(pady=10)

    logo_label.image = logo_image

    logo_window.mainloop()





def show_gui(protein_seq, zf_list, predictions, threshold=0.5):
    # Function to create the Tkinter GUI for displaying protein sequence and highlighting domains
    root = Tk()
    root.geometry("1200x800")
    root.title("Protein Sequence Viewer")
    root.resizable(0, 0)

    label = Label(root, text="Protein Sequence Viewer",
                  font=('Times New Roman', 20))
    label.pack(pady=10)

    text_area = scrolledtext.ScrolledText(root, width=100, height=15,
                                          wrap=WORD, font=("Courier", 14))
    text_area.pack(expand=True, fill="both", padx=20, pady=10)

    # Split the protein sequence into lines and add a line space above each line
    lines = protein_seq.split('\n')
    spaced_protein_seq = '\n\n'.join(lines)

    text_area.insert(END, spaced_protein_seq)

    offset = 7
    predictions = [f"{float(pred):.2f}" for pred in predictions]


    # Highlight all ZF domains and display their predictions
    for idx, (zf_seq, prediction) in enumerate(zip(zf_list, predictions),
                                               start=1):
        if (idx>1):
            start = spaced_protein_seq.find(zf_seq) + offset*(idx-1)
        else:
            start = spaced_protein_seq.find(zf_seq)

        while start!= -1:
            end = start + len(zf_seq)
            # Adjusted start and end for the text area, accounting for the offset
            # print(
            #     f"Found ZF domain {zf_seq} at position {start}-{end} (adjusted: {adjusted_start}-{adjusted_end})")

            if float(prediction) >= threshold:
                bg_color = "green"
                fg_color = "white"
            else:
                bg_color = "blue"
                fg_color = "white"

            text_area.tag_add(f"zf_{idx}", f"1.{start}",
                              f"1.{end}")
            text_area.tag_config(f"zf_{idx}", background=bg_color,
                                 foreground=fg_color)

            # Display prediction next to the zinc finger domain
            prediction_text = f"({prediction}) "
            prediction_pos = f"1.{end}"
            text_area.insert(prediction_pos, f"({prediction}) ",
                             f"prediction_{idx}")
            text_area.tag_add(f"prediction_{idx}", prediction_pos,
                              f"{prediction_pos}+{len(str(prediction)) + 3}c")
            text_area.tag_config(f"prediction_{idx}", foreground=bg_color)

            start = spaced_protein_seq.find(zf_seq, start + 1)
            # Update the offset for the next sequence to account for the inserted prediction text
            # offset += 6

    # Add legend for color coding
    legend_frame = Frame(root)
    legend_frame.pack(pady=10)

    green_label = Label(legend_frame, text="Prediction >= Threshold",
                        bg="green", fg="white", padx=10, pady=5)
    green_label.pack(side=LEFT, padx=10)

    blue_label = Label(legend_frame, text="Prediction < Threshold", bg="blue",
                       fg="white", padx=10, pady=5)
    blue_label.pack(side=LEFT, padx=10)

    # Create buttons for each zinc finger domain that passed the threshold
    green_count = 0

    for idx, (zf_seq, prediction) in enumerate(zip(zf_list, predictions),
                                               start=1):
        if float(prediction) >= threshold:
            green_count += 1
            file_path= f"pwm_per_zf/predictions_{green_count}.txt"
            pwm_button = Button(root, text=f"Display PWM {green_count}",
                                command=lambda  path=file_path: display_pwm_file(path))
            pwm_button.pack(pady=5)
            with open(file_path, 'r') as file:
                first_line = file.readline()
                # print(f"First row of {file_path}: {first_line.strip()}")

        # Place the "Show Logo" button on the right side
    show_logo_button = Button(root, text="Show Logo", command=show_logo,
                              bg="black", fg="white")
    show_logo_button.pack(pady=10, side=RIGHT, padx=10)
    root.mainloop()


def display_pwm_file(file_path):
    # all the results are in pwm_per_zf dir
    print("file_path in display:",file_path)
    # Function to display the contents of a PWM file
    pwm_root = Tk()
    pwm_root.geometry("700x500")
    pwm_root.title("PWM Viewer")

    text_area = scrolledtext.ScrolledText(pwm_root, width=100, height=30,
                                          wrap=WORD, font=("Courier", 12))
    text_area.pack(expand=True, fill="both", padx=20, pady=10)

    # Read the content of the file and display in the text area
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            pwm_content = file.read()
            text_area.insert(END, pwm_content)
    else:
        text_area.insert(END, f"Error: File '{file_path}' not found.")

        # Function to save PWM content to a file
    def save_pwm():
        save_path = filedialog.asksaveasfilename(defaultextension=".txt",
                                                 filetypes=[
                                                     ("Text files", "*.txt"),
                                                     ("All files", "*.*")])
        if save_path:
            with open(save_path, 'w') as save_file:
                save_file.write(pwm_content)

    save_button = Button(pwm_root, text="Save", command=save_pwm)
    save_button.pack(pady=10)

    pwm_root.mainloop()


def get_input_and_display():
    # Function to get input from user and display the GUI
    protein_seq = input_seq.get("1.0", "end-1c")
    threshold = float(
    input_threshold.get()) if input_threshold.get().strip() else 0.5  # Default threshold is 0.5 if not provided

    extracted_zf_list = pre_bindzf(protein_seq, 1)

    # Example list of extracted zinc finger domains and predictions
    # extracted_zf_list = [
    #     'GKFFSWRSNLTR', 'GKSFSRSSHLIG', 'GKSFSWFSHLVT', 'GKSFVHSSRLIR',
    #     'GKSFRQSTHLIL', 'GKSYSQRSHLVV', 'GKCFSRSSHLYS', 'GKSFSQSSALIV',
    #     'GKAFIRKNDLIK', 'GIIFSQNSPFIV', 'GTALVNTSNLIG'
    # ]

    predictions = run_bindzf(1)

    # predictions = [
    #     0.8, 0.7, 0.6, 0.6,
    #     0.4, 0.4, 0.9, 1.0,
    #     0.0, 0.1, 0.7
    # ]

    # Clear the input fields
    input_seq.delete("1.0", "end")
    input_threshold.delete(0, "end")

    run_pwmPredictor(1,threshold)
    pwm_per_predictions()
    logo()

    # Call function to create GUI with protein sequence, highlighted domains, and predictions
    show_gui(protein_seq, extracted_zf_list, predictions, threshold)



# Main Tkinter window setup
win = Tk()
win.geometry("700x400")
win.title("Protein Sequence Viewer Setup")
win.resizable(0, 0)

# Label and entry for protein sequence input
label_font = ('Ariel', 30)
label_seq = Label(win, text="Enter protein sequence:", font=label_font)
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
