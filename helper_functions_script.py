import re
import csv
import shutil

import numpy as np
import os
import pandas as pd
from shutil import copyfile
import logomaker
import matplotlib.pyplot as plt

import  main_PWMpredictor
from BindZF_predictor.code import main_bindzfpredictor_predict
import glob
import subprocess
import sys
import pickle



def find_zf_binding_domains(protein_seq):
    # Define the correct regular expression pattern with character classes for amino acids
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    protein_seq = protein_seq.replace(" ", "")
    zf_pattern = re.compile(
        r'(?P<zf_binding>[{0}]{{2}}C[{0}]{{2,4}}C(?P<zf_center>[{0}]{{12}})[H[{0}]{{3,5}}H)'.format(
            amino_acids))

    # Find all matches in the protein sequence
    matches = zf_pattern.finditer(protein_seq)
    zf_seq_list = []

    # Extract and print the matched ZF binding domains
    for match in matches:
        zf_domain = match.group('zf_binding')
        start, end = match.span()

        zf_center_bd = match.group('zf_center')
        start_center , end_center = match.span('zf_center')  # Span of the {{12}} part

        # Extract 40 amino acids on each side of the ZF binding domain
        left_context = protein_seq[max(0, start_center - 40):start_center]
        right_context = protein_seq[end_center:min(end_center + 40, len(protein_seq))]

        # Combine ZF domain and its context
        zf_with_context = left_context + zf_center_bd + right_context
        # list of tuples, where each tuple is the 12
        # center aa and the 92 length aa (center with neighbors)

        zf_seq_list.append((zf_center_bd, zf_with_context))

    return  zf_seq_list


def create_zf_dataset(zf_lst,num_protein):
    print("Creating ZF dataset...")
    # print(zf_lst)
    fieldnames = ['zf', 'label', 'seq', 'group']
    file_path = os.path.join('Results', 'zf_40_dataset1.csv')

    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for zf_sequence in zf_lst:
            # Fill 'label' and 'group' with default values of 0
            writer.writerow({'zf': zf_sequence[0], 'label': 0.0,
                             'seq': zf_sequence[1], 'group': 0})

    print(f"ZF dataset created successfully. Saved to {file_path}")


def create_pwm_dataset(num_protein,threshold=0.5):
    print("threshold in create_pwm_dataset: ", threshold)

    # Read TSV file into a pandas DataFrame without header
    tsv_file_path = f'results{num_protein}.tsv'
    tsv_df = pd.read_csv(tsv_file_path, delimiter='\t', header=None)

    # Add a column name to t+he TSV file
    tsv_df.columns = ['probabilities'] + list(tsv_df.columns[1:])


    csv_file_path = os.path.join('Results', 'zf_40_dataset1.csv')


    #csv_file_path = f'/Results/zf_40_dataset{num_protein}.csv'
    csv_df = pd.read_csv(csv_file_path)

    # Concatenate the first column of TSV with the first column of CSV
    result_df = pd.concat([csv_df.iloc[:, 0], tsv_df.iloc[:, 0]], axis=1)

    # Write the result DataFrame to a new CSV file
    result_csv_file_path = 'Data/PWMpredictor/PWM_perdictor_input.csv'
    result_df.to_csv(result_csv_file_path, index=False)

    print(f"Concatenation completed. Result saved to {result_csv_file_path}.")

     # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(result_csv_file_path)

    # Define the condition to filter rows based on the 'prob' column
    # df['probabilities_rounded'] = df['probabilities'].round(2)

    condition = df['probabilities'] >= threshold

    # Use the condition to filter rows and create a new DataFrame without those rows
    filtered_df = df[condition]

    # Save the filtered DataFrame back to a new CSV file or overwrite the original file
    filtered_df.to_csv(result_csv_file_path, index=False)

    print("the filter was done")
     # Load the updated CSV file into a pandas DataFrame
    df = pd.read_csv(result_csv_file_path)
     # Split the values in the "zf" column and create new columns 'AA1' through 'AA12'
    df[['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9', 'AA10',
        'AA11', 'AA12']] = df['zf'].apply(lambda x: pd.Series(list(x)))
    column_order = ['AA1', 'AA2', 'AA3',  'AA4', 'AA5',
                    'AA6', 'AA7', 'AA8', 'AA9', 'AA10',
                     'AA11', 'AA12']
    df =df[column_order]
    # Save the modified DataFrame back to the original CSV file or create a new CSV file
    df.to_csv(result_csv_file_path ,sep=' ' ,index=False)
    print("split to aa was done")

def clear_aa_columns():
  # Load the CSV file into a DataFrame
  file_path = 'Data/PWMpredictor/c_rc_df_copy.csv'

  df = pd.read_csv(file_path, sep=' ')  # Assuming the values are separated by space


  # Define the columns to clear
  columns_to_clear = ['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9', 'AA10', 'AA11', 'AA12','res_12']

  # Set all values in the specified columns to NaN
  df[columns_to_clear] = np.nan

  # Save the modified DataFrame back to the original CSV file (override)
  df.to_csv(file_path, index=False, sep=' ')

def override_aa_columns():
      # Load the original DataFrame
      file_path = 'Data/PWMpredictor/c_rc_df_copy.csv'

      original_df = pd.read_csv(file_path, sep=' ')

      # Load the DataFrame from another CSV file
      file_path2 = 'Data/PWMpredictor/PWM_perdictor_input.csv'
      second_df = pd.read_csv(file_path2, sep=' ')

      # Define the columns to replace
      columns_to_replace = ['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7',
                            'AA8', 'AA9', 'AA10', 'AA11', 'AA12']

      # Replace NaN values in the specified columns with values from the second DataFrame
      original_df[columns_to_replace] = original_df[columns_to_replace].fillna(
          second_df[columns_to_replace])

      # Concatenate 'AA1' to 'AA12' and store the result in 'res_12', replacing NaN values with an empty string
      original_df['res_12'] = original_df[
          ['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9',
           'AA10', 'AA11', 'AA12']].apply(
          lambda row: ''.join(row.dropna().astype(str)), axis=1).replace('',
                                                                         'NA')

      # Save the modified DataFrame back to the original CSV file (override)
      original_df.to_csv(file_path, index=False, sep=' ')

def filter_csv():
    pwm_input_file = "Data/PWMpredictor/PWM_perdictor_input.csv"
    c_rc_df_file = "Data/PWMpredictor/c_rc_df_copy.csv"
    output_file = "Data/PWMpredictor/c_rc_df_filtered.csv"

    # Read the number of lines in PWM_perdictor_input.csv
    with open(pwm_input_file, 'r') as f:
        num_lines = sum(1 for line in f)
    print("num of lines:", num_lines)

    # Read the content of c_rc_df_copy.csv up to the specified number of lines
    with open(c_rc_df_file, 'r') as f:
        content = ''.join([next(f) for _ in range(num_lines)])

    # Write the content to c_rc_df_filtered.csv
    with open(output_file, 'w') as f:
        f.write(content)


def pwm_format():

  with open("PWMpredictor/code/predictions.txt", "r") as file:
      probabilities = [float(line.strip()) for line in file]

  # Create a list of lists, where each inner list contains four probabilities
  formatted_probs = [probabilities[i:i+4] for i in range(0, len(probabilities), 4)]

  with open("PWMpredictor/code/predictions.txt", "w") as output_file:
      for row in formatted_probs:
          output_file.write("\t".join(map(str, row)) + "\n")


def pre_bindzf(protein_seq, num_protein):
    zf_lst = find_zf_binding_domains(protein_seq)
    if not zf_lst:
        print("No zinc finger binding domains found.")
        return

    create_zf_dataset(zf_lst, num_protein)

    extracted_zf_list = [item[0] for item in zf_lst]

    print("extracted_zf_list:\n", extracted_zf_list)

    return extracted_zf_list


def run_bindzf(num_protein):
    input_file_path = os.path.join('Results', 'zf_40_dataset1.csv')
    output_file_path = f"/Results/results{num_protein}.tsv"
    print(output_file_path)
    print("trying to get into model1")

    # Define the paths to your files
    # output_file_path = "BindZF_predictor/code/results1.tsv"
    output_file_path = "results1.tsv"
    model_file_path = "BindZF_predictor/code/model.p"
    encoder_file_path = "BindZF_predictor/code/encoder.p"

    python_executable = sys.executable

    # Construct the command to run your script with arguments
    command = [
        python_executable,
        "BindZF_predictor/code/main_bindzfpredictor_predict.py",
        # Path to your script
        "-in", input_file_path,
        "-out", output_file_path,
        "-m", model_file_path,
        "-e", encoder_file_path,
        "-r", "1"  # Assuming r=1 based on your previous usage
    ]

    # Execute the command
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")

    predictions = []
    try:
        with open(output_file_path, 'r') as file:
            for line in file:
                prediction_float = float(line.strip())  # Convert to float
                rounded_prediction = round(prediction_float,2)
                predictions.append(rounded_prediction)
    except FileNotFoundError:
        print(f"Error: File {output_file_path} not found.")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

    print("predictions:\n",predictions)
    return predictions

def run_pwmPredictor(num_protein, threshold=0.5):
    print("threshold in run_pwmPredictor: ", threshold)
    create_pwm_dataset(num_protein,threshold)
    src_path = 'Data/PWMpredictor/c_rc_df.csv'
    dest_path = 'Data/PWMpredictor/c_rc_df_copy.csv'
    copyfile(src_path, dest_path)

    clear_aa_columns()
    override_aa_columns()
    filter_csv()

    input_file_path = 'Data/PWMpredictor/c_rc_df_filtered.csv'
    output_file_path = 'PWMpredictor/code/predictions.txt'
    model_file_path = 'PWMpredictor/code/transfer_model100.h5'

    python_executable = sys.executable
    # Construct the command to run your script with arguments
    command = [
        python_executable,
        "main_PWMpredictor.py",
        "-i", input_file_path,
        "-o", output_file_path,
        "-m", model_file_path
    ]

    # Execute the command
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")

    pwm_format()

    file_path = 'PWMpredictor/code/predictions.txt'

    # Read the first 10 lines
    with open(file_path, 'r') as file:
        for _ in range(10):
            line = file.readline()
            if not line:
                break  # Exit the loop if EOF is reached
            print(
                line.strip())  # Print each line without trailing newline characters

    # Count the total number of lines in the file
    with open(file_path, 'r') as file:
        line_count = sum(1 for _ in file)

    print(f"Total number of lines: {line_count}")

    predictions_file_path = "PWMpredictor/code/predictions.txt"
    return predictions_file_path




def run_deepzf_for_protein(protein_seq, num_protein, threshold=0.5):
    pre_bindzf(protein_seq, num_protein)
    run_bindzf(num_protein)
    run_pwmPredictor(num_protein,threshold)
    pwm_per_predictions()
    logo()


def pwm_per_predictions():
    predictions_file_path = "PWMpredictor/code/predictions.txt"
    with open(predictions_file_path, 'r') as file:
        predictions = file.readlines()

    directory = "pwm_per_zf"


    # Ensure the directory exists or create it if it doesn't
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        # If directory exists, remove all files inside it
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")

    # Write chunks of 3 lines to separate files
    chunk_size = 3
    file_count = 1
    for i in range(0, len(predictions), chunk_size):
        chunk = predictions[i:i + chunk_size]
        file_name = os.path.join(directory, f'predictions_{file_count}.txt')
        with open(file_name, 'w') as f:
            for line in chunk:
                f.write(line)
        file_count += 1

    print(f"Files saved in directory '{directory}'.")



def logo():
    # Read data from the text file
    predictions_file_path = "PWMpredictor/code/predictions.txt"

    # Load the data from the file into a NumPy array
    matrix_data = np.loadtxt(predictions_file_path, delimiter='\t')

    # Convert NumPy array to pandas DataFrame
    df = pd.DataFrame(matrix_data, columns=['A', 'C', 'G', 'T'])

    # Print the DataFrame to verify
    print("DataFrame with updated column names:")
    print(df)

    # Create a Logo object
    logo = logomaker.Logo(df)

    # Customize the appearance if needed
    logo.style_spines(visible=False)

    # Save the logo to a PNG file
    plt.savefig('logo_output.png', format='png')