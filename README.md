# DeepZF pipeline

The deepZF pipeline enables biologists to input a specific protein sequence and receive as output the probability matrix of zinc finger binding preferences for each nucleotide triplet. The pipeline consists of several steps:

1. **Identifying Zinc Finger Motifs**: Using regular expressions to detect zinc finger motifs within the protein sequence.
![WhatsApp Image 2024-07-17 at 15 55 06](https://github.com/user-attachments/assets/db19df24-59ea-4cdd-acc9-573c4324916e)

2. **First Deep Learning Model (BindZFpredictor)**:
   - Predicting binding probabilities for identified zinc finger motifs.
   - Filtering the model's output using a specified binding threshold.
   
3. **Second Model (PWMpredictor)**:
   - Transferring filtered results to predict PWM (Position Weight Matrix) preferences.
   
4. **Concatenating Probability Matrices**:
   - Combining the probability matrices of all zinc fingers into a unified matrix.
   
5. **Converting to PWM Format**:
   - Converting the concatenated matrix into PWM format, which represents zinc finger binding preferences.
   
6. **Output**:
   - Providing the PWM matrix as the final output, reflecting zinc finger binding preferences derived from the input protein sequence.
   - In addition,the user will be able to see the LogoSeq for all fingers whose probability is above the given threshold.

![‏‏4](https://github.com/user-attachments/assets/725090db-aaa9-4e76-ad58-2b346caf29f8)
   - ![‏‏3](https://github.com/user-attachments/assets/de1a0965-c2fa-4931-899b-58e1a5bf3a0f)


**Usage**

To use the deepZF pipeline:
- Input a protein sequence of interest.
- Run the script, which will:
  - Identify zinc finger motifs.
  - Predict binding probabilities using specified models.
  - Filter and concatenate results.
  - Convert the concatenated matrix to PWM format.
- Receive the PWM matrix as output, representing zinc finger binding preferences.
