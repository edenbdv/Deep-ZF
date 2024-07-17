# DeepZF pipeline:

The deepZF pipeline enables biologists to input a specific protein sequence and receive as output the probability matrix of zinc finger binding preferences for each nucleotide triplet. The pipeline consists of several steps: identifying zinc finger motifs using regular expressions, passing them to the first deep learning model (BindZFpredictor), filtering the model's output using a binding threshold, transferring the filtered results to the second model (PWMpredictor), concatenating the probability matrices of all zinc fingers, converting the concatenated matrix to PWM format, and finally returning it as output to the biologist.

**Usage**

To use the deepZF pipeline:

Input a protein sequence of interest.
Run the script, which will:
Identify zinc finger motifs.
Predict binding probabilities using specified models.
Filter and concatenate results.
Convert the concatenated matrix to PWM format.
Receive the PWM matrix as output, representing zinc finger binding preferences.
