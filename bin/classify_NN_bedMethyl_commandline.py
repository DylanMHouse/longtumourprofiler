#!/usr/bin/env python
import argparse
import pandas as pd
from NN_model import NN_classifier

# Set up argument parser
parser = argparse.ArgumentParser(description="Predict methylation from BED file using NN model.")
parser.add_argument('--model', type=str, required=True, help="Path to the NN model file")
parser.add_argument('--bed', type=str, required=True, help="Path to the bedMethyl file")
parser.add_argument('--output_votes', type=str, required=True, help="Path to the output votes table")
parser.add_argument('--output_txt', type=str, required=True, help="Path to the summary txt file")

# Parse arguments
args = parser.parse_args()

# Load model and make predictions
NN = NN_classifier(args.model)
predictions, class_labels, n_features = NN.predict_from_bedMethyl(args.bed)

# Write predictions to table
df = pd.DataFrame({'class': class_labels, 'score': predictions, 'num_features': n_features})
df.to_csv(args.output_votes, sep='\t', index=False)

print(df)

# Write summary to txt file
summary = [
    'Number of features: ' + str(n_features),
    'Predicted Class: ' + str(class_labels[0]),
    'Score: ' + str(predictions[0])
]

with open(args.output_txt, 'w') as f:
    f.write("\n".join(summary))

