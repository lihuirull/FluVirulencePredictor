# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/11/12 19:53
import os
import sys
from pathlib import Path

import pandas as pd
from joblib import load

MAPPING_DICT = {0: 'Virulent', 1: 'Avirulent'}


def explode_markers(df, input_marker_path, output_directory, prefix):
    """
    Explode the markers from the CSV file into individual mutations and create a crosstab matrix.

    Parameters:
    - df: A DataFrame for exploding.
    - input_marker_path: Path to the CSV file containing marker data.
    - output_directory: Directory where the output files will be saved.
    - prefix: Prefix for the output filenames.

    Returns:
    - DataFrame with the crosstab matrix of markers.
    """
    # Extract the filename without extension
    input_marker_filename = os.path.basename(input_marker_path).rsplit(".", 1)[0]

    # Read the CSV file into a DataFrame
    df_original = df.copy()

    # Group and sum the 'Number of Virulence markers' by 'Strain ID'
    grouped = df.groupby('Strain ID')['Number of Virulence Markers'].sum()

    # Identify strains with no adaptive markers
    mask = df['Strain ID'].map(grouped) == 0
    no_adaptive_marker_count = mask.sum()
    if no_adaptive_marker_count > 0:
        print(f"There are {no_adaptive_marker_count} strains without any markers.")

    # Explode the 'Adaptive markers' column from a comma-separated string into a list
    df['Virulence Markers'] = df['Virulence Markers'].str.split(',')

    # Expand the list into a new DataFrame and merge it back
    df = df.explode('Virulence Markers')

    # Normalize the protein type for HA
    df['Protein Type'] = df['Protein Type'].apply(lambda x: 'HA' if "H" in x else x)

    # Add a new column with mutations and protein types combined
    df['Marker_Protein'] = df['Virulence Markers'] + '_' + df['Protein Type']

    # Create a new DataFrame with columns as possible 'marker_Protein' values, rows as 'Strain ID'
    df_matrix = pd.crosstab(df['Strain ID'], df['Marker_Protein'])

    # Replace counts with 1 (presence) or 0 (absence)
    df_matrix = df_matrix.applymap(lambda x: 1 if x > 0 else 0)

    # Reset index to turn 'Strain ID' into a column
    df_matrix.reset_index(inplace = True)

    # Add labels to the results DataFrame, keeping only strains with adaptive mutations
    df_label = pd.merge(df_matrix, df_original, on = 'Strain ID')

    # Drop duplicates based on 'Strain ID'
    df_label.drop(labels = ['Virulence Markers', 'Number of Virulence Markers', 'Protein Type'], axis = 1,
                  inplace = True)

    df_label.drop_duplicates(subset = "Strain ID", keep = "first", inplace = True)

    output_matrix_filename = f'{output_directory}/{prefix}{input_marker_filename}_matrix.csv'

    df_label.to_csv(output_matrix_filename, index = False)

    return df_label


def get_explode_marker_file(input_marker_path, output_directory = ".", prefix = ""):
    """
       Processes input marker file(s) by exploding the adaptive markers into individual columns
       indicating the presence or absence of each marker in each strain. This function handles
       both individual files and directories containing multiple marker files.

       Parameters:
       - input_marker_path (str/Path): The path to the input marker file or directory.
       - output_directory (str): The directory where the output files will be saved.
       - prefix (str): The prefix to be added to the output filenames.

       Returns:
       - DataFrame: The combined DataFrame of all processed marker files, or the single processed file.
    """
    os.makedirs(output_directory, exist_ok = True)

    if Path(input_marker_path).is_dir():
        all_df = pd.DataFrame()
        for root, _, files in os.walk(input_marker_path):
            for filename in files:
                file_path = os.path.join(root, filename)
                df = pd.read_csv(file_path)
                all_df = pd.concat([df, all_df])
        result_df = explode_markers(all_df, "combined_markers.csv", output_directory, prefix)

    elif Path(input_marker_path).is_file():
        df = pd.read_csv(input_marker_path)
        print(input_marker_path)
        result_df = explode_markers(df, input_marker_path, output_directory, prefix)
    else:
        print(f"Error: {input_marker_path} is not a valid file or directory", file = sys.stderr)
        return None
    return result_df


def predict_new_data(input_marker_path, model_path, threshold, data_path, output_directory = ".",
                     prefix = ""):
    """
    Predict class labels for new data using a trained model, threshold, and a set of top features.

    Parameters:
    - input_marker_path (str): Path to the input marker file or directory containing marker files.
    - model_path (str): Path to the saved model file.
    - threshold(str): Probability threshold for model prediction
    - output_directory (str): Directory where the output files will be saved.
    - prefix (str): Prefix for the output filenames.

    Returns:
    - DataFrame: A DataFrame with strain IDs and their predicted class labels.
    """
    # Process the input marker file first
    add_prefix = prefix + "_" if prefix else ""
    processed_data = get_explode_marker_file(input_marker_path, output_directory, add_prefix)
    processed_data.set_index("Strain ID", inplace = True)
    # Load the trained model, optimal threshold, and top features
    loaded_model = load(model_path)

    # Read the top features used by the model
    top_features = pd.read_csv(f'{data_path}/max_AUC.csv').columns.tolist()

    # Process the input data
    # Reindex the processed data to include all top features from the model
    # Fill any additional features not found in the data with zeros
    processed_data_top_features = processed_data.reindex(columns = top_features, fill_value = 0)
    print(processed_data_top_features)
    # Make sure that the processed data has the same number of features as the model expects
    if processed_data_top_features.shape[1] != len(top_features):
        raise ValueError(
            f"Processed data has {processed_data_top_features.shape[1]} features,"
            f" but the model expects {len(top_features)} features.")

    # Predict probabilities for the new data
    new_data_proba = loaded_model.predict_proba(processed_data_top_features)[:, 1]

    # Predict labels based on the optimal threshold
    predictions = (new_data_proba >= threshold).astype(int)
    predictions = [MAPPING_DICT[i] for i in predictions]

    # Create a DataFrame with the predictions and strain IDs
    prediction_results = pd.DataFrame({
        'Strain ID': processed_data_top_features.index,
        'Prediction': predictions,
        'Probability': new_data_proba
    }).reset_index(drop = True)

    prediction_results.to_csv(f"{output_directory}/{add_prefix}prediction.csv")

    return prediction_results

# if __name__ == '__main__':
# predictions = predict_new_data("marker_results", '../model/gnb_model.joblib', '../model/optimal_threshold.joblib',
#                  '../model/top_features.joblib', output_directory = "matrix_results", prefix = "")
# predictions = predict_new_data("test_protein_markers.csv", '../model/gnb_model.joblib', '../model/optimal_threshold.joblib',
#                                'model/top_features.joblib', output_directory = "matrix_results", prefix = "")
# print(predictions)
