#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 21:59:17 2025

@author: hmonteiro
"""

import numpy as np

def convert_draine_OptData_file(input_filename="callindex.out_sil.D03", output_filename="astrosil.nk"):
    """
    Converts a Draine's astronomical silicate optical constant file
    to a new file with wavelength, n, and k columns.

    Args:
        input_filename (str): The name of the input data file (e.g., "callindex.out_sil.D03").
        output_filename (str): The name of the output file (e.g., "astrosil_nk.txt").
    """
    try:
        # Read the header to find the number of wavelengths and data start line
        with open(input_filename, 'r') as f_in:
            lines = f_in.readlines()
            data_start_line = 0
            for i, line in enumerate(lines):
                if "wave(um)" in line:
                    data_start_line = i + 1
                    break
            if data_start_line == 0:
                raise ValueError("Could not find 'wave(um)' header line in the input file.")

        # Load the data using numpy.genfromtxt, skipping the header lines
        data = np.genfromtxt(input_filename, skip_header=data_start_line)

        # Extract the relevant columns
        # Column 0: wavelength (um)
        # Column 3: Re(n)-1
        # Column 4: Im(n) which is k
        wavelength_raw = data[:, 0]
        re_n_minus_1 = data[:, 3]
        im_n = data[:, 4]

        # Calculate n by adding 1 to Re(n)-1
        n_values = re_n_minus_1 + 1.0
        k_values = im_n # Im(n) is already k

        # Combine into a new array for sorting
        combined_data = np.column_stack((wavelength_raw, n_values, k_values))

        # Sort the combined_data array by the first column (wavelength) in ascending order
        # np.argsort returns the indices that would sort an array
        # We then use these indices to reorder the rows of combined_data
        sorted_indices = np.argsort(combined_data[:, 0])
        sorted_output_data = combined_data[sorted_indices]

        # Save to the new file
        header = "wavelength        n          k"
        np.savetxt(output_filename, sorted_output_data, fmt='%.6e', header=header, comments='# ')

        print(f"Successfully converted and sorted '{input_filename}' to '{output_filename}'.")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    
    # Example usage:
    convert_draine_OptData_file(input_filename="callindex.out_silD03", output_filename="astroSilD2003.nk")

