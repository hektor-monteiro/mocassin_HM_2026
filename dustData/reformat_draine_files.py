#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 19 19:36:51 2025

@author: hmonteiro
"""

import re
import numpy as np

def reformat_dust_data(input_file, output_file):
    """
    Reformats dust optical data from an input file to a specified output file.

    Args:
        input_file (str): Path to the input data file.
        output_file (str): Path to the output data file.
    """
    try:
        with open(input_file, 'r') as infile:
            lines = infile.readlines()
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return

    data = []
    radii = []
    nrad = None
    nwav = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if "NRAD" in line:
            match = re.search(r"(\d+)\s+[\d\.E+-]+\s+[\d\.E+-]+", line)
            if match:
                nrad = int(match.group(1))
            i += 1
            continue

        if "NWAV" in line:
            match = re.search(r"(\d+)\s+[\d\.E+-]+\s+[\d\.E+-]+", line)
            if match:
                nwav = int(match.group(1))
            i += 1
            continue
        
        radius_match = re.match(r"([\d\.E+-]+)\s*=\s*radius\(micron\)", line)
        if radius_match:
            radius = float(radius_match.group(1))
            radii.append(radius)
            i += 1
            header_line = lines[i].strip()
            if "w(micron)" in header_line:
                i += 1
                for _ in range(nwav):
                    # if i < len(lines):
                    #     data_line = lines[i].strip().split()
                    #     if len(data_line) == 5:
                    #         wavelength, q_ext, q_abs, q_sca, g_cos = map(float, data_line)
                    #         data.append([radius, wavelength, q_abs, q_sca, g_cos])
                    #     i += 1
                    # else:
                    #     break
                    data_line = lines[i].split()
                    wavelength, q_abs, q_sca, g_cos = map(float, data_line)
                    data.append([radius, wavelength, q_abs, q_sca, g_cos])
                    i += 1
            else:
                print(f"Warning: Expected header line not found after radius at line {i}.")
        else:
            i += 1

    if data:
        header = "radius(micron) w(micron) Q_abs Q_sca g=<cos>"
        np.savetxt(output_file, data, header=header, comments='', fmt='%.6E')
        print(f"Successfully reformatted data and saved to '{output_file}'.")
    else:
        print("No data found to reformat.")

if __name__ == "__main__":
    input_file = "/home/hmonteiro/Downloads/mocassin_HM_2025/mocassin_HM_2025/dustData/Sil_21"  # Replace with the name of your input file
    output_file = "/home/hmonteiro/Downloads/mocassin_HM_2025/mocassin_HM_2025/dustData/Sil_21.cat" # Replace with the desired name for your output file
    reformat_dust_data(input_file, output_file)