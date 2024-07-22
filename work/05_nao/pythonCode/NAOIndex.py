#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: pamelatueam
"""
import matplotlib.pyplot as plt
import pandas as pd

# Specify the path to your .txt file
file_path = 'work/05_nao/data/nao_index_data.txt'

# Read the file 
df = pd.read_csv(file_path, delimiter=r'\s+')  # Change delimiter if necessary (e.g., ',' for CSV)

# Print the DataFrame
# print(df)
colors = ['blue' if value < 0 else 'red' for value in df['Index']]

plt.bar(df['Year'], df['Index'], color=colors)

# Add labels and title
plt.xlabel('Year')
plt.ylabel('NAO Index')
plt.title('NAO Index DJFM 1864-2023' , fontsize=20)
plt.grid(True)

# Show the plot
plt.show()
