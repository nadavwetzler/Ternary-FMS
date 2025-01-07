# Ternary-FMS

A Ternary Diagram for Earthquake Focal Plane Solutions (FPS)

![FMS_ternary_map](https://github.com/user-attachments/assets/b18c32bf-ddf0-49eb-97b9-071b5c646644)

## Overview

**Ternary-FMS** provides a quick and efficient way to classify earthquake focal plane solutions (FPS) into four categories:
- **Strike-slip**
- **Normal**
- **Reverse**
- **Oblique**

The focal mechanisms are plotted on a ternary diagram for visualization. Additionally, the vector classification (`mClass`) returns the mechanism styles for further analysis or plotting.

## Features

1. Classify FPS into four categories based on dP, dT, and dB axes.
2. Generate a ternary diagram for FPS visualization.
3. Seamlessly handle focal mechanism solutions in the conventional strike, dip, and rake format.

## Usage

You can use **Ternary-FMS** in two ways:

1. **Plot the Ternary Diagram**  
   Use this mode to classify and plot FPS on a ternary diagram.  
   ```python
   cat = pd.read_csv('eq_until_2M_text.csv')
   dP, dT, dB = convert2tri(cat.Strike.values, cat.Dip.values, cat.Rake.values)
   mClass = TernaryFMS(dP, dT, dB, ax1)
2. **Only obtain the classes (no figure)**  
   Use this mode to classify and plot FPS on a ternary diagram.  
   ```python
   cat = pd.read_csv('eq_until_2M_text.csv')
   dP, dT, dB = convert2tri(cat.Strike.values, cat.Dip.values, cat.Rake.values)
   mClass = TernaryFMS(dP, dT, dB)
