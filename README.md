# Practicals-RatGroup
This repository contains code and tutorials made for the 2023 summer school on Advanced Computational Analysis for Behavioral and Neurophysiological Recordings. During that summer school, practicals were organized into two groups that worked on different datasets. The present repository corresponds to the group that worked on electrophysiological recordings performed in the hippocampus and the amygdala while the animal performed a navigation task and during subsequent periods of sleep.

The code provides a collection of MATLAB functions designed to facilitate the analysis of neural and behavioral data. It covers various analyses, including estimating tuning curves, fitting generalized linear models, decoding from spike train data, detecting cell assemblies, performing cross-correlation analyses, and conducting time-frequency analyses. The toolbox is associated with a series of tutorials that guide you through the analyses step-by-step.

## Table of Contents

- [Code](/Code)
- [Data](/Data)
- [Tutorials](/Tutorials)

## Code

The [Code](/Code) folder is organized into several sections:
- [Load](/Code/Load): Loading and preprocessing data from the example dataset.
- [Maps](/Code/Maps): Estimating tuning curves along one or two variables.
- [GLMs](/Code/GLMs): Fitting Generalized Linear Models and evaluating their significance.
- [Decoding](/Code/Decoding): Decoding variables from spike trains using a Bayesian approach.
- [Patterns](/Code/Patterns): Detecting cell assemblies using independent component analysis.
- [Corr](/Code/Corr): Estimating cross-correlations and trial-by-trial correlations.
- [TFMaps](/Code/TFMaps): Performing time-frequency analysis on continuous data.

Each section contains specific functions to perform the respective analysis. Refer to the individual folders and function documentation for detailed information on each analysis.

## Data

The [Data](/Data) folder contains the dataset used in that group during the summer school. This dataset consists of hippocampal and bla electrophysiological recording obtained from a rat running on a linear track. The dataset is organized into different experimental conditions, and the associated neural and behavioral data are used to demonstrate the functionality of the toolbox. For more details on the dataset and its structure of the experiment, refer to the [Data README](/Data/Readme.md).

## Tutorials

The [Tutorials](/Tutorials) folder contains a series of step-by-step tutorials that demonstrate how to use the functions provided in the [Code](/Code) folder. The tutorials cover a range of topics, including loading data, estimating tuning curves, fitting GLMs, decoding variables, detecting cell assemblies, computing correlations, and performing time-frequency analysis.
The tutorials were designed to help you gain a practical understanding of the toolbox's capabilities and how to apply them to your own data.

<br>

This project is licensed under the MIT License.

The code was developed by J. Fournier and Tulio Fernandez de Almeida. Tutorials were designed by Julien Fournier, Tulio Fernandez de Almeida, Olivier Peron, Mehdi Fallahnezhad, Gabrielle Girardeau and Nicolas Gervasi. 

---

