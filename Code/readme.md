# Matlab Toolbox for behavioral electrohysiology
This toolbox provides a collection of MATLAB functions designed to facilitate the analysis of neural and behavioral data. It covers estimating tuning curves, fitting generalized linear models, decoding spike train data, detecting cell assemblies, and performing time-frequency analysis. It was written in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
It is associated to a set of [tutorials](#Tutorials) which cover the aforementioned analyses step-by-step.

## Table of Contents
- [Load functions](#load-functions)
- [Mapping neural responses along one or two variables](#mapping-neural-responses-along-one-or-two-variables)
- [Generalized Linear Models (GLMs)](#generalized-linear-models-glms)
- [Bayesian decoding](#bayesian-decoding)
- [Cell Assembly Detection](#cell-assembly-detection)
- [Cross-Correlation Analysis](#cross-correlation-analysis)
- [Time-Frequency Analysis](#time-frequency-analysis)

## Load functions

The `Load` folder contains functions for loading and preprocessing data from the [example dataset](../Data) using functions like `SetLoadParams`, `LoaddataNav`, `LoaddataSpk`, and `LoaddataLfp`. These functions should be customized to your specific datasets while keeping the main logic intact. They return structures containing different types of data: `LoaddataNav` should return a structure (called `Nav` here) containing data about experimental conditions and explanatory variables, such as the appearance of a stimulus or the position of the animal; `LoaddataSpk` should returns a structure (called `Spk` here) with spike counts sampled at the sampling rate of data stored in `Nav` together with spike times and information about the recorded cells; `LoaddataLfp` should return a structure (called `Lfp` here) that contains signals that needs to be sampled at a higher sampling rate than data in `Nav`, for instance Lfp signals for which you would like to investigate the frequency content in the gamma range.

For more information and usage examples, refer to the related [documentation](/Code/Load) and [Tutorial1](../Tutorials).

## Mapping neural responses along one or two variables

The `Maps` section includes functions for estimating tuning curves along up to two variables, using variables stored in the structure `Nav` mentioned above and a set of responses, typically corresponding to spike trains stored in `Spk`. `SetMapsParams` defines options and `MapsAnalysis` computes tuning curves and their related statistics.

For more information and usage examples, refer to the related [documentation](/Code/Maps) and [Tutorial2.1](../Tutorials).

## Generalized Linear Models (GLMs)

In the `GLMs` section, you'll find functions to estimate GLMs with up to two variables and evaluate their significance. The predictors of the GLMs are variables selected from the structure `Nav` which are used to fit a set of responses, usually corresponding to spike trains stored in `Spk`. Use `SetGLMsParams` to define options and `GLMsAnalysis` to estimate tuning curves using a Poisson GLM model.

Learn more from the [function documentation](/Code/GLMs), [Tutorial2.2](../Tutorials) and [Tutorial_handson](../Tutorials).

## Bayesian decoding

The `Decoding` section provides functions for decoding up to two variables from neural spike trains using a Bayesian approach. Again, this is using variables loaded into the structure `Nav` and a set of responses, corresponding to spike trains stored in `Spk`. Use `SetDecParams` to define parameters and `DecodingAnalysis` to estimate decoded positions and errors.

For usage examples and more information, refer to the [function documentation](/Code/Decding), [Tutorial3.1](../Tutorials) and [Tutorial_handson](../Tutorials).

## Cell Assembly Detection

The `Patterns` section contains functions for detecting cell assemblies from neural spike response data using independent component analysis (ICA). It uses variables loaded into the structure `Nav` and detect coincidental activity among a set of responses, typically corresponding to spike trains stored in `Spk`. Use `SetPattParams` to define parameters and `PatternAnalysis` to identify cell assemblies.

Detailed usage examples can be found in the [function documentation](/Code/Patterns), [Tutorial3.2](../Tutorials) and [Tutorial_handson](../Tutorials).

## Cross-Correlation Analysis

In the `Corr` section, you will find functions to estimate cross-correlations between pairs of neural responses and extract noise/trial-by-trial correlations using a shuffling to estimate signal correlations related to some varaibles stored in `Nav`. Use `SetCrossCorrParams` to set parameters and `CrossCorrAnalysis` to compute cross-correlation results.

For detailed information, usage examples, and tutorials, see the [function documentation](/Code/Corr), [Tutorial3.3](../Tutorials) and [Tutorial_handson](../Tutorials).

## Time-Frequency Analysis

In the `TFMaps` folder, you'll find functions to perform time-frequency analysis of continuous data and compute power and coherence maps along one or wo explanatory variables stored in the `Nav` structure. Use `SetTFMapsParams` to define options and `TFMapsAnalysis` to compute maps from wavelet transform and wavelet coherence between a set of recording channels. 

For detailed usage examples, refer to the [function documentation](/Code/TFMaps), [Tutorial4.1](../Tutorials) and [Tutorial_handson](../Tutorials).

  
<br>

This project is licensed under the MIT License.
  
This toolbox was developed by J. Fournier and Tulio Fernandez de Almeida in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."


