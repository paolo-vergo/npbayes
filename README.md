# Hierarchical Nonparametric Bayesian Models to Smooth Functional Data
<br>
<p align="center">
  <img src="Img/final98_500_800.png">
</p>

<br> <br>

Welcome to the Hierarchical Nonparametric Bayesian Models for Smooth Functional Data project! This repository contains an implementation of a CPP Gibbs Sampler for a Hierarchical Dirichlet Mixture Model. It was developed as part of the Bayesian Statistics course at **PoliMi** in 2021, with special thanks to **Prof. Raffaele Argiento** for his guidance.

## Project Description

The goal of this project is to address the task of functional smoothing for grouped data. Specifically, we focus on estimating the step functional component of sport performance data using hierarchical nonparametric Bayesian models. By leveraging the power of Bayesian statistics and nonparametric modeling, we aim to effectively smooth the data and extract meaningful insights.

To accomplish this, we have implemented a CPP Gibbs Sampler for the Hierarchical Dirichlet Mixture Model. The sampler utilizes a Markov chain Monte Carlo (MCMC) algorithm to perform posterior inference and estimate the underlying step functional component. It provides a flexible and powerful approach for smoothing functional data. The model has the following statistical structure: <p align="center">
  <img src="Img/formulareadme.PNG">
</p>


In addition to the Hierarchical Dirichlet Mixture Model, this repository also includes implementations of Gibbs Samplers designed for slightly different models. Each model offers a unique perspective and approach to functional smoothing. For detailed formulations and purposes of each model, please refer to the specific README files within their respective folders.

## Repository organization
An index to move through folders:

* Multiple Throws: presented version
   * Data Generation: simulated data to test the sampler                       [R]
   * Results Analysis: script to analyze results produced by the sampler       [R]
   * Sampler: CPP implementation                                               [Cpp]
   * Readme : model specification and notes
* NIG: Normal-Inverse Gamma version
   * Data Generation: simulated data to test the sampler                       [R]
   * Results Analysis: script to analyze results produced by the sampler       [R]
   * Sampler: CPP implementation                                               [Cpp]
   * Readme : model specification and notes
* NN: Normal-Normal version
   * Data Generation: simulated data to test the sampler                       [R]
   * Results Analysis: script to analyze results produced by the sampler       [R]
   * Sampler: CPP implementation                                               [Cpp]
   * Readme : model specification and notes
* Shotput: 
   * Pre processing: script used to arrange data to feed them to the sampler   [R]
   * Post processing: script to analyze results produced by the sampler        [R]
   

Please note:

1. All simulated data resulting from the Data Generation scripts are already provided, allowing you to test each part of the code. Additionally, the results of the sampling procedure are included, enabling you to evaluate and explore the outputs.
2. The real Shotput dataset is not provided directly in the repository. However, you can still test the post-processing of the sampler results using the provided data structures in the Shotput folder.

## Requirements

The implementation in this repository has minimal requirements, making it easy to set up and run. Here are the key points:

- No special requirements are needed to run the code. All R packages used are commonly available ones, making them easily accessible.
- The CPP code only requires a standard compiler, enabling seamless execution.
- For one particular package (related to Wade's code), you can download it directly from GitHub as indicated in the relevant scripts where it is used.

## Interfacing

We have not provided a specific interface between the CPP samplers and the R data generation/results analysis scripts. As a general guideline, follow these steps:

1. Generate or arrange the data to be fed into the sampler and save the created data structures.
2. Run the samplers, passing the data structures created in Step 1 as inputs. The sampler will automatically save all the results in its working directory.
3. Load the sampler results in R to perform further analysis and interpretation, using the appropriate scripts provided.

## Authors

- **Mirko Giovio**: *Politecnico di Milano - MSc in Statistical Learning*
- **Riccardo Scaramuzza**: *Politecnico di Milano - MSc in Statistical Learning*
- **Daniele Venturini**: *Politecnico di Milano - MSc in Statistical Learning*
- **Paolo Vergottini**: *Politecnico di Milano - MSc in Statistical Learning*

## Acknowledgments

We extend our heartfelt thanks to **Prof. Raffaele Argiento** (*Cattolica University of Milan*) for his invaluable guidance and support throughout the development of this project. His expertise and mentorship have been instrumental in shaping this work and enabling its success.

This project demonstrates the power and potentiality of hierarchical nonparametric Bayesian models in smoothing functional data. By combining advanced statistical techniques with a flexible Gibbs Sampler implementation, we open new doors for extracting meaningful insights from grouped data. We hope that this repository inspires further research and exploration in the field of functional data analysis, and we invite you to delve into the code and unleash the full potential of these models.
