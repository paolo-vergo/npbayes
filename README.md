# Hierarchical-Nonparametric-Bayesian-Models-to-smooth-functional-data
Project carried out during Bayesian Statistics course of PoliMi, 2021.

This repository contains a CPP implementend version of Gibbs Sampler for the following Hierarchical Dirichlet Mixture Model:



![plot](Img/formulareadme.PNG)



It was designed in order to carry out a task of functional smoothing for grouped data, in particular to estimate the step functional component of sport performances 
data contained in the Shotput Dataset. Please refer to the Project Report for details about context and further research.

Repository contais also implementations of Gibbs Samplers designed for slightly different models. Refer to specific Readme files for their formulations and purposes.

# Repository organization
An index to move through folders:

* Multiple Throws: presented version
   * Data Generation: simulated data to test the sampler
   * Results Analysis: script to analyze results produced by the sampler
   * Sampler: CPP implementation 
   * Readme : model specification and notes
* NIG: Normal-Inverse Gamma version
   * Data Generation: simulated data to test the sampler
   * Results Analysis: script to analyze results produced by the sampler
   * Sampler: CPP implementation
   * Readme : model specification and notes
* NN: Normal-Normal version
   * Data Generation: simulated data to test the sampler
   * Results Analysis: script to analyze results produced by the sampler
   * Sampler: CPP implementation
   * Readme : model specification and notes
* Shotput: 
   * Pre processing: script used to arrange data to feed them to the sampler
   * Post processing: script to analyze results produced by the sampler
   
Please note that:

   







