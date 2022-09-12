# Near-Instantaneous-Battery-End-of-Discharge-Prognosis-via-Uncertain-Event-Likelihood-Functions
This tutorial is intended to indicate how to run supplementary files that replicate 
the content of the article entitled "Near-Instantaneous Battery End-of-Discharge 
Prognosis via Uncertain Event Likelihood Functions", published in the journal ISA 
Transactions, Elsevier.

In this directory, there should be three elements whose description is the
following:

- 'README.md':  Instructions for executing supplementary files.
- 'Markov Chain': Folder containing executable files for the case study
of Section 5.2.1., where exogenous inputs are characterized by a Markov Chain.
- 'Gaussian mixture': Folder containing executable files for the case study
of Section 5.2.2., where exogenous inputs are characterized by a Gaussian mixture.

Each folder in turn contains the following files:

a) 'near_instantaneous.m': MATLAB executable script.
b) 'ground_truth_eta01.mat': Loadable file in MATLAB.
c) 'ground_truth_eta09.mat': Loadable file in MATLAB.

Open the folder associated with the case study you want to verify.

In order to corroborate the results presented in the aforementioned scientific 
article, make sure that the files described in b), c) and d) are in the same 
directory.

The file 'near_instantaneous.m' should be opened in MATLAB. The file is subdivided 
into four sections:

1) 'Authorship': Author and contact information.
2) 'Definitions': Here the functions and parameters described in the article 
are defined. Reference is made to the corresponding equation or table number.
3) 'Method 1 (Standard): Certain event approach': Implementation of Method 1 
described in Section IV of the article.
4) 'Method 2 (Proposed): Uncertain event approach': Implementation of Method 2 
described in Section IV of the article.

To run the script 'near_instantaneous.m' file, follow the step-by-step instructions 
below:

- Go to section 2) and assign to the variable 'sigma_eta' a value of 0.1 or 0.9, 
depending on the case to be corroborated.
- Only execute section 2).
- Go to either section 3) or section 4), depending on the method you want to run. 
Once there, assign to the variable 'N' (number of Monte Carlo simulations) the 
integer value (greater than zero) that you want to test.
- Execute only section 3) or only section 4), depending on the method you want to 
check. Wait a few moments and a graph similar to the ones in Figs. 4 and 6 of the 
article will be displayed. Performance metrics for the method being tested will 
also be displayed on the command line.
