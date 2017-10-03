# Motion-Estimation-Compressed-Sensing-MRI

This repository contains MATLAB code for B-spline-based compressed sensing (SPLICS) method presented in the paper **J F P J Abascal, P Montesinos, E Marinetto, J Pascau, M Desco. Comparison of total variation with a motion estimation based compressed sensing approach for self-gated cardiac cine MRI in small animal studies. PLOS ONE 9(10): e110594, 2014.** DOI: http://dx.doi.org/10.1371/journal.pone.0110594

SPLICS generalizes spatiotemporal total variation (ST-TV) by modelling the motion between consecutive frames into the reconstruction framework. Motion is estimated using a nonrigid registration method based on hierarchical B-splines. SPLICS solves the following problem 

![](https://github.com/HGGM-LIM/Motion-Estimation-Compressed-Sensing-MRI/blob/master/SPLICSFormula.jpg)

where the first term corresponds to TV, T is the temporal sparsity operator, F is the Fourier transform, u is the reconstrcuted image and f is undersampled data.
The optimization problem is efficiently solved using the Split Bregman formulation.

## Code
This demo compares TV, spatiotemporal TV and SPLICS on cardiac cine MRI data. To run SPLICS method you need the following MATLAB toolbox: 
- FFD-based registration package. Dirk-Jan Kroon; B-spline grid, image and point based registration; 2012, retrieved from http://www.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration

## Repository files

The repository contains the following files:
Data and MATLAB functions

- **Demo_MotionEstimationCompressedSensingMRI_Sim.m:** Demo that loads data and shows how to use TV, ST-TV and SPLICS methods
- **SpatialTVSB.m:** TV method. Solved efficiently in the Fourier domain 
- **SpatioTemporalTVSB.m:** Spatiotemporal TV method. Solved efficiently in the Fourier domain 
- **SPLICS.m:** SPLICS method. Solved in the Image domain 
- **ComputeSplineRegTwoDirectionsPoolDouble.m:** Function to compute the registration between consecutive gates
- **dataCine8fr.mat:** Absolute image of retrospective cardiac cine small-animal data (8 frames, healthy rat) (Acquired data can be found at http://biig.uc3m.es/cardiac-cine-data)
- **TParameters_x5,x7,x10.mat:** Registration results needed for SPLICS method for different acceleration factors, precomputed with ComputeSplineRegTwoDirectionsPoolDouble.m
- **genPDF.m:** Function to generate a pdf with polynomial variable density sampling, by Michael Lustig 2007, downloaded from (SparseMRI V0.2), http://web.stanford.edu/~mlustig/SparseMRI.html, see M. Lustig, D.L
- **genSampling_LIM.m:** Monte-carlo algorithm to generate a sampling pattern. Modified from the original function by Michael Lustig 2007
- **maxSidePeakRatio.m:** Computes the maximum sidelobe to peak ratio of point spread function for undersampled data. Used within genSampling_LIM.m

If you use this code, please reference the paper JFPJ Abascal et al. Comparison of total variation with a motion estimation based compressed sensing approach for self-gated cardiac cine MRI in small animal studies. PLOS ONE 9(10): e110594, 2014. DOI: http://dx.doi.org/10.1371/journal.pone.0110594.   If you need to contact the author, please do so at  juanabascal78@gmail.com, paumsdv@gmail.com, desco@hggm.es
