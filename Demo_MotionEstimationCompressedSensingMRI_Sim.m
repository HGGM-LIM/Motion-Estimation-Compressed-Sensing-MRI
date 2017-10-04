% Demo_MotionEstimationCompressedSensingMRI_Sim.m
%
% If you use this code, please reference the paper JFPJ Abascal et al.
% Comparison of total variation with a motion estimation based compressed
% sensing approach for self-gated cardiac cine MRI in small animal studies.
% PLOS ONE 9(10): e110594, 2014. DOI:
% http://dx.doi.org/10.1371/journal.pone.0110594.  
%
% Code downloaded from repository: 
% https://github.com/HGGM-LIM/...
% -------------------------------------------------------------------------
%
% Demo for reconstructing cardiac cine MRI data with a B-spline-based
% compressed sensing method (SPLICS). SPLICS generalizes spatiotemporal
% total variation (ST-TV) by modelling the motion between consecutive
% frames into the reconstruction framework. SPLICS is efficiently solved
% using the Split Bregman formulation. 
% 
% TV: Solves min_u |grad_x,y u|_1  st. ||Fu-f||^2 < delta 
% using Goldstein'n code mrics.m downloaded from 
% (http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html), see Tom
% Goldstein and Stanley Osher. The Split Bregman Method for L1-Regularized
% Problems. SIAM J. Imaging Sci., 2(2), 323–343.  
%
% ST-TV: Solved using a modified version of TV that minimizes
% min_u betaxy|grad_x,y u|_1 + betat|grad_t u|_1 st. ||Fu-f||^2 < delta, 
% proposed in P Montesinos et al. Magn Reson Med., 72(2): 369–380, 2013.   
%
% SPLICS: It is a generalization of ST-TV that solves
% min_u alpha|grad_x,y u|_1 + beta|T u|_1 st. ||Fu-f||^2 < delta, 
% where T is a temporal sparsity operator that encodes motion between
% consecutive temporal frames
%
% Data corresponds to a prospective cardiac cine MRI study on a healthy rat 
% (one slice and 16 frames). 
%
% Undersampling is simulated using a modified version of Lustig's variable
% density pdf, downloaded from (SparseMRI V0.2)
% http://web.stanford.edu/~mlustig/SparseMRI.html, see M. Lustig, D.L
% Donoho and J.M Pauly "Sparse MRI: The Application of Compressed Sensing
% for Rapid MR Imaging" Magnetic Resonance in Medicine, 2007 Dec;
% 58(6):1182-1195.   
%
% Juan FPJ Abascal, Paula Montesinos
% Departamento de Bioingeniería e Ingeniería Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juanabascal78@gmail.com, paumsdv@gmail.com, desco@hggm.es
% 
% Load data: Simulated absolute image from the acquired retrospective
% cardiac cine data set with 8 frames.  
% The experimental acquired data set is available from
% http://biig.uc3m.es/cardiac-cine-data/
load('dataCine8fr','image0'); 

% Simulate data
data0   = fft2(image0);         

% Display images
figure; 
count       = 0;
for it = 1:size(image0,3)
    count       = count+1;
    subplot(3,3,count); imagesc(abs((image0(:,:,it)))); 
    axis image; axis off; colorbar; title(['image fr ' num2str(it)]); colormap gray;
end

N           = size(image0);     
% ------------------------------------------------------------------------
% Simulate data and undersampling pattern
% Pseudo-random pattern, for best results design an optimal pattern (see P
% Montesinos et al, RMR, 2013) 
%
% Parameters for Lustig's variable density pdf
rand('state',1);
DN          = [N(1),1];     
accFactor   = 7; % Acceleration factor 10, 7, 5

switch accFactor
    case 5
        sparsity    = 0.2; radius = 0; P = 4;            % x5
    case 7
        sparsity    = 0.15; radius = 0; P = 6;           % x7 
    case 10
        sparsity    = 0.1; radius = 0; P = 9;           % x10 
end 
pdf         = genPDF(DN,P,sparsity,2,radius,0);    % generates the sampling PDF
for it = 1:size(image0,3)
    temp    = genSampling_LIM(pdf,10,1);           % generates sampling pattern
    indR    = temp(:,1)==1;
    R       = zeros(N(1:2));
    R(indR,:) = 1;
%   figure; spy(R), 100*nnz(R(:,1))/N(1)        
    RAll(:,:,it)        = R;
end % it

data        = data0.*RAll;
% ------------------------------------------------------------------------
% IFFT
fr          = 3;
u_ifft      = ifft2(data(:,:,fr));
% ------------------------------------------------------------------------
% STATIC Spatial Total Variation reconstruction using Split Bregman
% Code download from
% http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html
mu          = 1;
lambda      = 1;
gamma       = 1e-4;
nInner      = 1;
nBreg       = 35;

% Goldstein's spatial TV using the Split Bregman formulation
% u_tv = mrics(RAll(:,:,1),data(:,:,1), mu, lambda, gamma, nInner, nBreg);
% 
% SpatialTVSB.m: same as mrics.m but it computes the solution error
% norm
% Reconstruction of one slice only
[u_tv,err_tv] = SpatialTVSB(RAll(:,:,fr),data(:,:,fr), mu, lambda, gamma, nInner, nBreg,image0(:,:,fr));

% ------------------------------------------------------------------------
% SpatioTemporal Total Variation with larger temporal sparsity
% Dynamic reconstruction
betaxy      = 0.5;
betat       = 0.5;
mu          = 1;
lambda      = 1;
gamma       = 1;
nInner      = 1;
nBreg       = 35;
[u_ttv,err_tot_ttv,err_ttv] = SpatioTemporalTVSB(RAll,data,N,betaxy,betat,mu,lambda,gamma,nInner,nBreg,image0);
% ------------------------------------------------------------------------
% SPLICS
% 
% Motion Estimation
% Registration step
if 0
    % To compute the registration, download the FFD registration software
    % and run this part (it takes ~2min)
    
    % Compute the registration step
%     matlabpool(4); % comment if matlabpool not available
    
    % Parameters for the registration step
    Options.Penalty     = 1e-4;
    Options.Registration= 'NonRigid';
    Options.MaxRef      = 3;
    Options.Similarity  = 'sd';
    
    % Estimate the temporal operator, computed by the registration of
    % consecutive gates from previous reconstruction (we used PBR)
    [GridAll,SpacingAll,GridDAll,SpacingDAll] = ComputeSplineRegTwoDirectionsPoolDouble(u_ttv,Options);
    TParameters.GridAll         = GridAll;
    TParameters.SpacingAll      = SpacingAll;
    TParameters.GridDAll        = GridDAll;
    TParameters.SpacingDAll     = SpacingDAll;
%      matlabpool close; 
    save('TParameters_x15','TParameters');
else
    % Load registration result already computed
    switch accFactor
        case 5
            load('TParameters_x5','TParameters');
        case 7
            load('TParameters_x7','TParameters');
        case 10
            load('TParameters_x10','TParameters');            
    end
end

mu          = 1;
lambda      = 1;
alpha       = 0.5; % TV weight
beta        = 0.5; % temporal weight
gamma       = 1e-4; % L2 stability term
nBreg       = 35;

% matlabpool(4); % uncomment if matlabpool available
scale       = sqrt(prod(N(1:2)));
G           = @(x)(fft2(x)/scale);
GT          = @(x)(ifft2(x)*scale);

[u_sp,err_tot_sp,err_sp] = SPLICS(TParameters,G,GT,scale,data,RAll,N,mu,lambda,gamma,alpha,beta,nBreg,image0);
% u = PRIMOR_CT(TParameters,G,f,R,N,uref,mu,lambda,gamma,alpha,beta,nBreg)
% [u,er
% [umocs,errmocs] = PRIMOR_CT(TParameters,G,GT,scale,data,RAll,N,uref,mu,lambda,gamma,alpha,beta,nBreg,image0);
% Reconstructed image and auxiliary variables are displayed for TV and
% prior terms, for some iteration numbers. The number of nonzero
% coefficients on the respective transformed domains are given as a
% precentage
% matlabpool close;
% ------------------------------------------------------------------------
% SpatioTemporal Total Variation with larger temporal sparsity
% Comparison of results
indX    = 60:130;
indY    = 100:170;
figure; 
tmp     = (image0(:,:,fr));
subplot(3,2,1); imagesc(abs(tmp(indX,indY))); axis image; 
axis off; colormap gray; title('Full data'); ca = caxis;
tmp     = (u_ifft);
subplot(3,2,2); imagesc(abs(tmp(indX,indY))); axis image; 
axis off; colormap gray; title('IFFT2'); caxis(ca);
tmp     = (u_tv);
subplot(3,2,3); imagesc(abs(tmp(indX,indY))); axis image; 
axis off; colormap gray; title(['STV , ' num2str(100*sparsity) '% undersampling' ]);
tmp     = (u_ttv(:,:,fr)); caxis(ca);
subplot(3,2,4); imagesc(abs(tmp(indX,indY))); axis image; 
axis off; colormap gray; title(['STTV, ' num2str(100*sparsity) ' % undersampling' ]); 
caxis(ca);
tmp     = (u_sp(:,:,fr)); caxis(ca);
subplot(3,2,5); imagesc(abs(tmp(indX,indY))); axis image; 
axis off; colormap gray; title(['SPLICS, ' num2str(100*sparsity) ' % undersampling' ]); 
caxis(ca);

figure; 
count       = 0;
for it = 1:size(image0,3)
    count       = count+1;
    subplot(3,3,count); imagesc(abs((u_ttv(:,:,it)))); 
    axis image; axis off; colorbar; title('ST-TV'); colormap gray;
end

figure; 
count       = 0;
for it = 1:size(image0,3)
    count       = count+1;
    subplot(3,3,count); imagesc(abs((u_sp(:,:,it)))); 
    axis image; axis off; colorbar; title('SPLICS'); colormap gray;
end

% Solution error 
figure; 
subplot(2,1,1); 
plot(err_tv); hold on; plot(err_ttv(:,fr),'--r'); plot(err_sp(:,fr),'-.m'); 
xlabel('Iteration number'); ylabel('Relative solution error norm'); 
legend('S-TV','ST-TV','SPLICS');
title(['Relative sol err frame ' num2str(fr)]);
subplot(2,1,2); 
plot(err_tot_ttv,'--r'); hold on; plot(err_tot_sp,'-.m'); 
xlabel('Iteration number'); ylabel('Solution error norm'); 
legend('ST-TV','SPLICS');
title('Total sol err');

%