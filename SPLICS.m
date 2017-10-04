
function [u,varargout] = SPLICS(TParameters,G,GT,scale,f,R,N,mu,lambda,gamma,alpha,beta,nBreg,varargin)
% [u,err_tot,errAll] = SPLICS(TParameters,G,GT,scale,f,R,N,mu,lambda,gamma,alpha,beta,nBreg,varargin)
%
% Inputs: 
%
% TParameters = result from the registration between gates (see
% ComputeSplineRegTwoDirectionsPoolDouble.m). It is a structure with fields
% TParameters.GridAll, TParameters.SpacingAll, TParameters.GridDAll,
% TParameters.SpacingDAll 
% G, GT     = Forward and adjoint operators
% f         = data, Nx x Ny x time
% R         = undersampling pattern, same size as f. Entries are 1 when
% data is sampled and zero otherwise
% N         = image size n_x x n_y x n_t (number of pixels in spatial and
% temporal dimensions)
% mu        = 1, weight of data constraint. Increase for faster
% convergence, decrease for noisy data (1 usually works fine)
% lambda    = 1
% gamma     = 1e-4 to 1, L2-norm stability
% alpha     = 0.3, TV sparsity parameter
% alpha     = 0.7, Temporal sparsity parameter 
% nBreg     = Bregman iterations. This is the most critical parameter,
% chose regarding the noise in the data (lower number of iterations for
% noisier data). It can be selected by checking the convergence (errAll)
% and providing a fixed number of iterations. If a noise estimate is
% available, the discrepancy principle can be used as stoppoing criterion
%
%
% Outputs: 
%
% u         = reconstructed image of size n_x x n_y x n_t 
% err_tol   = total solution error norm at each iteration, provided the
% target image is given
% errAll    = relative solution error norm for each frame at each
% iteration, size nBreg x number of gates
%
%
% Requirements: 
%
% To run SPLICS download FFD-based registration software. 
% FFD-based registration package: For the motion estimation step we used
% the FFD-based registration package available in MATLAB Central (Dirk-Jan
% Kroon; B-spline grid,  image and point based registration; 2012,
% retrieved from 
% http://www.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration), 
% Add all directories to the path and run compile_c_files.m for faster
% implementation. 
%
%
% REF: If you use this code, please reference the paper JFPJ Abascal et al.
% Comparison of total variation with a motion estimation based compressed
% sensing approach for self-gated cardiac cine MRI in small animal studies.
% PLOS ONE 9(10): e110594, 2014. DOI:
% http://dx.doi.org/10.1371/journal.pone.0110594.  
%
% This code is an extension to the temporal dimension 
% of spatial TV Goldstein'n code mrics.m downloaded from  
% (http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html), see Tom
% Goldstein and Stanley Osher. The Split Bregman Method for L1-Regularized
% Problems. SIAM J. Imaging Sci., 2(2), 323–343.  
% 
%
% Juan FPJ Abascal
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juanabascal78@gmail.com, desco@hggm.es

h       = figure;
hw      = waitbar(0);

tolKrylov   = 1e-1;  % Krylov convergence criterion, the smaller the value
                     % the higher the precission for solving the linear system     

dimIm       = N;
rows        = N(1);
cols        = N(2);
numTime     = N(3);
Nf          = size(f);

R           = reshape(R,Nf);
f           = f.*R;

% Normalize data
normFactor  = getNormalizationFactor(f,f);
f           = normFactor*f;

if nargin >= 14
    uTarget     = varargin{1};
    errAll      = zeros(nBreg,numTime);
    err_tot     = zeros(nBreg,1);
    uTarget     = uTarget*(normFactor*scale);
    for iw = 1:numTime
        uTargetNorm(iw) = norm(reshape(uTarget(:,:,iw),[],1));
    end
end % nargin

% Reserve memory for the auxillary variables
f0      = f;
u       = zeros(rows,cols,numTime);
x       = zeros(rows,cols,numTime);
y       = zeros(rows,cols,numTime);

bx      = zeros(rows,cols,numTime);
by      = zeros(rows,cols,numTime);

dx      = zeros(rows,cols,numTime);
dy      = zeros(rows,cols,numTime);

p       = zeros(rows,cols,numTime);
bp      = zeros(rows,cols,numTime);

murf    = zeros(rows,cols,numTime);

% Forward and backward spline-based transformations
% Forward (Grid,Spacing) must be the transformation of moving frame i to
% match the next frame, and backward (GridD,SpacingD) transform frame i to
% match the frame i-1
GridAll     = TParameters.GridAll;
SpacingAll  = TParameters.SpacingAll;
GridDAll    = TParameters.GridDAll;
SpacingDAll = TParameters.SpacingDAll;

% Backprojection
for ip = 1:numTime
    murf(:,:,ip)    = mu*GT(f(:,:,ip));
end

%  Do the reconstruction
for outer = 1:nBreg;
    figure(hw); waitbar(outer/nBreg);    
    rhs_p   = lambda*OperatorLt(p-bp);
    rhs_tv  = lambda*(Dxt(x-bx)+Dyt(y-by));
    rhs     = murf + rhs_p + rhs_tv + gamma*u;
    u       = reshape(krylov(rhs(:)),N);
    
    % Derivatives
    dx      = Dx(u);
    dy      = Dy(u);
    dp      = OperatorL(u);
    
    % update x, y, p
    [x,y]   = shrink2(dx+bx, dy+by,alpha/lambda);
    p       = shrink1(dp+bp, beta/lambda);
    
    % update bregman parameters
    bx      = bx+dx-x;
    by      = by+dy-y;
    bp      = bp+dp-p;
    
    % Bregman iteration for the data constraint
    for iw = 1:numTime
        fForw         = G(u(:,:,iw)).*R(:,:,iw);
        f(:,:,iw)     = f(:,:,iw) + f0(:,:,iw)-fForw;
        murf(:,:,iw)  = mu*GT(f(:,:,iw));  
    end
    
    % Solution error norm
    if nargin >= 14
        % Compute the error
        for iw = 1:numTime
            % Relative solution error per frame
            errThis = norm(reshape(uTarget(:,:,iw)-u(:,:,iw),[],1))./uTargetNorm(iw);
            errAll(outer,iw) = errThis;
        end
        % Total solution error
        err_tot(outer)  = norm(uTarget(:)-u(:));
    end
    
    if any([outer ==1, outer == 5, rem(outer,10)==0])
        % Display image and auxiliary variables for TV and prior terms. The
        % number of nonzero coefficients on the respective transformed
        % domains are given as a precentage
        figure(h);
        subplot(2,2,1);
        imagesc(abs(u(:,:,1))/(normFactor*scale)); 
        title(['SPLICS Iter ' num2str(outer)]); colorbar; axis image; axis off;
        subplot(2,2,2);
        imagesc(abs(x(:,:,1))); axis image; axis off;
        title(['x, spatial sparsity % ' num2str(100*nnz(x(:,:,1))/(prod(N(1:2))))]);       
        colormap gray;
        subplot(2,2,3);
        imagesc(abs(p(:,:,1))); axis image; axis off;
        title(['p, temporal sparsity % ' num2str(100*nnz(p(:,:,1))/(prod(N(1:2))))]);       
        colormap gray;
        drawnow;
    end % rem
            
end % outer
close(hw);

if (nargout >= 1)
    varargout{1} = err_tot;
end % 
if (nargout >= 2)
    varargout{2} = errAll;
end % 

% undo the normalization so that results are scaled properly
u = u/(normFactor*scale);

% =====================================================================
    function uL     = OperatorL(uThis)
        % Apply the transform operator Lu_ijk=u_ijk-u_i'j'k-1 that
        % associates a pixel ij in the frame k with a unique pixel (using
        % Nearest Neighbour interpolation) in the frame k, the
        % interpolation matrix is obtained from the spline-registratio
        % toolbox.
        % NN(:,:,1) transforms a pixel ij in frame 1 to another i'j' in
        % frame 2
        %
        % This is as Tu_k = u_k-T_{k-1}(u_{k-1})
        % Each frame is substracted the forward tranformtation of the
        % previous frame
        
        uTAll       = InterpolateBasedSplines(uThis,GridAll,SpacingAll);
        uL          = zeros(dimIm);
        
        % Compare to previous frame
        for it = 1:dimIm(3)
            if it == 1
                LuThis         = uTAll(:,:,end);
            else
                LuThis         = uTAll(:,:,it-1);
            end
            
            uL(:,:,it)         = uThis(:,:,it)-LuThis;
        end
    end

    function uLt    = OperatorLt(uThis)
        %         % Apply the transpose of the transform operator
        %         % Lu_ijk=u_ijk-u_i'j'k-1
        %
        % This is as Ttu_k = u_k-Du_k+1
        % Each frame is substracted the backward tranformtation of the next
        % frame
        
        uTAll       = InterpolateBasedSplines(uThis,GridDAll,SpacingDAll);
        uLt         = zeros(dimIm);
        
        % Compare to previous frame
        for it = 1:dimIm(3)
            if it == dimIm(3)
                LuThis      = uTAll(:,:,1);
            else
                LuThis      = uTAll(:,:,it+1);
            end
            
            uLt(:,:,it)     = uThis(:,:,it)-LuThis;
        end
    end

    function uTAll = InterpolateBasedSplines(uThis,GridAll,SpacingAll)
        % Transform each image corresponding to each frame, following the
        % spline-based transformation given by GridAll and SpacingAll. This
        % should transform each frame to be similar to the next one
        
        uTAll           = zeros(dimIm);
        for ih = 1:numTime
            % Tranform u:
            IT             = uThis(:,:,ih);
            uTAll(:,:,ih)  = bspline_transform(GridAll(:,:,:,ih),IT,SpacingAll(:,:,ih));
        end % ih
        
    end

    function normFactor = getNormalizationFactor(R,f)
        
        normFactor = 1/norm(f(:)/size(R==1,1));
        
    end

    function d = Dx(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,2:cols,:) = u(:,2:cols,:)-u(:,1:cols-1,:);
        d(:,1,:) = u(:,1,:)-u(:,cols,:);
    end

    function d = Dxt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,1:cols-1,:) = u(:,1:cols-1,:)-u(:,2:cols,:);
        d(:,cols,:) = u(:,cols,:)-u(:,1,:);
    end

    function d = Dy(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
        d(1,:,:) = u(1,:,:)-u(rows,:,:);
    end

    function d = Dyt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(1:rows-1,:,:) = u(1:rows-1,:,:)-u(2:rows,:,:);
        d(rows,:,:) = u(rows,:,:)-u(1,:,:);
    end

    function [xs,ys] = shrink2(x,y,lambda)
        s = sqrt(x.*conj(x)+y.*conj(y));
        ss = s-lambda;
        ss = ss.*(ss>0);
        s = s+(s<lambda);
        ss = ss./s;
        xs = ss.*x;
        ys = ss.*y;
    end

    function xs = shrink1(x,lambda)
        s = abs(x);
        xs = sign(x).*max(s-lambda,0);
    end

    % Krylov solver subroutine to solve the linear system
    % X = GMRES(A,B,RESTART,TOL,MAXIT,M)
    % bicgstab(A,b,tol,maxit)
    function dx = krylov(r)
        %dx = gmres (@jtjx, r, 30, tolKrylov, 100);
        [dx,flag] = bicgstab(@jtjx, r, tolKrylov, 100);
    end

    % Callback function for matrix-vector product (called by krylov)
    function b = jtjx(sol)
        solMat  = reshape(sol,N);
        
        % TV part
        btv     = lambda*(Dyt(Dy(solMat)) + Dxt(Dx(solMat)));
        
        % Temporal operator part
        bP      = lambda*OperatorLt(OperatorL(solMat));
        
        % Data constraint part
        bF      = zeros(N);
        for iq = 1:N(3)
            tmp            = R(:,:,iq).*(G(solMat(:,:,iq)));
            bF(:,:,iq)     = mu*GT(tmp);
        end        
        
        b       = btv(:) + bP(:) + bF(:) + gamma*sol;
    end

end

%
