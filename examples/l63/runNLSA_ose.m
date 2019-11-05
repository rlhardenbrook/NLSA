function [ model, In, Out ] = runNLSA( experiment, iProc, nProc, ifPlot );
%
% This function creates an NLSA model and  executes the various NLSA steps 
% for data generated by the dynamical system on the 2-torus in torusData.m  
%
% Each step saves the results on disk, so one can resume a partially completed
% calculation by commenting out the steps in this script which have been 
% already exectuted. 
%
% Similarly, if one changes the NLSA model parameters specified in the function
% torusNLSAModel_ose, it is only necessary to repeat the steps that are
% affected by the parameter changes. (E.g., if the diffusion maps bandwidth 
% parameter is changed, the steps up to the distance symmetrization can be
% skipped.)
%
% Input arguments:
%
% experiment:   a string identifier for the NLSA model, passed to
%               the function torusNLSAModel_ose
%
% iProc, nProc: These arguments provide rudimentary paralelization features for 
%               certain steps in the code that support it. Setting nProc > 1 
%               means that the computation is divided into nProc batches. These
%                batches can be executed in parallel by launching nProc 
%               instances of Matlab and running this function with iProc set to
%               1 for instance #1, 2 for instance #2, ...     
%               
% ifPlot:       Set to true to make basic eigenfunction scatterplots
%
% To display the optimal bandwidth from automatic bandwidth for the diffusion 
% operator, run the following command (only available for gl_mb diffusion 
% operators):
%
% [ epsilonOpt, Info ] = computeOptimalBandwidth( model );
%
% To recover the computed eigenfunctions and the projected and reconstructed 
% data, run the following commands:
%
% Diffusion eigenfunctions:
% phi = getDiffusionEigenfunctions( model ); 
% 
% Projected target data onto the diffusion eigenfunctions:
% a = getProjectedData( model );
%
% Reconstructed data: 
% x = getReconstructedData( model );
%
% Diffusion eigenfunctions (out-of-sample-extension): 
% phiO = getOseDiffusionEigenfunctions( model );
%
% Reconstructed data (out-of-sample extension):
% xO = getOseReconstructedData( model );
%  
% Modified 2016/02/02

%% Default input arguments
if nargin == 0
    experiment = 'test'; 
end
if nargin <= 1 
    iProc = 1;
    nProc = 1;
end

if nargin <= 3
    ifPlot = false;
end

%% NLSA Model
disp( experiment )
[ model, In, Out ] = l63NLSAModel_ose( experiment ); 


%% NLSA steps
disp( 'Takens delay embedding' ); computeDelayEmbedding( model )

% The next step is only needed for velocity-dependent distances such as 
% the "at" and "cone" distances
disp( 'Phase space velocity' ); computeVelocity( model )

% The next step is only needed if the target data are different from the 
% source data
disp( 'Takens delay embedding, target data' ); computeTrgDelayEmbedding( model )

fprintf( 'Pairwise distances, %i/%i\n', iProc, nProc ); 
computePairwiseDistances( model, iProc, nProc )

disp( 'Distance symmetrization' ); symmetrizeDistances( model )

% The next step is only needed for automatic bandwidth selection
disp( 'Kernel sum' ); computeKernelDoubleSum( model )

disp( 'Diffusion eigenfunctions' ); computeDiffusionEigenfunctions( model )

disp( 'Projection of target data onto diffusion eigenfunctions' );
%computeProjection( model );

disp( 'Reconstruction of the projected data' )
%computeReconstruction( model );

disp( 'Takens delay embedding, out-of-sample data' )
computeOutDelayEmbedding( model )

% The next step is only needed for velocity-dependent distances
disp( 'Phase space velocity, out-of-sample data' )
computeOutVelocity( model )

fprintf( 'OSE pairwise distances, %i/%i\n', iProc, nProc ); 
computeOsePairwiseDistances( model, iProc, nProc )

disp( 'OSE kernel normalization' ); computeOseKernelNormalization( model );

disp( 'OSE kernel degree' ); computeOseKernelDegree( model );

disp( 'OSE operator' ); computeOseDiffusionOperator( model );

disp( 'OSE eigenfunctions' ); computeOseDiffusionEigenfunctions( model );

%disp( 'Reconstruction of out-of-sample data in Takens delay-coordinate space' )
%computeOseEmbData( model )

%disp( 'Reconstruction of out-of-sample data in physical space' )
%computeOseReconstruction( model )

%% BASIC PLOTS
if ifPlot
    x  =  getSrcData( model );                     % in-sample data
    xO =  getOutData( model );                     % out-of-sample data
    phi = getDiffusionEigenfunctions( model );     % in-sample eigenfunctions
    phiO = getOseDiffusionEigenfunctions( model ); % out-of-sample extensions

    idxPhiPlot = 5; % eigenfunction to plot

    % The following only works for data embedded in R3

    % In-sample eigenfunctions
    figure;
    scatter3( x( 1, In.nXB + 1 : end - In.nXA ), ...
              x( 2, In.nXB + 1 : end - In.nXA ), ...
              x( 3, In.nXB + 1 : end - In.nXA ), 10, phi( :, idxPhiPlot ) )
    view( 20, 80 )
    title( [ 'phi' int2str( idxPhiPlot ) ', in sample' ] )

    % Out-of-sample extended eigenfunctions
    figure; 
    scatter3( xO( 1, Out.nXB + 1 : end - Out.nXA ), ...
              xO( 2, Out.nXB + 1 : end - Out.nXA ), ...
              xO( 3, Out.nXB + 1 : end - Out.nXA ), 10, phiO( :, idxPhiPlot ) )
    view( 20, 80 )
    title( [ 'phi' int2str( idxPhiPlot ) ', out-of-sample' ] )
end