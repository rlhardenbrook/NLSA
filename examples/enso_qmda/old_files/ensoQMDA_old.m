%% QUANTUM MECHANICAL DATA ASSIMILATION OF THE EL NINO SOUTHERN OSCILLATION
%
% Observed data is Indo-Pacific SST.
% Predicted data is Nino 3.4 index.
%
% Modified 2021/07/05
    
%% SCRIPT EXECUTION OPTIONS
% Data extraction from netCDF files
ifDataSrc    = false; % source data (training phase)
ifDataTrg    = false; % target data (training phase)
ifDataObs    = false; % observed data (training phase)
ifDataOutObs = false; % observed data (data assimilation phase) 
ifDataOutTrg = false; % target data (data assimilation phase)

% NLSA
ifNLSA    = true; % perform NLSA (source data, training phase)
ifNLSAObs = true; % perform NLSA (observed data, training phase)
ifNLSAOut = true; % perform out-of-sample extension (data assimilation phase)

% Data assimilation
ifDATrainingData  = true; % read training data for data assimilation
ifDATestData      = true; % read test data for data assimilation
ifKoopmanOp       = true; % compute Koopman operators
ifObservableOp    = true; % compute quantum mechanical observable operators  
ifAutotune        = true;  % tune the observation kernel
ifFeatureOp       = false; % compute quantum mechanical feature operators
ifFeatureOp_sfunc = true; % feature operators using kernel shape function 
ifDA              = true; % perform data assimilation
ifDAErr           = true; % compute forecast errors

% Data assimilation options
ifVB              = true; % use variable-bandwidth kernel  
ifSqrtm           = false; % apply square root to matrix-valued feature map

% IO options
ifSaveData        = true;  % save DA output to disk
ifLoadData        = false;  % load data from disk
ifPrintFig        = true; % print figures to file

% Plotting/movie options
ifPlotExpStd = false; % running forecast plot of expectation and standard dev.
ifPlotPrb    = true; % running probability forecast
ifPlotErr   = true; % plot forecast errors
ifPlotEnt   = false; % entropy plots 
ifMovieProb = false; % histogram movie
ifShowObs   = false; % show observations in plots


%% NLSA PARAMETERS
% ERSSTv4 reanalysis
% NLSA.dataset           = 'ersstV4';
%NLSA.trainingPeriod    = { '197001' '202002' };
%NLSA.climatologyPeriod = { '198101' '201012' };
%NLSA.sourceVar         = { 'IPSST' };
%NLSA.targetVar         = { 'Nino3.4' };
%NLSA.embWindow         = 48; 
%NLSA.kernel            = 'cone';
%NLSA.ifDen             = false;

% CCSM4 control run
NLSA.dataset           = 'ccsm4Ctrl';
NLSA.trainingPeriod    = { '000101' '109912' };
NLSA.testPeriod        = { '110001' '130012' };
NLSA.climatologyPeriod = NLSA.trainingPeriod;
% NLSA.sourceVar         = { 'IPSST' };
NLSA.sourceVar         = { 'Nino3.4' };
NLSA.obsVar            = { 'IPSST' };
NLSA.targetVar         = { 'Nino3.4' 'Nino4' 'Nino3' 'Nino1+2' };
% NLSA.embWindow         = 9;
NLSA.embWindow         = 11; 
NLSA.embWindow         = 23; 
% NLSA.embWindow         = 12; 
NLSA.kernel            = 'l2';
NLSA.ifDen             = true;
% NLSA.embWindow         = 13; 
% NLSA.kernel            = 'cone';
% NLSA.ifDen             = false;

experiment = experimentStr( NLSA );

%% SETUP GLOBAL PARAMETERS 
switch experiment

case 'ccsm4Ctrl_000101-109912_110001-130012_IPSST_IPSST_emb9_l2_den'
    idxY  = [ 1 ]; % predicted  components (Nino 3.4 index) 
    idxR  = 1;     % realization (ensemble member) in assimilation phase   
    nL    = 3000;  % eigenfunctions used for operator approximation
    % nL    = 4000;  % eigenfunctions used for operator approximation
    nLObs = 3000;   % eigenfunctions used as observational feature vectors
    nQ    = 31;     % quantization levels 
    nTO   = 1;     % timesteps between obs
    nTF   = 12;    % number of forecast timesteps (must be at least nTO)
    qObs  = 0.1; 
    epsilonScl = 1; % bandwidth scaling factor 
    
    shape_fun = @bump; % kernel shape function
    % shape_fun = @rbf;  % kernel shape function

    % Number of parallel workers
    NPar.koopmanOp    = 0; % Koopman operator calculation  
    NPar.observableOp = 0; % Observable operator calculation 
    NPar.featureOp    = 0; % Feature operator calculation
    NPar.qmda         = 0; % Main QMDA loop
    NPar.da           = 0; % Main QMDA loop

    Plt.tLim = { '120005' '122005' }; % time limit to plot
    Plt.idxTF = [ 0 : 3 : 12 ] + 1; % lead times for running forecast
    Plt.idxY = 1; % estimated components for running probability forecast
    Plt.yQLim = [-3 3 ];

case 'ccsm4Ctrl_000101-109912_110001-130012_IPSST_IPSST_emb11_l2_den'
    idxY  = [ 1 ]; % predicted  components (Nino 3.4 index) 
    idxR  = 1;     % realization (ensemble member) in assimilation phase   
    nL    = 3000;  % eigenfunctions used for operator approximation
    % nL    = 4000;  % eigenfunctions used for operator approximation
    nLObs = 3000;   % eigenfunctions used as observational feature vectors
    nQ    = 31;     % quantization levels 
    nTO   = 1;     % timesteps between obs
    nTF   = 12;    % number of forecast timesteps (must be at least nTO)
    qObs  = 0.1; 
    epsilonScl = 1; % bandwidth scaling factor 
    
    shape_fun = @bump; % kernel shape function
    % shape_fun = @rbf;  % kernel shape function

    % Number of parallel workers
    NPar.koopmanOp    = 0; % Koopman operator calculation  
    NPar.observableOp = 0; % Observable operator calculation 
    NPar.featureOp    = 0; % Feature operator calculation
    NPar.qmda         = 0; % Main QMDA loop
    NPar.da           = 0; % Main QMDA loop

    Plt.tLim = { '120005' '122005' }; % time limit to plot
    Plt.idxTF = [ 0 : 3 : 12 ] + 1; % lead times for running forecast
    Plt.idxY = 1; % estimated components for running probability forecast
    Plt.yQLim = [-3 3 ];

case 'ccsm4Ctrl_000101-109912_110001-130012_IPSST_IPSST_emb12_l2_den'
    idxY  = [ 1 ]; % predicted  components (Nino 3.4 index) 
    idxR  = 1;     % realization (ensemble member) in assimilation phase   
    nL    = 3000;  % eigenfunctions used for operator approximation
    % nL    = 4000;  % eigenfunctions used for operator approximation
    nLObs = 3000;   % eigenfunctions used as observational feature vectors
    nQ    = 31;     % quantization levels 
    nTO   = 1;     % timesteps between obs
    nTF   = 12;    % number of forecast timesteps (must be at least nTO)
    qObs  = 0.1; 
    epsilonScl = 1; % bandwidth scaling factor 
    
    shape_fun = @bump; % kernel shape function
    % shape_fun = @rbf;  % kernel shape function

    % Number of parallel workers
    NPar.koopmanOp    = 0; % Koopman operator calculation  
    NPar.observableOp = 0; % Observable operator calculation 
    NPar.featureOp    = 0; % Feature operator calculation
    NPar.qmda         = 0; % Main QMDA loop
    NPar.da           = 0; % Main QMDA loop

    Plt.tLim = { '120005' '122005' }; % time limit to plot
    Plt.idxTF = [ 0 : 3 : 12 ] + 1; % lead times for running forecast
    Plt.idxY = 1; % estimated components for running probability forecast
    Plt.yQLim = [-3 3 ];

case 'ccsm4Ctrl_000101-109912_110001-130012_IPSST_IPSST_emb13_cone'
    idxY  = [ 1 ]; % predicted  components (Nino 3.4 index) 
    idxR  = 1;     % realization (ensemble member) in assimilation phase   
    nL    = 3000;  % eigenfunctions used for operator approximation
    nLObs = 3000;   % eigenfunctions used as observational feature vectors
    nQ    = 31;     % quantization levels 
    nTO   = 1;     % timesteps between obs
    nTF   = 12;    % number of forecast timesteps (must be at least nTO)
    qObs  = 0.1; 
    epsilonScl = 1; % bandwidth scaling factor 
    
    % shape_fun = @bump; % kernel shape function
    shape_fun = @rbf;  % kernel shape function

    % Number of parallel workers
    NPar.koopmanOp    = 0; % Koopman operator calculation  
    NPar.observableOp = 0; % Observable operator calculation 
    NPar.featureOp    = 0; % Feature operator calculation
    NPar.qmda         = 0; % Main QMDA loop

    Plt.tLim = { '120005' '122005' }; % time limit to plot
    Plt.idxTF = [ 0 : 3 : 12 ] + 1; % lead times for running forecast
    Plt.idxY = 1; % estimated components for running probability forecast

case 'ccsm4Ctrl_000101-109912_110001-130012_Nino3.4_IPSST_emb11_l2_den'
    idxY  = [ 1 ]; % predicted  components (Nino 3.4 index) 
    idxR  = 1;     % realization (ensemble member) in assimilation phase   
    nL    = 1500;  % eigenfunctions used for operator approximation
    % nL    = 3000;  % eigenfunctions used for operator approximation
    nLObs = 3000;   % eigenfunctions used as observational feature vectors
    nQ    = 31;     % quantization levels 
    nTO   = 1;     % timesteps between obs
    nTF   = 12;    % number of forecast timesteps (must be at least nTO)
    qObs  = 0.1; 
    epsilonScl = 1; % bandwidth scaling factor 
    
    shape_fun = @bump; % kernel shape function
    % shape_fun = @rbf;  % kernel shape function

    % Number of parallel workers
    NPar.koopmanOp    = 0; % Koopman operator calculation  
    NPar.observableOp = 0; % Observable operator calculation 
    NPar.featureOp    = 0; % Feature operator calculation
    NPar.qmda         = 0; % Main QMDA loop

    Plt.tLim = { '120005' '122005' }; % time limit to plot
    Plt.idxTF = [ 0 : 3 : 12 ] + 1; % lead times for running forecast
    Plt.idxY = 1; % estimated components for running probability forecast
    Plt.yQLim = [-3 3 ];

case 'ccsm4Ctrl_000101-109912_110001-130012_Nino3.4_IPSST_emb23_l2_den'
    idxY  = [ 1 ]; % predicted  components (Nino 3.4 index) 
    idxR  = 1;     % realization (ensemble member) in assimilation phase   
    % nL    = 1500;  % eigenfunctions used for operator approximation
    nL    = 3000;  % eigenfunctions used for operator approximation
    nLObs = 3000;   % eigenfunctions used as observational feature vectors
    nQ    = 31;     % quantization levels 
    nTO   = 1;     % timesteps between obs
    nTF   = 12;    % number of forecast timesteps (must be at least nTO)
    qObs  = 0.1; 
    epsilonScl = 1; % bandwidth scaling factor 
    
    shape_fun = @bump; % kernel shape function
    % shape_fun = @rbf;  % kernel shape function

    % Number of parallel workers
    NPar.koopmanOp    = 0; % Koopman operator calculation  
    NPar.observableOp = 0; % Observable operator calculation 
    NPar.featureOp    = 0; % Feature operator calculation
    NPar.qmda         = 0; % Main QMDA loop

    Plt.tLim = { '120005' '122005' }; % time limit to plot
    Plt.idxTF = [ 0 : 3 : 12 ] + 1; % lead times for running forecast
    Plt.idxY = 1; % estimated components for running probability forecast
    Plt.yQLim = [-3 3 ];

otherwise
    error( [ 'Invalid experiment '  experiment ] )
end

disp( [ 'EXPERIMENT: ' experiment ] )

%% BATCH PROCESSING
iProc = 1; % index of batch process for this script
nProc = 1; % number of batch processes

%% EXTRACT SOURCE DATA
if ifDataSrc
    for iVar = 1 : numel( NLSA.sourceVar )
        msgStr = sprintf( 'Reading source data %s...', ...
                          NLSA.sourceVar{ iVar } );
        disp( msgStr ) 
        t = tic;
        importData( NLSA.dataset, NLSA.trainingPeriod, ...
                    NLSA.sourceVar{ iVar }, NLSA.climatologyPeriod ) 
        toc( t )
    end
end

%% EXTRACT TARGET DATA
if ifDataTrg
    for iVar = 1 : numel( NLSA.targetVar )
        msgStr = sprintf( 'Reading target data %s...', ...
                           NLSA.targetVar{ iVar } );
        disp( msgStr ) 
        t = tic;
        importData( NLSA.dataset, NLSA.trainingPeriod, ...
                    NLSA.targetVar{ iVar }, NLSA.climatologyPeriod ) 
        toc( t )
    end
end

%% EXTRACT OBSERVED DATA
if ifDataObs
    for iVar = 1 : numel( NLSA.obsVar )
        msgStr = sprintf( 'Reading observed data %s...', ...
                          NLSA.obsVar{ iVar } );
        disp( msgStr ) 
        t = tic;
        importData( NLSA.dataset, NLSA.trainingPeriod, ...
                    NLSA.obsVar{ iVar }, NLSA.climatologyPeriod ) 
        toc( t )
    end
end

%% EXTRACT OUT-OF-SAMPLE OBSERVED DATA
if ifDataOutObs
    for iVar = 1 : numel( NLSA.obsVar )
        msgStr = sprintf( 'Reading out-of-sample observed data %s...', ...
                          NLSA.obsVar{ iVar } );
        disp( msgStr ) 
        t = tic;
        importData( NLSA.dataset, NLSA.testPeriod, ...
                    NLSA.obsVar{ iVar }, NLSA.climatologyPeriod ) 
        toc( t )
    end
end

%% EXTRACT OUT-OF-SAMPLE TARGET DATA
if ifDataOutTrg
    for iVar = 1 : numel( NLSA.targetVar )
        msgStr = sprintf( 'Reading out-of-sample target data %s...', ...
                           NLSA.targetVar{ iVar } );
        disp( msgStr ) 
        t = tic;
        importData( NLSA.dataset, NLSA.testPeriod, ...
                    NLSA.targetVar{ iVar }, NLSA.climatologyPeriod ) 
        toc( t )
    end
end

%% BUILD NLSA MODELS, DETERMINE BASIC ARRAY SIZES
% In and InObs are data structures containing the NLSA parameters for the 
% source and observed data, respectively, in the training phase. OutObs is a
% data structure containing the NLSA parameters for the observed data in the 
% data assimilation (out-of-sample) phase.
%
% nY is the number of response variables (Nino indices).
%
% nSE is the number of training samples available after Takens delay
% embedding.
%
% nSO is the number of verification samples available after Takens delay
% embedding. 
%
% nE is equal to half of the delay embedding window length employed in the 
% training data. 
%
% idxT1 is the initial time stamp in the training data ("origin") where delay 
% embedding is performed.
%
% idxT1Obs is the initial time stamp in the observed data ("origin") where delay
% embedding is performed.
%
% nShiftTakens is the temporal shift applied to align the data from the 
% observations NLSA model (modelObs) with the center of the Takens embedding 
% window eployed in the source model (model).  We compute nShiftTakens using
% nE, idxT1 and idxTObs.
%
% idxTObs are the time indices in the test phase at which observations are
% made.
%
% nDA is the number of data assimilation cycles, i.e., number of elements of
% idxTObs.

[ model, In, modelObs, InObs, OutObs ] = ensoQMDA_nlsaModel( experiment );

nY  = numel( idxY );
nSE = getNTotalEmbSample( model );    
nSB = getNXB( model.embComponent );
nSO = getNTotalOutEmbSample( modelObs );

nE           = floor( ( getEmbeddingWindow( model.embComponent ) - 1 ) / 2 );
idxT1        = getOrigin( model.trgEmbComponent( 1 ) ); 
idxT1Obs     = getOrigin( modelObs.trgEmbComponent( 1 ) ); 
nShiftTakens = idxT1 - nE - idxT1Obs;
% nShiftTakens = 8; % Temporary change for testing
% nShiftTakens = 0; % Temporary change for testing


idxTObs = 1 : nTO : nSO;
nDA     = numel( idxTObs ); 

if nShiftTakens < 0
    error( [ 'Negative Takens delay embedding window shift. ' ...
             'Increase the delay embedding origin of the source data, ' ...
             'or decerase the delay embedding origin of the observed data.' ] )
end

%% PARALLEL POOL
% Create parallel pool if running NLSA and the NLSA model has been set up
% with parallel workers. This part can be commented out if no parts of the
% NLSA code utilizing parallel workers are being executed. 
%
% In.nParE is the number of parallel workers for delay-embedded distances
% In.nParNN is the number of parallel workers for nearest neighbor search
if ifNLSA
    if isfield( In, 'nParE' ) && In.nParE > 0
        nPar = In.nParE;
    else
        nPar = 0;
    end
    if isfield( In, 'nParNN' ) && In.nParNN > 0
        nPar = max( nPar, In.nParNN );
    end
    if nPar > 0
        poolObj = gcp( 'nocreate' );
        if isempty( poolObj )
            poolObj = parpool( nPar );
        end
    end
end


%% OUTPUT DIRECTORY
outDir = fullfile( pwd, 'figs', experiment );
if ifPrintFig && ~isdir( outDir )
    mkdir( outDir )
end



%% PERFORM NLSA FOR SOURCE DATA
% Output from each step is saved on disk.
if ifNLSA

    disp( 'Takens delay embedding for source data...' ); t = tic; 
    computeDelayEmbedding( model )
    toc( t )

    % The following step is needed only if we are using velocity-dependent
    % kernels.
    if isa( model.embComponent, 'nlsaEmbeddedComponent_xi' )
      disp( 'Phase space velocity (time tendency of data)...' ); t = tic; 
      computeVelocity( model )
      toc( t )
    end

    % The following steps are needed only if we are using variable-bandwidth
    % kernels.
    if isa( model, 'nlsaModel_den' )
        fprintf( 'Pairwise distances for density data, %i/%i...\n', ...
                  iProc, nProc ); 
       t = tic;
       computeDenPairwiseDistances( model, iProc, nProc )
       toc( t )

       disp( 'Distance normalization for kernel density steimation...' );
       t = tic;
       computeDenBandwidthNormalization( model );
       toc( t )

       disp( 'Kernel bandwidth tuning for density estimation...' ); t = tic;
       computeDenKernelDoubleSum( model );
       toc( t )

       disp( 'Kernel density estimation...' ); t = tic;
       computeDensity( model );
       toc( t )

       disp( 'Takens delay embedding for density data...' ); t = tic;
       computeDensityDelayEmbedding( model );
       toc( t )
    end

    fprintf( 'Pairwise distances (%i/%i)...\n', iProc, nProc ); t = tic;
    computePairwiseDistances( model, iProc, nProc )
    toc( t )

    disp( 'Distance symmetrization...' ); t = tic;
    symmetrizeDistances( model )
    toc( t )

    if isa( model.diffOp, 'nlsaDiffusionOperator_gl_mb' ) ...
            || isa( model.diffOp, 'nlsaDiffusionOperator_gl_mb_bs' )
        disp( 'Kernel bandwidth tuning...' ); t = tic;
        computeKernelDoubleSum( model )
        toc( t )
    end

    disp( 'Kernel eigenfunctions...' ); t = tic;
    computeDiffusionEigenfunctions( model )
    toc( t )
end

%% PERFORM NLSA FOR OBSERVED DATA
if ifNLSAObs

    % Execute NLSA steps. Output from each step is saved on disk.

    disp( 'Takens delay embedding for observed data...' ); t = tic; 
    computeDelayEmbedding( modelObs )
    toc( t )

    disp( 'Takens delay embedding for target data...' ); t = tic; 
    computeTrgDelayEmbedding( modelObs )
    toc( t )
    
    % The following step is needed only if we are using velocity-dependent
    % kernels.
    if isa( modelObs.embComponent, 'nlsaEmbeddedComponent_xi' )
        disp( 'Phase space velocity (time tendency of data)...' ); t = tic; 
        computeVelocity( modelObs )
        toc( t )
    end

    % The following steps are needed only if we are using variable-bandwidth
    % kernels.
    if isa( modelObs, 'nlsaModel_den' )
        fprintf( 'Pairwise distances for density data, %i/%i...\n', ...
                  iProc, nProc ); 
        t = tic;
        computeDenPairwiseDistances( modelObs, iProc, nProc )
        toc( t )

        disp( 'Distance normalization for kernel density steimation...' );
        t = tic;
        computeDenBandwidthNormalization( modelObs );
        toc( t )

        disp( 'Kernel bandwidth tuning for density estimation...' ); t = tic;
        computeDenKernelDoubleSum( modelObs );
        toc( t )

        disp( 'Kernel density estimation...' ); t = tic;
        computeDensity( modelObs );
        toc( t )

        disp( 'Takens delay embedding for density data...' ); t = tic;
        computeDensityDelayEmbedding( modelObs );
        toc( t )
    end

    fprintf( 'Pairwise distances (%i/%i)...\n', iProc, nProc ); t = tic;
    computePairwiseDistances( modelObs, iProc, nProc )
    toc( t )

    disp( 'Distance symmetrization...' ); t = tic;
    symmetrizeDistances( modelObs )
    toc( t )

    if isa( modelObs.diffOp, 'nlsaDiffusionOperator_gl_mb' ) ...
          || isa( modelObs.diffOp, 'nlsaDiffusionOperator_gl_mb_bs' )
        disp( 'Kernel bandwidth tuning...' ); t = tic;
        computeKernelDoubleSum( modelObs )
        toc( t )
    end

    disp( 'Kernel eigenfunctions...' ); t = tic;
    computeDiffusionEigenfunctions( modelObs )
    toc( t )
end

%% DO OUT-OF-SAMPLE EXTENSION FOR OBSERVED DATA
if ifNLSAOut
    
    disp( 'Takens delay embedding for out-of-sample observed data...' ); 
    t = tic;
    computeOutDelayEmbedding( modelObs )
    toc( t )

    disp( 'Takens delay embedding for out-of-sample target data...' ); t = tic; 
    computeOutTrgDelayEmbedding( modelObs )
    toc( t )

    % The following step is needed only if we are using velocity-dependent
    % kernels.
    if isa( modelObs.outEmbComponent, 'nlsaEmbeddedComponent_xi' )
        disp( 'Phase space velocity for out-of-sample data...' ); t = tic; 
        computeOutVelocity( modelObs )
        toc( t )
    end

    % The following steps are needed only if we are using variable-bandwidth
    % kernels.
    if isa( modelObs, 'nlsaModel_den' )
        fprintf( 'OSE pairwise distances for density data... %i/%i\n', ...
                  iProc, nProc ); 
        t = tic;
        computeOseDenPairwiseDistances( modelObs, iProc, nProc )
        toc( t )

        disp( 'OSE density bandwidth normalization...' ); t = tic;
        computeOseDenBandwidthNormalization( modelObs )
        toc( t )

        disp( 'OSE kernel density estimation...' ); t = tic;
        computeOseDensity( modelObs )
        toc( t )

        disp( 'OSE density delay embedding...' ); t = tic; 
        computeOseDensityDelayEmbedding( modelObs );
        toc( t )
    end

    fprintf( 'OSE pairwise distances... %i/%i\n', iProc, nProc ); t = tic; 
    computeOsePairwiseDistances( modelObs, iProc, nProc )
    toc( t )

    disp( 'OSE kernel normalization...' ); t = tic;
    computeOseKernelNormalization( modelObs )

    disp( 'OSE kernel degree...' ); t = tic;
    computeOseKernelDegree( modelObs )
    toc( t )

    disp( 'Nystrom operator...' ); t = tic;
    computeOseDiffusionOperator( modelObs )
    toc( t )

    disp( 'OSE diffusion eigenfunctions...' ); t = tic;
    computeOseDiffusionEigenfunctions( modelObs )
    toc( t )
end


%% READ TRAINING DATA FOR DATA ASSIMILATION
if ifDATrainingData
    disp( 'Retrieving data assimilation training data...' ); 
    t = tic;
    
    % Eigenfunctions from source data 
    [ phi, mu, lambda ] = getDiffusionEigenfunctions( model ); 

    % Eigenfunctions from observed data
    [ phiObs, muObs, lambdaObs ] = getDiffusionEigenfunctions( modelObs );

    % Target data. yL2 is the empirical L2 norm. 
    y = getData( modelObs.trgEmbComponent, [], [], idxY ); 
    yL2 = std( y, 0, 2 );

    % Align phiObs and y with the center of the embedding window
    iStart = nShiftTakens + 1;
    iEnd   = iStart + nSE - 1;
    phiObs = phiObs( iStart : iEnd, : );
    muObs  = muObs( iStart : iEnd ); 
    y      = y( :, iStart : iEnd );
    toc( t )
end

%% READ TEST DATA FOR DATA ASSIMILATION
if ifDATestData
    disp( 'Retrieving data assimilation test data...' ); 
    t = tic;

    % Eigenfunctions 
    phiOut = getOseDiffusionEigenfunctions( modelObs ); 

    % Test (out-of-sample) target data
    yOut = getData( modelObs.outTrgEmbComponent, [], [], idxY ); 
    yObs = yOut( :, idxTObs );

    % Initialization/verification timestamps
    iStart  = idxT1Obs;
    iEnd    = idxT1Obs + nSO - 1;
    tNumOut = getOutTime( modelObs ); % serial time numbers for true signal
    tNumOut = tNumOut( iStart : iEnd ); % shift by delay window 
    tNumObs = datemnth( tNumOut( 1 ), idxTObs - 1 )'; % initialization times  
    tNumVer = repmat( tNumObs, [ 1 nTF + 1 ] ); % verification times
    tNumVer = datemnth( tNumVer, repmat( 0 : nTF, [ nDA, 1 ] ) )';

    % Forecast lead times
    tF = 0 : nTF; 

    toc( t )
end
    
%% KOOPMAN OPERATORS
if ifKoopmanOp
    disp( sprintf( 'Computing Koopman operators for %i timesteps...', nTF ) )
    t = tic;
    U = koopmanOperator( 1 : nTF, phi( :, 1 : nL ), mu, NPar.koopmanOp );
    toc( t )
end

%% MULTIPLICATION OPERATORS FOR TARGET OBSERVABLES
if ifObservableOp
    str = sprintf( [ 'Computing multiplication operators for ' ...
                     '%i observables...' ], nY );  
    disp( str )
    t = tic;

    % Compute multiplication operators representing y and y^2
    M  = multiplicationOperator( y', phi( :, 1 : nL ), mu, ...
                                 NPar.observableOp );
    % M2 = multiplicationOperator( (y .^ 2)', phi( :, 1 : nL ), mu, ...
    %                              NPar.observableOp );
    M2 = M ^ 2;
    
    if nQ > 0 
        b    = linspace( 0, 1, nQ + 1 ); % CDF bins
        yQ = zeros( nQ + 1, nY );
        specIdx = zeros( nQ, 2, nY ); % to rename idxQ
        MEval = zeros( nL, nY );
        MEvec = zeros( nL, nL, nY );
        parfor( iY = 1 : nY, NPar.observableOp )

            % Obtain quantization batch limits as uniform partition of the 
            % empirical inverse CDF of observably iY.
            yQ( :, iY ) = ksdensity( y( iY, : ), b, 'Function', 'icdf' );

            % Compute the spectrum of the multiplication operator representing
            % observable iY.
            [ evecM, evalM ] = eig( M( :, :, iY ), 'vector' );
            [ MEval( :, iY ), idx ] = sort( evalM, 'ascend' ); 
            MEvec( :, :, iY ) = evecM( :, idx, iY );

            % Partition spectrum into batches
            inds = discretize( MEval( :, iY ), yQ( :, iY ) ); 
            specIdxLoc = zeros( nQ, 2 );
            for iQ = 1 : nQ
                specIdxLoc( iQ, 1 ) = find( inds == iQ, 1, 'first' );  
                specIdxLoc( iQ, 2 ) = find( inds == iQ, 1, 'last' );  
            end
            specIdx( :, :, iY ) = specIdxLoc;
        end

        yQ( 1, : )   = min( y' );
        yQ( end, : ) = max( y' ); 

        dyQ = yQ( 2 : end, : ) - yQ( 1 : end - 1, : );

        % Partition the spectra into nQ quantization batches, and ccmpute the 
        % batch-conditional means of the eigenvalues.
        % specPartition = nlsaPartition( 'nSample', nL, 'nBatch', nQ );
        % specIdx = getBatchLimit( specPartition );
        % yQ = fmap( specPartition, @mean, MEval );
        % dyQ = MEval( specIdx( :, 2 ), : ) - MEval( specIdx( :, 1 ), : );
    end

    toc( t )
end

%% KERNEL TUNING
if ifAutotune
    disp( 'Reading observed data for kernel tuning...' )
    t = tic;
    x = getData( modelObs.embComponent );
    if ifVB
        [ ~, rhoInfo ] = computeOptimalBandwidth( modelObs.density );
        rho = getDensity( modelObs.density ) .^ (- 1 / rhoInfo.dEst );
    end
    toc( t )

    disp( 'Computing pairwise distances for kernel tuning...' )
    t = tic;
    d2 = dmat( x );
    if ifVB
        d2 = rho .\ dmat( x ) ./ rho';
    end
    toc( t )

    disp( 'Performing kernel tuning...' )
    t = tic ;
    [ epsilonOpt, epsilonInfo ] = autotune( d2( : ), shape_fun, ...
                                            exponent=[0 20], nBatch=2, n=500);  
    toc( t )
end


%% FEATURE OPERATORS USING KERNEL SHAPE FUNCTION
if ifFeatureOp_sfunc
    disp( 'Reading out-of-sample observed data for feature operators...' )
    t = tic;
    xO = getData( modelObs.outEmbComponent );
    xO = xO(:, idxTObs );
    if ifVB
        rhoO = getDensity( modelObs.oseDensity ) .^ (- 1 / rhoInfo.dEst );
    end
    toc( t )

    str = sprintf( [ 'Computing feature vectors for ' ...
                     '%i out-of-sample states...' ], nSO );  
    disp( str )
    t = tic;
    iStart = nShiftTakens + 1;
    iEnd   = iStart + nSE - 1;
    if ifVB
        if ifSqrtm
            k = shape_fun( ...
                rho .\ dmat( x, xO ) ./ rhoO' ...
                / ( epsilonOpt * epsilonScl ) ^ 2 ) ;
        else
            k = sqrt( ...
                shape_fun( ...
                rho .\ dmat( x, xO ) ./ rhoO' ...
                / ( epsilonOpt * epsilonScl ) ^ 2 ) );
        end
        k = k( iStart : iEnd, : );
    else
        k = sqrt( ...
            shape_fun( ...
            dmat( x( :, iStart : iEnd ), xO ) ...
            / ( epsilonOpt * epsilonScl ) ^ 2 ) );
    end
    if any( all( k == 0, 1) )
        error( 'Zero feature vectors detected.' )
    end
    toc( t )

    disp( 'Computing feature operators...' )
    t = tic;
    if ifSqrtm
        K = sqrtMultiplicationOperator( k, phi( :, 1 : nL ), mu, ...
                                        NPar.featureOp );
    else
        K = multiplicationOperator( k, phi( :, 1 : nL ), mu, NPar.featureOp );
    end
    toc( t )
end

%% FEATURE OPERATORS USING KERNEL EIGENFUNCTIONS 
if ifFeatureOp
    str = sprintf( [ 'Computing feature operators for ' ...
                     '%i out-of-sample states...' ], nSO );  
    disp( str )
    t = tic;
    k = featureVector( phiOut( idxTObs, 1 : nLObs ), phiObs( :, 1 : nLObs ), ...
                       ( lambdaObs( 1 : nLObs ) .^ qObs )', NPar.featureOp );
    K = multiplicationOperator( k, phi( :, 1 : nL ), mu, NPar.featureOp );
    toc( t )
end



%% DATA ASSIMILATION
if ifDA 
    % xi0 is the initial wavefunction. We set xi0 to the constant wavefunction 
    % representing the invariant distrubution of the system, i.e., only the 
    % leading spectral expansion coefficient 
    xi0      = zeros( 1, nL );
    xi0( 1 ) = 1;

    if nQ > 0 
        [ yExp, yStd, yPrb ] = qmda( K, U, M, M2, MEvec, specIdx, xi0, ...
                                     nTO, NPar.qmda );
    else
        [ yExp, yStd ] = qmda( K, U, M, M2, [], [], xi0, nTO, NPar.qmda );
    end
    yExp0 = yExp( :, :, 1 );
    yStd0 = yExp( :, :, 1 );
    yExp = yExp( :, :, 2 : end );
    yStd = yStd( :, :, 2 : end );
    if nQ > 0
        yPrb0 = yPrb( :, :, :, 1 ); % probability mass
        yPrb = yPrb( :, :, :, 2 : end ); % probability mass
        yDen0 = yPrb0 ./ dyQ; 
        yDen = yPrb ./ dyQ;              % probability density 
    end

    if ifSaveData
        dataFile = fullfile( outDir, 'forecast.mat' );
        outVars = { 'tNumVer' 'tNumObs' 'tNumOut' 'yOut' ...
                    'yExp' 'yStd' 'yExp0' 'yStd0' };  
        if nQ > 0
            outVars = [ outVars { 'yQ' 'yPrb' 'yDen' 'yPrb0' 'yDen0' } ];
        end
        save( dataFile, outVars{ : }, '-v7.3' )
    end
end

%% DATA ASSIMILATION ERRORS
% yRMSE is an array of size [ nY, nTF + 1 ] such that yRMSE( i, j ) is the 
% root mean square error of the mean forecast yExp for target variable 
% idxY( i ) at lead time j - 1. 
%
% yPC is as yRMSE but stores the anomaly correlation coefficient between the
% true signal and mean forecast.
%
% yRMSE_est is the estimated forecast error based on the predicted standard 
% deviation yStd.
if ifDAErr
    disp( 'Computing forecast skill scores...' ); 
    t = tic;
    [ ~, yRMSE, yPC ] = forecastError( yOut, yExp );
    yRMSE = yRMSE ./ yL2;
    yRMSE_est = sqrt( mean( yStd .^ 2, 3 ) ) ./ yL2;
    toc(t)

    if ifSaveData
        dataFile = fullfile( outDir, 'forecast_error.mat' );
        save( dataFile, 'tF', 'yRMSE', 'yRMSE_est', 'yPC' )
    end
end

%% RUNNING FORECAST PLOTS OF EXPECTATION AND STANDARD DEVIATION
if ifPlotExpStd

    nTPlt = numel( Plt.idxTF ); % number of lead times to plot
    nYPlt = numel( Plt.idxY ); % number of target variables to plot
    nTickSkip = 12;           % number of months between time axis ticks

    % Figure parameters
    Fig.nTileX     = nYPlt;
    Fig.nTileY     = nTPlt;
    Fig.units      = 'inches';
    Fig.figWidth   = 9;    
    Fig.deltaX     = .55;
    Fig.deltaX2    = .7;
    Fig.deltaY     = .5;
    Fig.deltaY2    = .3;
    Fig.gapX       = .3;
    Fig.gapY       = .6;
    Fig.gapT       = 0.32; 
    Fig.aspectR    = ( 3 / 4 ) ^ 5;
    Fig.fontName   = 'helvetica';
    Fig.fontSize   = 10;
    Fig.tickLength = [ 0.01 0 ];
    Fig.visible    = 'on';
    Fig.nextPlot   = 'add'; 

    % Set up figure and axes
    [ fig, ax, axTitle ] = tileAxes( Fig );

    % Loop over target variables
    for iY = 1 : nYPlt

        % Loop over lead times
        for iT = 1 : nTPlt 

            % Serial date numbers for verification times
            tPlt      = tNumVer( ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : );
            tPlt      = tPlt( : )';
            ifTruePlt = tNumOut >= tPlt( 1 ) & tNumOut <= tPlt( end );
            ifObsPlt  = tNumObs >= tPlt( 1 ) & tNumObs <= tPlt( end ); 
            tLabels   = cellstr( datestr( tPlt , 'mm/yy' ) ); 
        
            % Assemble plot data for predicted mean and standard deviation
            yExpPlt = yExp( Plt.idxY( iY ), ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : );
            yExpPlt = yExpPlt( : );
            yStdPlt = yStd( Plt.idxY( iY ), ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : );
            yStdPlt = yStdPlt( : );

            set( gcf, 'currentAxes', ax( iY, iT ) )

            % Plot predicted mean and standard deviation
            [ hExp, hErr ] = boundedline( tPlt, yExpPlt, yStdPlt, 'b-' );

            % Plot true signal 
            plot( tNumOut( ifTruePlt ), ...
                  yOut( Plt.idxY( iY ), find( ifTruePlt ) ), ...
                  'r-' )
            % Plot observations
            plot( tNumObs( ifObsPlt ), ...
                  yObs( Plt.idxY( iY ), find( ifObsPlt ) ), ...
                  'r*', 'linewidth', 2 )

            grid on
            set( gca, 'xTick', tPlt( 1 : nTickSkip : end ), ...
                      'xTickLabel', tLabels( 1 : nTickSkip : end' ), ...
                      'xLimSpec', 'tight' )
            ylim( [ -3 3 ] )
            ylabel( NLSA.targetVar( idxY( Plt.idxY( iY ) ) ) )
            title( sprintf( 'Lead time = %i months', Plt.idxTF( iT ) - 1 ) )
            if iT == nTPlt
                xlabel( 'verification time' )
            end
        end
    end
end


%% RUNNING PROBABILITY FORECAST 
if ifPlotPrb

    if ifLoadData
        dataFile = fullfile( outDir, 'forecast.mat' );
        disp( sprintf( 'Loading forecast data from file %s...', ...
              dataFile ) )
        load( dataFile, 'tNumVer', 'tNumObs', 'tNumOut', 'yOut', 'yExp', ...
                        'yStd', 'yQ', 'yDen' ) 
    end

    nTPlt = numel( Plt.idxTF ); % number of lead times to plot
    nYPlt = numel( Plt.idxY ); % number of target variables to plot
    nTickSkip = 12;           % number of months between time axis ticks

    % Determine time indices to plot within verification time interval
    tNumLimPlt = datenum( Plt.tLim, 'yyyymm' );
    [idxTLimPlt, ~] = find( tNumObs == tNumLimPlt', 2 );
    idxTPlt = idxTLimPlt( 1 ) : idxTLimPlt( 2 );


    % Figure parameters
    Fig.nTileX     = nYPlt;
    Fig.nTileY     = nTPlt;
    Fig.units      = 'inches';
    Fig.figWidth   = 9;    
    Fig.deltaX     = .55;
    Fig.deltaX2    = .9;
    Fig.deltaY     = .75;
    Fig.deltaY2    = .25;
    Fig.gapX       = .2;
    Fig.gapY       = .7;
    Fig.gapT       = .25; 
    Fig.aspectR    = ( 3 / 4 ) ^ 5;
    Fig.fontName   = 'helvetica';
    Fig.fontSize   = 10;
    Fig.tickLength = [ 0.007 0 ];
    Fig.visible    = 'on';
    Fig.nextPlot   = 'add'; 

    % Set up figure and axes
    [ fig, ax, axTitle ] = tileAxes( Fig );

    % Loop over target variables
    for iY = 1 : nYPlt

        % Loop over lead times
        for iT = 1 : nTPlt 

            iTF = (1 : nTO) + Plt.idxTF( iT ) - 1;

            % Serial date numbers and labels for verification times
            tPlt      = tNumVer( iTF, idxTPlt );
            tPlt      = tPlt( : )';
            ifTruePlt = tNumOut >= tPlt( 1 ) & tNumOut <= tPlt( end );
            ifObsPlt  = tNumObs >= tPlt( 1 ) & tNumObs <= tPlt( end ); 
            tLabels   = cellstr( datestr( tPlt , 'mm/yy' ) ); 
        
            % Form grid of forecast times and bin centers. We add "ghost" 
            % points in the end to ensure that everything gets plotted using
            % pcolor.

            tGrd = [ tPlt tPlt( end ) + 1 ];
            yGrd = yQ( :, iY );
            yGrd( 1 ) = Plt.yQLim( 1 );
            yGrd( end ) = Plt.yQLim( end );

            % yGrd = [ yQ( :, iY ); yQ( end, iY ) + 1 ];

            % Assemble plot data for mean and probability distribution
            yExpPlt = yExp( Plt.idxY( iY ), iTF, idxTPlt );
            yExpPlt = yExpPlt( : );
            yDenPlt = yDen( :, Plt.idxY( iY ), iTF, idxTPlt );
            yDenPlt = reshape( yDenPlt, nQ, numel( tPlt ) );
            yDenPlt = [ yDenPlt; yDenPlt( end, : ) ];
            yDenPlt = [ yDenPlt yDenPlt( :, end ) ];

            set( gcf, 'currentAxes', ax( iY, iT ) )

            % Plot predicted probability density
            % hPrb = pcolor( tGrd, yGrd, yDenPlt ); 
            hPrb = pcolor( tGrd, yGrd, log10( yDenPlt ) ); 
            set( hPrb, 'edgeColor', 'none' )
            %contourf( tGrd, yGrd, log10( yDenPlt ), 7 )

            % Plot predicted mean 
            hExp = plot( tPlt, yExpPlt, 'k-', 'lineWidth', 2 );

            % Plot true signal 
            hTrue = plot( tNumOut( ifTruePlt ), ...
                          yOut( Plt.idxY( iY ), find( ifTruePlt ) ), ...
                          'r-', 'lineWidth', 2 );
            % Plot observations
            if ifShowObs
                plot( tNumObs( ifObsPlt ), ...
                      yObs( Plt.idxY( iY ), find( ifObsPlt ) ), ...
                      'r*', 'linewidth', 2 )
            end

            set( gca, 'xTick', tPlt( 1 : nTickSkip : end ), ...
                      'xTickLabel', tLabels( 1 : nTickSkip : end' ), ...
                      'xLimSpec', 'tight' )
            ylim( [ -2.8 3 ] )
            set( gca, 'cLim', [ -3 0 ] )
            % set( gca, 'cLim', [ 0 5 ] )
            % colormap(gray)
            axPos = get( gca, 'position' );
            hC = colorbar( 'location', 'eastOutside' );
            cPos = get( hC, 'position' );
            cPos( 3 ) = .5 * cPos( 3 );
            cPos( 1 ) = cPos( 1 ) + .08;
            set( gca, 'position', axPos )
            set( hC, 'position', cPos )
            ylabel( hC, 'log_{10}\rho' )
            
            % ylabel( NLSA.targetVar( idxY( Plt.idxY( iY ) ) ) )
            ylabel( 'Nino 3.4' )
            title( sprintf( 'Lead time \\tau = %i months', Plt.idxTF( iT ) - 1 ) )

            if iT == nTPlt
                xlabel( 'verification time' )
            end

            % Add legend
            if iT == 1 
                if ifShowObs
                    hL = legend( [ hTrue hObs hExp ], ...
                                  'true', 'observations', 'mean forecast', ...
                                  'location', 'northWest' );
                else
                    hL = legend( [ hTrue hExp ], ...
                                  'true', 'mean forecast', ...
                                  'location', 'northWest' );
                end
                sL = hL.ItemTokenSize;
                sL( 1 ) = .5 * sL( 1 );
                hL.ItemTokenSize = sL;
                lPos = get( hL, 'position' );
                lPos( 2 ) = lPos( 2 ) + 0.04;
                set( hL, 'position', lPos )
            end
        end
    end

    title( axTitle, 'Running forecasts: CCSM4 ENSO' )

    % Print figure
    if ifPrintFig
        figFile = fullfile( outDir, 'nino34_prob.png' );
        print( fig, figFile, '-dpng', '-r300' ) 
    end
end

%% PLOT ERROR SCORES
if ifPlotErr

    if ifLoadData
        dataFile = fullfile( outDir, 'forecast_error.mat' );
        disp( sprintf( 'Loading forecast error data from file %s...', ...
              dataFile ) )
        load( dataFile, 'tF', 'yRMSE', 'yRMSE_est', 'yPC' )
    end

    % Set up figure and axes 
    Fig.nTileX     = 2;
    Fig.nTileY     = nY;
    Fig.units      = 'inches';
    Fig.figWidth   = 8; 
    Fig.deltaX     = .4;
    Fig.deltaX2    = .4;
    Fig.deltaY     = .5;
    Fig.deltaY2    = .25;
    Fig.gapX       = .3;
    Fig.gapY       = .3;
    Fig.gapT       = 0.05; 
    Fig.aspectR    = 3 / 4;
    Fig.fontName   = 'helvetica';
    Fig.fontSize   = 10;
    Fig.tickLength = [ 0.02 0 ];
    Fig.visible    = 'on';
    Fig.nextPlot   = 'add'; 

    [ fig, ax, axTitle ] = tileAxes( Fig );

    % Plot skill scores
    for iY = 1 : nY 

        % Normalized RMSE
        set( gcf, 'currentAxes', ax( 1, iY ) )
        plot( tF, yRMSE( iY, : ), 'linewidth', 1.5 )
        plot( tF, yRMSE_est( iY, : ), 'r-', 'lineWidth', 1.5 )

        grid on
        ylim( [ 0 1 ] )
        xlim( [ tF( 1 ) tF( end ) ] )
        set( gca, 'xTick', tF( 1 ) : 3 : tF( end ), 'yTick', 0 : .2 : 1.2 )  
        % ylabel( 'Normalized RMSE' )
        legend( 'Normalized RMSE', 'Forecast standard deviation', ...
                'location', 'southEast' )
        if iY == nY
            xlabel( 'Lead time \tau (months)' )
        else
            set( gca, 'xTickLabel', [] )
        end
        text( 0.3, 0.95, '(a)' )
        
        % Anomaly correlation
        set( gcf, 'currentAxes', ax( 2, iY ) )
        plot( tF, yPC( iY, : ), 'linewidth', 1.5 )

        grid on
        ylim( [ 0 1 ] )
        xlim( [ tF( 1 ) tF( end ) ] )
        set( gca, 'xTick', tF( 1 ) : 3 : tF( end ), 'yTick', 0 : .2 : 1.2, ...
                  'yAxisLocation', 'right' )   
        % ylabel( 'Anomaly correlation' )
        legend( 'Anomaly correlation', 'location', 'southWest' )
        if iY == nY
            xlabel( 'Lead time \tau (months)' )
        else
            set( gca, 'xTickLabel', [] )
        end
        text( 11, 0.95, '(b)' )
    end

    % Add figure title
    title( axTitle, 'Forecast skill: CCSM4 ENSO' )

    % Print figure
    if ifPrintFig
        figFile = fullfile( outDir, 'nino34_err.png' );
        print( fig, figFile, '-dpng', '-r300' ) 
    end
end




if ifPlotEnt

    nTPlt = numel( Plt.idxTF );
    nYPlt = numel( Plt.idxY );
    nTickSkip = 12;

    Fig.figWidth   = 6;    % in inches
    Fig.deltaX     = .55;
    Fig.deltaX2    = .85;
    Fig.deltaY     = .5;
    Fig.deltaY2    = .3;
    Fig.gapX       = .3;
    Fig.gapY       = .6;
    Fig.nSkip      = 5;
    Fig.cLimScl    = 1;

    nTileX = nYPlt;
    nTileY = nTPlt;

    panelX = ( Fig.figWidth - Fig.deltaX - Fig.deltaX2 - ( nTileX -1 ) * Fig.gapX ) / nTileX;
    panelY = panelX * ( 3 / 4 )^3;

    posn     = [ 0 ...
                 0 ...
                 nTileX * panelX + ( nTileX - 1 ) * Fig.gapX + Fig.deltaX + Fig.deltaX2 ...
                 nTileY * panelY + ( nTileY - 1 ) * Fig.gapY + Fig.deltaY + Fig.deltaY2 ];


    fig = figure( 'units', 'inches', ...
                  'paperunits', 'inches', ...
                  'position', posn, ...
                  'paperPosition', posn, ...
                  'color', 'white', ...
                  'doubleBuffer', 'on', ...
                  'backingStore', 'off', ...
                  'defaultAxesTickDir', 'out', ...
                  'defaultAxesNextPlot', 'replace', ...
                  'defaultAxesBox', 'on', ...
                  'defaultAxesFontSize', 8, ...
                  'defaultTextFontSize', 8, ...
                  'defaultAxesTickDir',  'out', ...
                  'defaultAxesTickLength', [ 0.01 0 ], ...
                  'defaultAxesFontName', 'helvetica', ...
                  'defaultTextFontName', 'helvetica', ...
                  'defaultAxesLayer', 'top' );

    ax = zeros( nTileX, nTileY );

    for iAx = 1 : nTileX
        for jAx = 1 : nTileY
            ax( iAx, jAx ) = axes( ...
                    'units', 'inches', ...
                    'position', ...
                    [ Fig.deltaX + ( iAx - 1 ) * ( panelX + Fig.gapX ), ...
                      Fig.deltaY + ( nTileY - jAx ) * ( panelY + Fig.gapY ), ...
                      panelX, panelY ] );
        end
    end

    axPos = reshape( get( ax, 'position' ), [ nTileX nTileY ] );

    for iY = 1 : nYPlt
        for iT = 1 : nTPlt

            % verification time serial numbers
            tPlt = tNumVer( ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : );
            tPlt = tPlt( : )';
            ifTruePlt = tPlt <= tNumOut( end ); 
            idxTObsPlt = find( tNumObs >= tPlt( 1 ) & tNumObs <= tPlt( end ) ); 
            tLabels = cellstr( datestr( tPlt , 'mm/yy' ) ); 

            dKLPlt = squeeze( ...
                 dKL( Plt.idxY( iY ), ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : ) ); 
            dKLPlt = dKLPlt( : );

            eKLPlt = squeeze( ...
                 eKL( Plt.idxY( iY ), ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : ) ); 
            eKLPlt = eKLPlt( : );

            set( gcf, 'currentAxes', ax( iY, iT ) )
            plot( tPlt( ifTruePlt ), eKLPlt( ifTruePlt ), 'g-' )
            hold on
            plot( tPlt, dKLPlt, 'b-' ) 
            plot( tPlt, sign( tPlt ) * log2( nQ ), 'm-' )
            axis tight
            %legend( 'E', 'D' ) 
            %legend boxoff
            yLm = get( gca, 'yLim' );
            for iObs = 1 : numel( idxTObsPlt )
                plot( tNumObs( idxTObsPlt( iObs ) ) * [ 1 1 ], yLm, 'r:' )
            end
            title( sprintf( 'Lead time \\tau = %i months', Plt.idxTF( iT ) - 1 ) )
            if iY == 1 
                ylabel( 'Nino 3.4 index' )
            end
            if iT == nTPlt
                xlabel( 'verification time' )
            end
            set( ax, 'tickDir', 'out', 'tickLength', [ 0.005 0 ] )
            set( gca, 'xTick', tPlt( 1 : nTickSkip : end ), 'xTickLabel', tLabels( 1 : nTickSkip : end' ) )
        end
    end

    print -dpng -r300 figNino34Ent.png
end

if ifMovieProb


    nTPlt = numel( Plt.idxTF );
    nYPlt = numel( Plt.idxY );
    
    % bin coordinates
    aPlt = cell( 1, nYPlt );
    daPlt = cell( 1, nYPlt );
    for iY = 1 : nYPlt 
        aPlt{ iY } = aY( Plt.idxY( iY ), : ); 
        aPlt{ iY }( 1 ) = min( y( Plt.idxY( iY ), : ), [], 2 );
        aPlt{ iY }( end ) = max( y( Plt.idxY( iY ), : ), [], 2 );
        daPlt{ iY } = ( aPlt{ iY }( 2 : end ) - aPlt{ iY }( 1 : end - 1 ) )';
    end

    % verification time serial numbers
    tPlt = cell( 1, nTPlt );
    ifTruePlt = cell( 1, nTPlt );
    ifObsPlt = cell( 1, nTPlt );
    tLabels = cell( 1, nTPlt );
    for iT = 1 : nTPlt
        tPlt{ iT } = tNumVer( ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : );
        tPlt{ iT } = tPlt{ iT }( : )';
        ifTruePlt{ iT } = tNumOut >= tPlt{ iT }( 1 ) ...
                        & tNumOut <= tPlt{ iT }( end );
        ifObsPlt{ iT } = tNumObs >= tPlt{ iT }( 1 ) ...
                       & tNumObs <= tPlt{ iT }( end ); 
        tLabels{ iT } = cellstr( datestr( tPlt{ iT } , 'mm/yy' ) ); 
    end
    nFrame = numel( tPlt{ 1 } );
    
    % probability values
    pPlt = cell( nTPlt, nYPlt );
    pLim = zeros( nTPlt, nYPlt );
    for iY = 1 : nYPlt
       for iT = 1 : nTPlt 
            pPlt{ iT, iY } = real( squeeze( ...
               p( :, Plt.idxY( iY ), ( 1 : nTO ) + Plt.idxTF( iT ) - 1, : ) ) ); 
            pPlt{ iT, iY } = reshape( pPlt{ iT, iY }, [ nQ nTO * nDA ] ); 
            pPlt{ iT, iY } = bsxfun( @rdivide, pPlt{ iT, iY }, daPlt{ iY } );
            pLim( iT, iY ) = max( pPlt{ iT, iY }( : ) );
        end
    end
        
    Mov.figWidth   = 300;    % in pixels
    Mov.deltaX     = 35;
    Mov.deltaX2    = 15;
    Mov.deltaY     = 35;
    Mov.deltaY2    = 40;
    Mov.gapX       = 30;
    Mov.gapY       = 40;
    Mov.visible    = 'on';
    Mov.fps        = 12;

    nTileX = nYPlt;
    nTileY = nTPlt;

    panelX = ( Mov.figWidth - Mov.deltaX - Mov.deltaX2 - ( nTileX -1 ) * Mov.gapX ) / nTileX;
    panelY = panelX * 3 / 4;

    posn     = [ 0, ...
                 0, ...
                 nTileX * panelX + ( nTileX - 1 ) * Mov.gapX + Mov.deltaX + Mov.deltaX2, ...
                 nTileY * panelY + ( nTileY - 1 ) * Mov.gapY + Mov.deltaY + Mov.deltaY2 ];

    writerObj = VideoWriter( 'movieNino34_density.avi' );
    writerObj.FrameRate = Mov.fps;
    writerObj.Quality = 100;
    open( writerObj );

    fig = figure( 'units', 'pixels', ...
              'paperunits', 'points', ...
              'position', posn, ...
              'paperPosition', posn, ...
              'visible', Mov.visible, ...
              'color', 'white', ...
              'doubleBuffer', 'on', ...
              'backingStore', 'off', ...
              'defaultAxesTickDir', 'out', ...
              'defaultAxesNextPlot', 'replace', ...
              'defaultAxesBox', 'on', ...
              'defaultAxesFontSize', 8, ...
              'defaultTextFontSize', 8, ...
              'defaultAxesTickDir',  'out', ...
              'defaultAxesTickLength', [ 0.02 0 ], ...
              'defaultAxesFontName', 'helvetica', ...
              'defaultTextFontName', 'helvetica', ...
              'defaultAxesLineWidth', 1, ...
              'defaultAxesLayer', 'top' );

    ax = zeros( nTileX, nTileY );

    for iAx = 1 : nTileX
        for jAx = 1 : nTileY
            ax( iAx, jAx ) = axes( ...
                'units', 'pixels', ...
                'position', ...
                [ Mov.deltaX + ( iAx - 1 ) * ( panelX + Mov.gapX ), ...
                  Mov.deltaY + ( nTileY - jAx ) * ( panelY + Mov.gapY ), ...
                  panelX, panelY ] );
        end
    end

    axTitle = axes( 'units', 'pixels', 'position', [ Mov.deltaX, Mov.deltaY, ...
                              nTileX * panelX + ( nTileX - 1 ) * Mov.gapX, ...
                              nTileY * panelY + ( nTileY - 1 ) * Mov.gapY + 15 ], ...
                    'color', 'none', 'box', 'off' );

    for iFrame = 1 : nFrame
        for iY = 1 : nYPlt 
            for iT = 1 : nTPlt 
                set( gcf, 'currentAxes', ax( iY, iT ) ) 
                if iFrame - Plt.idxTF( iT ) >= 0    
                    histogram( 'binEdges', aPlt{ iY }, ...
                               'binCounts', pPlt{ iT, iY }( :, iFrame - Plt.idxTF( iT ) + 1 ) );
                    hold on
                    if ifTruePlt{ iT }( iFrame )
                        plot( yO( Plt.idxY( iY ), iFrame ) * [ 1 1 ], ...
                              [ 0 1.1 ] * pLim( iT, iY ), 'r-' )
                    end
                    if ifObsPlt{ iT }( iFrame )
                        plot( yO( Plt.idxY( iY ), iFrame ), ...
                              1, 'r*' )
                    end

                    set( gca, 'xLim', [ aPlt{ iY }( 1 ) - 0.1 * ( aPlt{ iY }( end ) - aPlt{ iY }( 1 ) ), ...
                    aPlt{ iY }( end ) + 0.1 * ( aPlt{ iY }( end ) - aPlt{ iY }( 1 ) ) ], ... 
                    'yLim', [ 0 1 ] )  
                else
                    axis off
                end
                title( sprintf( 'Lead time %i months', Plt.idxTF( iT ) - 1 ) ) 
            end
        end
        
        title( axTitle, [ 'Verification time ' tLabels{ 1 }{ iFrame } ] )
        axis( axTitle, 'off' )

        frame = getframe( gcf );
        writeVideo( writerObj, frame )

        for iY = 1 : nYPlt
            for iT = 1 : nTPlt
                cla( ax( iY, iT ), 'reset' )
            end
        end
        cla( axTitle, 'reset' )
    end

    close( writerObj )
end


%% STRING IDENTIFIER FOR DATA ANALYSIS EXPERIMENT
function s = experimentStr( P )
s = strjoin_e( { P.dataset ...
                 strjoin_e( P.trainingPeriod, '-' ) ... 
                 strjoin_e( P.testPeriod, '-' ) ...  
                 strjoin_e( P.sourceVar, '_' ) ...
                 strjoin_e( P.obsVar, '_' ) ...
                 sprintf( 'emb%i', P.embWindow ) ...
                 P.kernel }, ...
               '_' );

if isfield( P, 'ifDen' )
    if P.ifDen
        s = [ s '_den' ];
    end
end
end
