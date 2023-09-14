% KOOPMAN ANALYSIS OF FSD DATA
%
% Available data: 2000 to 2022
% To retrieve NLSA basis functions: phi = getDiffusionEigenfunctions(model);
% To retrieve Koopman eigenfunctions: z = getKoopmanEigenfunctions(model);
%
% Modified 2023-09-14

%% SCRIPT EXECUTION OPTIONS
ifTrainingData          = true; % read training data in appropriate format for NLSA
ifNLSA                  = true; % perform NLSA
ifKoopman               = true; % compute Koopman eigenfunctions
ifKoopmanReconstruction = false; % perform Koopman reconstruction


%% DATA & NLSA PARAMETERS
NLSA.yrLim     = [2000 2022]; % analysis years
NLSA.dayLim    = [65 269]; % analysis days in each year
NLSA.var       = 'alpha'; % variables used for training
NLSA.embWindow = 60; % delay embedding window
NLSA.kernel    = 'l2'; % kernel type
NLSA.den       = true; % set true to use variable-bandwidth kernel

dataFunc  = @importData_alpha; % data import function
modelFunc = @nlsaModel_alpha; % function to build NLSA model

experiment = experimentStr(NLSA);          % string identifier
disp(['NLSA EXPERIMENT: ' experiment])


%% BATCH PROCESSING
iProc = 1; % current process
nProc = 1; % number of processes


%% DATA RECONSTRUCTION
Rec.iC = 1; % component to reconstruct
Rec.iB = 1; % batch to reconstruct


%% READ DATA
if ifTrainingData
    disp(['Reading training data using function ' func2str(dataFunc) '...'])
    t = tic;
    dataFunc(NLSA);
    toc(t)
end


%% BUILD NLSA MODEL
disp('Building NLSA model...')
t = tic;false
model = modelFunc(NLSA);
toc(t)


%% PERFORM NLSA
% Output from each step is saved on disk.
if ifNLSA
    disp('Takens delay embedding for source data...'); t = tic;
    computeDelayEmbedding(model)
    toc(t)

    % The following step is needed only if we are using velocity-dependent
    % kernels.
    if isa(model.embComponent, 'nlsaEmbeddedComponent_xi')
        disp('Phase space velocity (time tendency of data)...'); t = tic;
        computeVelocity(model)
        toc(t)
    end

    % The following step is only needed if query partition is employed for
    % source data.
    if ~isempty(model.embComponentQ)
        disp('Forming query partition for source data...'); t = tic;
        computeDelayEmbeddingQ(model)
        toc(t)
    end

    % The following step is only needed if test partition is employed for source
    % data.
    if ~isempty(model.embComponentT)
        disp('Forming test partition for source data...'); t = tic;
        computeDelayEmbeddingT(model)
        toc(t)
    end

    % The following steps are needed only if we are using variable-bandwidth
    % kernels.
    if isa(model, 'nlsaModel_den')
        fprintf('Pairwise distances for density data, %i/%i...\n', ...
                  iProc, nProc);
        t = tic;
        computeDenPairwiseDistances(model, iProc, nProc)
        toc(t)

        disp('Distance normalization for kernel density estimation...');
        t = tic;
        computeDenBandwidthNormalization(model);
        toc(t)

        disp('Kernel bandwidth tuning for density estimation...'); t = tic;
        computeDenKernelDoubleSum(model);
        toc(t)

        disp('Kernel density estimation...'); t = tic;
        computeDensity(model);
        toc(t)

        % The next step is only needed if a query partition was used for the
        % density data.
        if ~isempty(model.denEmbComponentQ)
            disp('Density splitting...'); t = tic;
            computeDensitySplitting(model);
            toc(t)
        end

        disp('Takens delay embedding for density data...'); t = tic;
        computeDensityDelayEmbedding(model);
        toc(t)
        % The following step is only needed if query partition is employed for
        % density data.
        if ~isempty(model.denEmbComponentQ)
            disp('Forming query partition for density data...'); t = tic;
            computeDensityDelayEmbeddingQ(model)
            toc(t)
        end

        % The following step is only needed if test partition is employed for
        % density data.
        if ~isempty(model.denEmbComponentT)
            disp('Forming test partition for density data...'); t = tic;
            computeDensityDelayEmbeddingT(model)
            toc(t)
        end
    end

    fprintf('Pairwise distances (%i/%i)...\n', iProc, nProc); t = tic;
    computePairwiseDistances(model, iProc, nProc)
    toc(t)

    disp('Distance symmetrization...'); t = tic;
    symmetrizeDistances(model)
    toc(t)

    disp('Kernel bandwidth tuning...'); t = tic;
    computeKernelDoubleSum(model)
    toc(t)

    disp('Kernel eigenfunctions...'); t = tic;
    computeDiffusionEigenfunctions(model)
    toc(t)
end


%% COMPUTE EIGENFUNCTIONS OF KOOPMAN GENERATOR
if ifKoopman
    disp('Koopman eigenfunctions...'); t = tic;
    computeKoopmanEigenfunctions(model, 'ifLeftEigenfunctions', true)
    toc(t)
end


%% PERFORM KOOPMAN RECONSTRUCTION
if ifKoopmanReconstruction

    disp('Takens delay embedding...'); t = tic;
    computeTrgDelayEmbedding(model, Rec.iC)
    toc(t)

    disp('Projection of target data onto Koopman eigenfunctions...'); t = tic;
    computeKoopmanProjection(model, Rec.iC)
    toc(t)

    disp('Reconstruction of target data from Koopman eigenfunctions...');
    t = tic;
    computeKoopmanReconstruction(model, Rec.iB, Rec.iC)
    toc(t)
end
