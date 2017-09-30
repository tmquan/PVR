function y=make_map(data_x,data_y,data_z,dict_x,dict_y,dict_z,dict_n,lam,data_dir,save_dir,original)
%% 
    addpath(genpath('.'));

    xsz = [data_x,data_y,data_z];
    dsz = [ dict_x,dict_y,dict_z];
    %% Seed the randomness
    rng(2016);

    plan.elemSize = [xsz,  1];
    plan.dataSize = [xsz,  1]; % For example
    plan.atomSize = [dsz,  1];
    plan.dictSize = [dsz, dict_n];
    plan.blobSize = [xsz, dict_n];
    plan.iterSize = [xsz, dict_n]; 

    % plan.elemSize = [256, 256,  256,  1];
    % plan.dataSize = [256, 256,  256,  1]; % For example
    % plan.atomSize = [ 11,  11,  11,   1];
    % plan.dictSize = [ 11,  11,  11,  16];
    % plan.blobSize = [256, 256,  256, 16];
    % plan.iterSize = [256, 256,  256, 16]; 

    %% Initialize the plan
    plan.alpha  = params; % See param.m
    plan.gamma  = params; 
    plan.delta  = params;
    plan.theta  = params;
    plan.omega  = params;
    plan.lambda = params; 
    plan.sigma  = params; 
    plan.rho    = params; 


    plan.lambda.Value		= lam; %10; 1; 0.1; 0.01; 0.001; 
    plan.weight         	= 100;
    plan.sigma.Value		= 100;
    plan.rho.Value			= 100;
    plan.sigma.AutoScaling 	= 0;
    plan.rho.AutoScaling 	= 0;

    %% Solver initialization
    plan.Verbose = 1;
    plan.MaxIter = 100;
    plan.AbsStopTol = 1e-6;
    plan.RelStopTol = 1e-6;

    %% Initialize the dictionary
    D0 = zeros(plan.dictSize, 'single'); % Turn on if using single precision
    D0 = rand(plan.dictSize);
    % size(S0)
    plan.dataSize
    % S0 = reshape(S0, plan.dataSize);

    %% Run the CSC algorithm
    isTrainingDictionary=1;

    imageFiles = dir(data_dir);
    imageNames = {imageFiles(arrayfun(@(x) ~x.isdir, imageFiles)).name};
    S0 = [];
    for k=1:1
        S = imreadtif([data_dir imageNames{k} ]);
        % S = resize(S, xsz);
        S = scale1(S);
        S0 = cat(4, S, S0);
    end
    
    if( isempty(save_dir))
        save_dir = 'maps/'
    end
    if(exist(save_dir, 'dir'))
        rmdir(save_dir, 's');
    end
    mkdir(save_dir);
    mkdir([save_dir 'first_learn/']);
    mkdir([save_dir 'second_learn/']);
    
    [resD] = ecsc_cpu(D0, S0, plan, isTrainingDictionary,save_dir,original);

    %%
    close all;
    %plan.lambda.Value	= 0.01;
    S_deploy = S0(:,:,:,1);
    [resX] = ecsc_cpu(resD.G, S_deploy, plan, 0,save_dir,original);

    %% Save the maps
    [s, d, y, gy] = saveMaps(S_deploy, resX.G, resX.Y, plan, 'image', save_dir);
return