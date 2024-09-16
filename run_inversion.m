% clear workspace and command window
clear
clc

% point to current directory for saving files
cur_dir = cd;

% add path to topotoolbox;
addpath(genpath('C:\Users\sfgallen\OneDrive - Colostate\Documents\GitHub\topotoolbox'));

% pick a file identifier for output file naming (called fileTag here)
fileTag = 'test_model';

% define timestep (in years) of forward model used in Bayesian inversion
dt = 10000;

% define burnin and number of model interations (total interation is the
% sum of both values)
burn_in = 5e2;
n_runs = 5e3;
niter = burn_in+n_runs;

%% load prepped stream data and unpack
% data was prepared using TopoToolbox (TTB). All data is topologically
% sorted
stream_data = load('stream_data.mat');
stream_data = stream_data.stream_data;
S = stream_data.S;              % TTB STREAMobj
Sz = stream_data.Sz;            % elevations of river network from DEM in meters
pSz = stream_data.pSz;          % pre-incision network elevations 
S_DA = stream_data.S_DA ;       % upstream drainage area in m^2
Serr = 10.*ones(size(Sz));       % error on stream network elevations in meters

% make sure pSz doesn't flow backward to start. should result in minor
% changes from the interp. surface, but shouldn't be too different
pSz = check_z(S,pSz);

%% set up the parameters for the forward model
% here I assume that the uplift rate is known (I set it here i 0.01 mm/yr)
% as is the timing that incision started (e.g. the age of the youngest
% volcanic flow, here I just assume it is 1 Ma).

% known parameters
U = 0.01e-3;             % uplift rate in m/yr
run_time = 1e6;         % age of pre-incision surface in years

% this leaves unknowns of K, m and n in the stream power model. Note, if
% m/n is assumed to be known, you can reduce the number of paramaters to 2.
% we define the range of priors here
Kexp = [-50,-4];                    % log10(K) -- erodibility.
n = [0.5, 10];                      % slope exponent on stream power model
m = [0.1, 10];                      % drainage area exponent on stream power model

% place priors in a matrix for to make things easier
prior_bounds = [Kexp; n; m];
pranges = [prior_bounds(1,2)-prior_bounds(1,1), prior_bounds(2,2)-prior_bounds(2,1), prior_bounds(3,2)-prior_bounds(3,1)];


%% set up matrices to catch relevant data before entering time MCMC loop
% initialize matrix to catch parameters accepted during the MCMC
params = nan(burn_in+n_runs,3);

% allocate memory to catch likelihoods, acceptance rate 
like_data = nan(burn_in+n_runs,3);

% make sure new matlab session use different random numbers
rng shuffle

%% set up the initial model
% generate initial model parameters 
params(1,1) = -15.5+1e-1*randn(1);           % Kexp
params(1,2) =  4.6+1e-1*randn(1);           % n
params(1,3) =  2.54+1e-1*randn(1);         % m

K_t = (10^(params(1,1))); % steam erodibility
n_t = params(1,2);      % slope exponent
m_t = params(1,3);      % DA exponent

% stream incision model
[Z0] = river_incision_forward_model(S,S_DA,U,pSz,K_t,m_t,n_t,run_time,dt);
% plot(S.distance,Sz,'k.'); hold on
% plot(S.distance,pSz,'b.'); hold on
% plot(S.distance,Z0,'g.'); hold on

%% calculate relevant log-likehoods for the model parameters and results
% log-likelihood of model parameter given the priors
lpcurrent = logprior_uniform(params(1,:),prior_bounds);

% calculate the log-likelihood of the model given the data
llcurrent = mod_loglikelihood(Sz,Z0,Serr);

%% cast relevant variables before entering the loop
current = params(1,:);
lMAP = -Inf;
mMAP = current;
nacc = 0;

%% Use Richards auto-stepsize algorithum from David Egholm 
% initiate variables
acount = 0;            % accepted samples after burn-in
bcount = 0;            % accepted samples before burn-in
rcount = 0;            % rejected samples
accrat = 0;            % acceptance ratio

accfac = 1e-3;    % ???????????? accceptance factor?
k = 0.001;  % stepsize 
status = nan(1,niter);

%% time to run the MCMC loop
h = waitbar(0,'running MCMC...');

for i = 2:niter
    
    % acceptance ratio
    if (i > burn_in)
        accrat = (sum(abs(status((i-100):(i-1))))+1)/100;
    elseif (bcount < burn_in)
        accrat = 0.1;
    else
        accrat = 0.3;
    end
    
    % burnin
    if (bcount < 0.5*burn_in)
        if ((i > 100)&&(bcount < 2))
            k = 0.5;         % 0.1
        elseif ((i > 10)&&(bcount < 2))
            k = 0.05*i;  % 0.001
        end
        
    elseif (bcount < burn_in)
        k = k*((1-accfac) + accfac*accrat/0.2);
    elseif (acount < niter)
        k = k*((1-accfac) + accfac*accrat/0.3);
    end
    
    if (bcount < burn_in)
        Temp = 1.0 + 20.0*(burn_in-bcount)/burn_in;
    else
        Temp = 1.0;
    end
    
    % generate new parameters from the old ones
    step = 0.5*randn(1,length(current)).*pranges;
    candidate = current+k*step;
    
    % run forward model with candidate parameters
    K_t = (10^(candidate(1))); % steam erodibility
    n_t = candidate(2);      % slope exponent
    m_t = candidate(3);      % DA exponent
    
    % stream incision model
    [Z0] = river_incision_forward_model(S,S_DA,U,pSz,K_t,m_t,n_t,run_time,dt);
    
    % calculate model, priors and step for log of acceptance ratio
    lpcandidate = logprior_uniform(candidate,prior_bounds);
    llcandidate = mod_loglikelihood(Sz,Z0,Serr);
    
    % transition probabilities
    lr1 = logproposal(candidate,current,k);
    lr2 = logproposal(current,candidate,k);
    
    % add all the probabilities together for the acceptance ratio
    logalpha = lpcandidate + llcandidate + lr1 - lpcurrent - llcurrent - lr2;
    
     % Take the minimum of the log(alpha) and 0.
    if (logalpha > 0)
        logalpha = 0;
    end
    
    % Generate a U(0,1) random number and take its logarithm.
    logt = log(rand());
    
    % Accept or reject the step.
    if (logt < logalpha)
        
        % Accept the step.
        current = candidate;
        
        %accepted after burnin
        if (bcount > burn_in)
            status(i) = 1;
            acount = acount + 1;
            %if accepted during burnin
        else
            status(i) = -1;
            bcount = bcount + 1;
        end
        
        if i > burn_in
            nacc = nacc+1;
        end
        
        % Update the MAP solution if this one is better.
        if ((lpcandidate + llcandidate) > lMAP)
            lMAP = lpcandidate + llcandidate;
            mMAP = candidate;
        end
        
        % update current likelihoods if needed
        lpcurrent = lpcandidate;
        llcurrent = llcandidate;
      
    else
        % reject the step
        status(i) = 0;
        rcount = rcount + 1;
    end
    
    % store the chain paths and other data
    params(i, :) = current;
    like_data(i,1) = llcurrent;
    like_data(i,2) = lpcurrent; 


    % calculate acceptance ratio if past burnin
    if i > burn_in
        accrate = nacc / (n_runs);
        like_data(i,3) = nacc/(i)*100;
    end

    waitbar(i/(niter),h)
end
close(h)

% save the data
save([cur_dir,'/params_',fileTag,'.mat'],'params');
save([cur_dir,'/mMAP_',fileTag,'.mat'],'mMAP');
save([cur_dir,'/like_dat_',fileTag,'.mat'],'like_data');