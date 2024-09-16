% this is just a basic script to plot the chain paths, likelihood path,
% acceptance rate history, best-fit results, and basin average erosion rate
% history
% 
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date modified: 07/27/2023

% clear workspace and command window
clear
clc

addpath(genpath('C:\Users\sfgallen\Documents\GitHub\topotoolbox'));

% point to current directory for saving files
cur_dir = cd;

% use the same file identifier used for output file naming in the MCMC
% script
fileTag = 'test_model';

% define the timestep in years
dt = 10000;

% declare known variables
U = 0.01e-3;             % uplift rate in m/yr
run_time = 1e6;         % age of pre-incision surface in years

% declare burn in and nruns
burn_in = 5e3;
nruns = 5e4;

%% load the results
% parameter values
params = load(['params_',fileTag,'.mat']);
params = params.params;

% likelihood history
like_data = load(['like_dat_',fileTag,'.mat']);
like_data = like_data.like_data;

% maximum a posteriori (MAP) solution
mMAP = load(['mMAP_',fileTag,'.mat']);
mMAP = mMAP.mMAP;

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

%% run the best-fit forward model
% erodibility 
K_t = 10.^mMAP(1);

% slope and drainage area exponents
n_t = mMAP(2);
m_t = mMAP(3);

% stream incision model
[Z0,t_v,ero] = river_incision_forward_model_w_ero(S,S_DA,U,pSz,K_t,m_t,n_t,run_time,dt);
%% plot the results
% best-fit forward model
figure(1)
plot(S.distance./1e3,pSz,'k.','color',[0.4 0.9 0.4],'markersize',10); hold on
plot(S.distance./1e3,Sz,'k.','color',[0.4 0.4 0.9],'markersize',10); hold on
plot(S.distance./1e3,Z0,'k.','color',[0.2 0.2 0.2],'markersize',10); hold on
xlabel('Distance (km)'); ylabel('Elevation (m)');
legend('pre-incision','modern','best-fit model','Location','northwest')

%
figure(2)
[~,ax] = plotmatrix(params(burn_in+1:end,:));
ylabel(ax(1,1),'log_1_0(K)')
ylabel(ax(2,1),'n')
ylabel(ax(3,1),'m')
xlabel(ax(3,1),'log_1_0(K)')
xlabel(ax(3,2),'n')
xlabel(ax(3,3),'m')

% chain paths
figure(3)
max_it = length(params(:,1));
y_labs = {'log10(K)','n','m'};
for jj = 1:3
    subplot(3,1,jj)
    plot(params(:,jj),'k-'); hold on
    ylabel(y_labs{jj})
    xlim([0 max_it])
end
xlabel('Iteration')

% likelihood and acceptance rate history
figure(4)
subplot(2,1,1)
plot(like_data(:,1),'k-');
xlim([0 max_it])
ylabel('log-likelihood');
subplot(2,1,2)
plot(like_data(:,3),'k-');
xlim([0 max_it])
ylabel('acceptance rate (%) after burn-in');
xlabel('Iteration');

% plot basin average erosion history
figure(5)
plot(t_v./1e6,ero.*1000,'b-','linewidth',2)
xlabel('Model time (Myr)');
ylabel('Basin average erosion rate (mm/yr)')