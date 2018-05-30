%% 1. Creation of the Synthetic Data
% This section creates the data. |gen| struct store the parameter and
% |data_generation.m| compute everything

clc; clear all; addpath('../functions','R2'); dbstop if error 

% Grid size
gen.xmax = 100; %total length in unit [m]
gen.ymax = 15; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.nx = 500;
gen.ny = 25;

% Generation Method: All method generate with FFTMA a gaussian field.
% 'Normal'              with normal distribution \sim N(gen.mu,gen.std)
% 'LogNormal'   
% 'fromRho':            log transform it with the parameter defined below 
% 'fromK':              generate with FFTMA a field of Hyraulic conductivity and log transform it with the parameter defined below 
gen.method              = 'fromLogPhi';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 3;          % number of well or number of point
gen.covar(1).model      = 'exponential';
gen.covar(1).range0     = [4 40];
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
gen.mu                  = -1.579;%0.27; % parameter of the first field. 
gen.std                 = sqrt(0.22361);%.05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.f ={};
gen.Rho.f.res_matrix    = 0;
gen.Rho.elec.spacing    = 2; % in unit [m] | adapted to fit the fine grid
gen.Rho.elec.bufzone    = 2; % number of electrod to skip 
gen.Rho.elec.config_max = 6000; % number of configuration of electrode maximal 
gen.Rho.i.grid.nx       = gen.nx; % | adapted to fit the fine grid
gen.Rho.i.grid.ny       = gen.ny; % log-spaced grid  | adapted to fit the fine grid
gen.Rho.i.res_matrix    = 0; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.i.max_iterations= 10;

% Other parameter
gen.plotit              = 0;%true;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = '600x40';
gen.seed                = 8;

% Run the function. 
% fieldname = data_generation(gen);
[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);
% fieldname = 'GEN-600x40_2017-12-21_15-44';

%%

filename='GEN-600x40_2018-05-26_16-29';
load(['data_gen/' filename])

f=gen.Rho.f;
i=gen.Rho.i;

% We get 2% error between the calculated and observed of the inversion
misfit = sqrt( mean (((i.output.pseudoCalc - i.output.pseudoObs) ./ (i.output.pseudoObs)).^2 ))

% We get 2.1% error between the observed of the inversion and calculated of the forward
misfit = sqrt( mean (((f.output.pseudo - i.output.pseudoObs) ./ (i.output.pseudoObs)).^2 ))

% We get 1% error between the calculated of the inversion and observed of the inversion
misfit = sqrt( mean (((f.output.pseudo - i.output.pseudoCalc) ./ (i.output.pseudoCalc)).^2 ))


f={};
f.res_matrix        = gen.Rho.f.res_matrix;
f.grid              = gen.Rho.f.grid;    
f.header            = 'Forward';  % title of up to 80 characters
f.job_type          = 0;
f.filepath          = ['data_gen/IO-file/'];
f.readonly          = 0;
f.alpha_aniso       = gen.Rho.f.alpha_aniso;
f.elec_spacing      = gen.Rho.f.elec_spacing;
f.elec_id           = gen.Rho.f.elec_id;
f.rho               = flipud(i.output.res);
f.num_regions       = numel(f.rho);
f.rho_min           = gen.Rho.f.rho_min;
f.rho_avg           = gen.Rho.f.rho_avg;
f.rho_max           = gen.Rho.f.rho_max;

f                   = Matlat2R2(f,gen.Rho.elec); % write file and run forward modeling
misfit = sqrt( mean (((f.output.pseudo - gen.Rho.i.output.pseudoObs) ./ (gen.Rho.i.output.pseudoObs)).^2 ))
misfit = sqrt( mean (((f.output.resistance - gen.Rho.f.output.resistance) ./ (gen.Rho.f.output.resistance)).^2 ))
