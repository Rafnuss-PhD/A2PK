%% DATA CREATION
clc; clear all; addpath(genpath('./.')); dbstop if error 
% This section gather all possible way to create the data. |gen| struct
% store the parameter and |data_generation.m| compute everything

% Grid size
gen.xmax = 240; %total length in unit [m]
gen.ymax = 20; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.sx = 10;
gen.sy = 7;

% Generation Method: All method generate with FFTMA a gaussian field.
% 'Normal'              with normal distribution \sim N(gen.mu,gen.std)
% 'LogNormal'   
% 'fromRho':            log transform it with the parameter defined below 
% 'fromK':              generate with FFTMA a field of Hyraulic conductivity and log transform it with the parameter defined below 
gen.method              = 'fromPhi';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 4;          % number of well or number of point
gen.covar(1).model      = 'exponential';
gen.covar(1).range0     = [2.7 27];
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
gen.mu                  = 0.27; % parameter of the first field. 
gen.std                 = .05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 240;
gen.Rho.grid.ny           = 20; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 3; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = true;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'test';
gen.seed                = 23456;

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);

%% How to compute the Resolution Matrix
clear all; addpath(genpath('./.')); dbstop if error 
load('result-A2PK/GEN-test_2017-09-26_17-22.mat'); %alpha = 27339.756;

% Show raw result of inversion
[X2,Y2] = meshgrid(Sigma.x_raw,Sigma.y_raw);
figure(1); clf; c_axis=[ min(sigma_true(:)) max(sigma_true(:))];
subplot(2,1,1); surface(grid_gen.X,grid_gen.Y,sigma_true,'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis(c_axis); 
subplot(2,1,2); surface(X2,Y2,Sigma.d_raw,'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis(c_axis); 
colorbar('south')

% COmpute coordinate of the raw data
y = zeros(1,numel(gen.Rho.i.yy)-1);
for i=2:numel(gen.Rho.i.yy)-1
    y(i) = 2*gen.Rho.i.yy(i) - y(i-1);
end
x = zeros(1,numel(gen.Rho.i.xx)-1);
[~,id] = min(gen.Rho.i.xx.^2);
for i=(id+1):numel(gen.Rho.i.xx)-1
    x(i) = gen.Rho.i.xx(i) + gen.Rho.i.xx(i) - x(i-1) ;
end
x(1:id-1) = sort(max(grid_gen.x)-x((end-id+2):end));
[X,Y] = meshgrid(x,y);

% Compute the sensitivity matrix
sens = gen.Rho.i.output.J*(gen.Rho.i.output.Wd*gen.Rho.i.output.Wd')*gen.Rho.i.output.J';

figure(2); clf; c_axis=[ min(gen.Rho.i.output.sen(:)) max(gen.Rho.i.output.sen(:))];
subplot(2,1,1); surface(X,Y,reshape(diag(sens),numel(gen.Rho.i.yy)-1,numel(gen.Rho.i.xx)-1),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; xlim([0 240]); ylim([0 20]); %caxis(c_axis);
subplot(2,1,2); surface(X2,Y2,flipud(gen.Rho.i.output.sen),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis(c_axis); 


% Compute the Resolution Matrix
Res = (sens + alpha * gen.Rho.i.output.R)^(-1) * sens;

% Plot few stuff
x=120; y=15; [~,i] = min( sqrt( (x-X(:)).^2 + (y-Y(:)).^2));
figure(3); clf; 
subplot(2,1,1); surface(X,Y,reshape(Res(i,:),numel(gen.Rho.i.yy)-1,numel(gen.Rho.i.xx)-1),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on
hold on; plot(X(i),Y(i),'.r','MarkerSize',69);  colorbar('south'); xlim([0 240]); ylim([0 20])
subplot(2,1,2); surface(X,Y,reshape(Res(i,:),numel(gen.Rho.i.yy)-1,numel(gen.Rho.i.xx)-1),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on


figure(4); clf; 
subplot(2,1,1); surface(X,Y,reshape(diag(Res),numel(gen.Rho.i.yy)-1,numel(gen.Rho.i.xx)-1),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; colorbar('south'); xlim([0 240]); ylim([0 20])
subplot(2,1,2); surface(X,Y,reshape(diag(Res),numel(gen.Rho.i.yy)-1,numel(gen.Rho.i.xx)-1),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on;

figure(5); clf; 
subplot(2,1,1); surface(X,Y,reshape(sum(Res,2),numel(gen.Rho.i.yy)-1,numel(gen.Rho.i.xx)-1),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis([0 1]); colorbar('south'); xlim([0 240]); ylim([0 20])
subplot(2,1,2); surface(X,Y,reshape(sum(Res,2),numel(gen.Rho.i.yy)-1,numel(gen.Rho.i.xx)-1),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis([0 1]);



%%
clear all; addpath(genpath('./.')); dbstop if error 
load('result-A2PK/GEN-test_2017-09-26_17-22.mat'); %alpha = 27339.756;

[kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
parm.nscore=1;
Nscore = nscore(kern, parm, 0);
Sec.x=Sigma.x_raw; Sec.y=Sigma.y_raw; [Sec.X , Sec.Y] = meshgrid(Sigma.x_raw, Sigma.y_raw); 
Sec.d = reshape(Nscore.forward(Sigma.d_raw(:)) ,numel(Sec.y),numel(Sec.x));
Prim.d = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);
Prim.x = grid_gen.x; Prim.y = grid_gen.y; Prim.X = grid_gen.X; Prim.Y = grid_gen.Y;

c_axis=[ min(Prim.d(:)) max(Prim.d(:)) ];
subplot(2,1,1);imagesc(grid_gen.x, grid_gen.y, Prim.d); caxis(c_axis); title('zt True field');
subplot(2,1,2); imagesc(Sec.x, Sec.y, Sec.d); caxis(c_axis); title('Inverted field'); axis tight; box on


% Built the matrix G which link the true variable Prim.d to the measured coarse scale d
G = zeros(numel(Sec.d), numel(Prim.d));
for i=1:numel(Sec.d)
    f = griddedInterpolant({Sec.y,Sec.x},reshape(Sigma.Res(i,:),numel(Sec.y),numel(Sec.x)),'Linear');
    res_t = f({Prim.y,Prim.x});
    G(i,:) = res_t(:) ./sum(res_t(:));
end

i=1000;
figure(1); clf; 
subplot(2,1,1);imagesc(Sec.x,Sec.y,reshape(Sigma.Res(i,:),numel(Sec.y),numel(Sec.x))); 
subplot(2,1,2);imagesc(Prim.x,Prim.y,reshape(G(i,:), numel(Prim.y), numel(Prim.x))); 


% Compute coarse scale d
Test_Sec_d = reshape(G * Prim.d(:), numel(Sec.y), numel(Sec.x));
figure(2); clf; 
subplot(2,1,1);imagesc(Sec.x,Sec.y,Sec.d); 
subplot(2,1,2);imagesc(Sec.x,Sec.y,Test_Sec_d); 

imagesc()

% Simulate sampling
Prim_pt = sampling_pt(Prim,Prim.d,2,0);

% Direct Simulation
parm.k.covar = gen.covar;
parm.seed_U = 'shuffle';
parm.n_real=1;

[Res] = ESA2P(Prim.y,Prim.x,Prim_pt,Sec.d,G,parm);


%%

% Define the true case
[kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
parm.nscore=1;
Nscore = nscore(kern, parm, 0);
Sec.x=Sigma.x_raw; Sec.y=Sigma.y_raw; 
Sec.d = reshape(Nscore.forward(Sigma.d_raw(:)) ,numel(Sec.y),numel(Sec.x));
Prim.d = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);
c_axis=[ min(Prim.d(:)) max(Prim.d(:)) ];
subplot(2,1,1);imagesc(grid_gen.x, grid_gen.y, Prim.d); caxis(c_axis); title('zt True field');
subplot(2,1,2); imagesc(Sec.x, Sec.y, Sec.d); caxis(c_axis); title('Inverted field'); axis tight; box on

% Built the matrix G which link the true variable Prim.d to the measured coarse scale d
G = zeros(numel(Sec.X),Prim.nxy);
for ij=1:Prim.nxy
     [~,id_z]=min((Prim.X(ij)-Sec.X(:)).^2 + (Prim.Y(ij)-Sec.Y(:)).^2);
     G(id_z,ij)=1;%/(NSigma.dx*NSigma.dy);
end
for ij = 1:numel(Sec.X)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end

% Compute coarse scale d
Sec.d = reshape(G * Prim.d(:), numel(Sec.y), numel(Sec.x));

% Simulate sampling
Prim_pt = sampling_pt(Prim,Prim.d,2,0);


figure(1);clf
c_axis=[ min(Prim.d(:)) max(Prim.d(:))];
%subplot(2,1,1);surf(Prim.x, Prim.y, Prim.d,'EdgeColor','none'); caxis(c_axis); title('Prim.d True field');;view(2); axis tight; set(gca,'Ydir','reverse'); box on
%subplot(2,1,2);surf(Sec.x, Sec.y, Sec.d,'EdgeColor','none'); caxis(c_axis); title('d Average of True field');view(2); axis tight; set(gca,'Ydir','reverse'); box on
subplot(1,2,1); hold on; imagesc(Prim.x, Prim.y, Prim.d); caxis(c_axis); title('Prim.d True field');
scatter(Prim_pt.x,Prim_pt.y,[],Prim_pt.d,'filled','MarkerEdgeColor','k'); axis tight equal; box on
subplot(1,2,2);imagesc(Sec.x, Sec.y, Sec.d); caxis(c_axis); title('d Average of True field'); axis tight equal; box on


%% Direct Simualtion
parm.k.covar = covar;
parm.seed_U = 'shuffle';
parm.n_real=100;

[Res] = ESA2P(Prim.y,Prim.x,Prim_pt,Sec.d,G,parm);

for i_real =1:parm.n_real
    r=Res(:,:,i_real);
    s(i_real) = var(r(:))
    error_pt(i_real)=sqrt(sum((r(Prim_pt.id(:))-Prim_pt.d).^2));
    Res_d(:,:,i_real) = reshape(G * r(:), numel(Sec.y), numel(Sec.x));
    error_d(i_real) = sqrt(sum(sum((Res_d(:,:,i_real) - Sec.d ).^2)));
end

figure(3); clf
c_axis=[ min(Prim.d(:)) max(Prim.d(:))];
subplot(parm.n_real+1,2,1);imagesc(Prim.x, Prim.y, Prim.d); caxis(c_axis); title('Prim.d True field'); hold on;
scatter(Prim_pt.x,Prim_pt.y,[],Prim_pt.d,'filled','MarkerEdgeColor','k'); axis tight equal; box on
subplot(parm.n_real+1,2,2);imagesc(Sec.x, Sec.y, Sec.d); caxis(c_axis); title('d Average of True field'); axis tight equal; box on
for i_real =1:parm.n_real
    subplot(parm.n_real+1,2,i_real*2+1);imagesc(Prim.x, Prim.y,Res(:,:,i_real)); caxis(c_axis); title('Res');view(2); axis tight equal;box on
    hold on; scatter(Prim_pt.x,Prim_pt.y,[],Prim_pt.d,'filled','MarkerEdgeColor','k'); axis tight equal; box on
    subplot(parm.n_real+1,2,i_real*2+2);imagesc(Sec.x, Sec.y, Res_d(:,:,i_real)); caxis(c_axis); title('Res d Average'); axis tight equal; box on
end




%% Sequential Simulation

parm.k.covar = covar;
parm.seed_path = 'shuffle';
parm.seed_search = 'shuffle';
parm.seed_U = 'shuffle';

parm.k.method = 'sbss';
parm.k.lookup = false;
parm.k.nb = 30;
parm.k.wradius = 30;

parm.mg = 0;

parm.n_real=10;

Prim_sim=Prim;
Prim_sim.d=[];

Res = SSA2P(Prim.y,Prim.x,Prim_pt,Sec.d,G,parm);

for i_real =1:parm.n_real
    r=Res(:,:,i_real);
    s(i_real) = var(r(:));
   % Res(:,:,i_real) = Res(:,:,i_real) ./ var(r(:));
    Res_d(:,:,i_real) = reshape(G * r(:), numel(Sec.y), numel(Sec.x));
end


% Plot
figure(2);clf
c_axis=[ min(Prim.d(:)) max(Prim.d(:))];
subplot(parm.n_real+1,2,1);imagesc(Prim.x, Prim.y, Prim.d); caxis(c_axis); title('Prim.d True field');;view(2); axis tight; set(gca,'Ydir','reverse'); box on
hold on; scatter(Prim_pt.x,Prim_pt.y,[],Prim_pt.d,'filled','MarkerEdgeColor','k'); axis tight equal; box on
subplot(parm.n_real+1,2,2);imagesc(Sec.x, Sec.y, Sec.d); caxis(c_axis); title('d Average of True field');view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
for i_real =1:parm.n_real
    subplot(parm.n_real+1,2,i_real*2+1);imagesc(Prim.x, Prim.y,Res(:,:,i_real)); caxis(c_axis); title('Res');view(2); axis tight; set(gca,'Ydir','reverse'); box on
 hold on; scatter(Prim_pt.x,Prim_pt.y,[],Prim_pt.d,'filled','MarkerEdgeColor','k'); axis tight equal; box on
    subplot(parm.n_real+1,2,i_real*2+2);imagesc(Sec.x, Sec.y, Res_d(:,:,i_real)); caxis(c_axis); title('Res d Average');view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
end



figure; clf;hold on
a=variogram_gridded_perso(Prim.d);
plot(a)
for i_real=1:parm.n_real
a=variogram_gridded_perso(Res(:,:,i_real));
plot(a)
end
