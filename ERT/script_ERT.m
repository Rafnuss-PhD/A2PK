%% Script to use A2PK with ERT
% This script will generate the dataset and the result for the paper ... 


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
gen.ny = 30;

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

% Electrical inversion
gen.f ={};
gen.f.res_matrix    = 0;
gen.f.n_plus        = 10;    
gen.elec.spacing    = 2; % in unit [m] | adapted to fit the fine grid
gen.elec.bufzone    = 2; % number of electrod to skip 
gen.elec.config_max = 6000; % number of configuration of electrode maximal 
gen.i.grid.nx       = 100; % | adapted to fit the fine grid
gen.i.grid.ny       = 15; % log-spaced grid  | adapted to fit the fine grid
gen.i.n_plus        = 10;    
gen.i.res_matrix    = 3; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.i.max_iterations= 10;

% Other parameter
gen.plotit              = 1;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'AsSumbitted';
gen.seed                = 4;

% Run the function. 
data_generation(gen);
% [fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);



%% 2. Area-to-point Kriging
% This section load the synthetic data create the covariance matrices and
% compute the kriging weight to finally produce the Kriging map
clc; clear all;
addpath('../functions','R2');
fieldname = 'GEN-AsSumbitted_2018-07-09_16-34';
load(['result/' fieldname]);

% Add some nice colormap
addpath('C:\Users\rnussba1\Documents\MATLAB\Colormaps\')
colormap(viridis())

% Normal Score based on empirical CDF
% [kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
% parm.nscore=1;
% Nscore = nscore(kern, parm, 0);
% Sec=Sigma; 
% Sec.d = reshape(Nscore.forward(Sec.d(:)) ,numel(Sec.y),numel(Sec.x));
% Prim.d = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);
% Prim.x = grid_gen.x; Prim.y = grid_gen.y; Prim.X = grid_gen.X; Prim.Y = grid_gen.Y;

% Normal Score based on known distribution of Prim and Sec
Nscore.forward = @(x) (log(x./43)/1.4-gen.mu)./gen.std;
Nscore.inverse = @(x) 43.*exp(1.4*(x.*gen.std+gen.mu));
Sec=Sigma; 
Sec.d = Nscore.forward(Sigma.d);
Prim.d_full = Nscore.forward(1000./gen.f.rho);
Prim.d = reshape(Prim.d_full(gen.f.grid.inside),grid_gen.ny,grid_gen.nx);
Prim.x = grid_gen.x; Prim.y = grid_gen.y; 
[Prim.X, Prim.Y] = meshgrid(Prim.x, Prim.y);

xy_axis=[min(grid_gen.x) max(grid_gen.x) min(grid_gen.y) max(grid_gen.y)];

% Figure of intial data
figure(1); clf;c_axis=[ min(Prim.d(:)) max(Prim.d(:)) ];
subplot(2,1,1);surface(Prim.x, Prim.y, Prim.d,'EdgeColor','none','facecolor','flat'); axis(xy_axis); caxis(c_axis); title('zt True field'); set(gca,'Ydir','reverse');
subplot(2,1,2); surface(Sec.x, Sec.y, Sec.d,'EdgeColor','none','facecolor','flat'); set(gca,'Ydir','reverse'); caxis(c_axis); title('Inverted field'); axis(xy_axis); box on;colormap(viridis())
% export_fig -eps 'Prim_and_sec'

figure(2); clf; c_axis_n=log10([ min(sigma_true(:)) max(sigma_true(:)) ]);
subplot(2,1,1);surface(grid_gen.x, grid_gen.y, log10(sigma_true),'EdgeColor','none','facecolor','flat'); 
caxis(c_axis_n); title('True Electrical Conductivity \sigma^{true}');axis equal tight; box on; xlabel('x');ylabel('y'); set(gca,'Ydir','reverse');
subplot(2,1,2); surface(Sigma.x, Sigma.y, log10(Sigma.d),'EdgeColor','none','facecolor','flat'); 
caxis(c_axis_n); title('Inverted Electrical Conductivity U \sigma^{est}'); axis(xy_axis); box on; xlabel('x');ylabel('y');set(gca,'Ydir','reverse'); colorbar('southoutside')
colormap(viridis())
% export_fig -eps 'ElecCondRawData'

% Built the matrix G which link the true variable Prim.d to the measured coarse scale d
G = zeros(numel(Sec.d), numel(Prim.d));
Sec.dnorm = zeros(size(Sec.d));
for i=1:numel(Sec.d)
    Res = reshape(Sigma.res(i,:),numel(Sec.y),numel(Sec.x));
    %Res = Res./sum(Res(:));
    % IN
    Res_in = reshape(Res(gen.i.grid.inside), gen.i.grid.ny,gen.i.grid.nx);
    F = griddedInterpolant({Sec.y,Sec.x},Res,'linear');
    res_t = F({Prim.y,Prim.x});
    G(i,:) = res_t(:) ./sum(res_t(:)) .* sum(Res_in(:));
    % OUT
    Sec.dnorm(i) = Res(~gen.i.grid.inside)'*Sec.d(~gen.i.grid.inside);
end

% U = zeros(numel(Sec.d), numel(Prim.d));
% for i=1:numel(Sec.d)
%     Res = zeros(numel(Sec.y),numel(Sec.x));
%     Res(i)=1;
%     F = griddedInterpolant({Sec.y,Sec.x},Res,'linear');
%     res_t = F({Prim.y,Prim.x});
%     U(i,:) = res_t(:) ./sum(res_t(:));
% end
% U(isnan(U))=0;
% 
% 
% G2 =Sigma.res*U;
% imagesc(reshape(G * Prim.d(:), numel(Sec.y), numel(Sec.x)))
% figure; imagesc(reshape(G2 * Prim.d(:), numel(Sec.y), numel(Sec.x)))
% figure; imagesc(reshape((G2*Prim.d(:)-G * Prim.d(:))./(G * Prim.d(:)), numel(Sec.y), numel(Sec.x)))

% View Resolution matrix
% OLD large grid: i=[1838 4525 8502]; th=.05; x_lim=[16 25; 35 65; 90 99]; y_lim=[-.12 7; -.12 19.57; -.12 8]; ccaxis=[-.02 .02; -.002 .002; -.01 .01];
% OLD grid: i=[971 3138 5534]; th=.05; x_lim=[8 20; 35 65; 80 100]; y_lim=[-.007 7; -.007 14.37; -.007 10]; ccaxis=[-.05 .05; -.01 .01; -.02 .02];
i=[726 1789 3162]; th=.05; x_lim=[7.5 21.5; 35.5 63.5; 80.5 99.5]; y_lim=[0.25 8.1447; 0.25 14.5; 0.25 10.7105]; ccaxis=[-.05 .05; -.01 .01; -.02 .02];
figure(3); clf; colormap(viridis())
subplot(3,numel(i),[1 numel(i)]);hold on; surface(Sec.X,Sec.Y,reshape(diag(Sigma.res),numel(Sec.y),numel(Sec.x)),'EdgeColor','none','facecolor','flat');   axis(xy_axis); 
scatter3(Sec.X(i),Sec.Y(i),100*ones(numel(i),1),'sr','filled')
for ii=1:numel(i)
    rectangle('Position',[x_lim(ii,1) y_lim(ii,1) range(x_lim(ii,:)) range(y_lim(ii,:))],'EdgeColor','r','linewidth',1)
end
view(2);  axis(xy_axis); set(gca,'Ydir','reverse'); box on; xlabel('x');ylabel('y'); title('Diagonal of the resolution matrix R'); colorbar('southoutside')
for ii=1:numel(i)
    subplot(3,numel(i),numel(i)+ii);hold on; surface(Sec.X,Sec.Y,reshape(Sigma.res(i(ii),:),numel(Sec.y),numel(Sec.x)),'EdgeColor','none','facecolor','flat');  
    scatter3(Sec.X(i(ii)),Sec.Y(i(ii)),100,'sr','filled')
    view(2); axis tight equal; set(gca,'Ydir','reverse'); box on; axis equal tight; xlabel('x');ylabel('y'); title('Resolution of the red dot R(i,:)'); xlim(x_lim(ii,:)); ylim(y_lim(ii,:));
    caxis(ccaxis(ii,:))
    colorbar('southoutside')
    
    subplot(3,numel(i),2*numel(i)+ii);hold on; surface(Prim.X,Prim.Y,reshape(G(i(ii),:),numel(Prim.y),numel(Prim.x)),'EdgeColor','none','facecolor','flat');  
    scatter3(Sec.X(i(ii)),Sec.Y(i(ii)),100,'sr','filled')
    view(2); axis tight equal; set(gca,'Ydir','reverse'); box on; axis equal tight; xlabel('x');ylabel('y'); title('Resolution of the red dot R(i,:)'); xlim(x_lim(ii,:)); ylim(y_lim(ii,:));
    caxis(ccaxis(ii,:))
    colorbar('southoutside')
    
%     subplot(2,numel(i),2*numel(i)+ii);hold on; surface(Prim.X,Prim.Y,reshape(G(i(ii),:), numel(Prim.y), numel(Prim.x)),'EdgeColor','none','facecolor','flat'); scatter3(Sec.X(i(ii)),Sec.Y(i(ii)),100,'sr','filled')
%     view(2); axis tight equal; set(gca,'Ydir','reverse'); box on;  axis equal;  xlabel('x');ylabel('y'); xlim(x_lim(ii,:)); ylim(y_lim(ii,:))
end
% export_fig -eps 'Res'


% Compute coarse G-transformed of the true fine scale and compare it to the
% inverted field
Test_Sec_d = reshape(G * Prim.d(:), numel(Sec.y), numel(Sec.x))+Sec.dnorm;
figure(4); clf; colormap(viridis())
subplot(3,1,1);surface(Sec.X,Sec.Y,Sec.d,'EdgeColor','none','facecolor','flat'); view(2); set(gca,'Ydir','reverse'); caxis([-3 3]); axis(xy_axis); box on; xlabel('x');ylabel('y'); title('Inverted Electrical Conductivity \Sigma^{est}')
subplot(3,1,2);surface(Sec.X,Sec.Y,Test_Sec_d,'EdgeColor','none','facecolor','flat'); view(2); set(gca,'Ydir','reverse');  caxis([-3 3]);axis(xy_axis); box on; xlabel('x');ylabel('y');% colorbar('southoutside'); title('G-transform of the True Electrical Conductivity Gz^{true}')
subplot(3,1,3);surface(Sec.X,Sec.Y,Test_Sec_d-Sec.d,'EdgeColor','none','facecolor','flat'); view(2); set(gca,'Ydir','reverse'); axis(xy_axis); box on; xlabel('x');ylabel('y'); colorbar('southoutside'); title('Error'); caxis('auto')
% export_fig -eps 'SigmavsGztrze'


% Generate a sampling
Prim_pt = sampling_pt(Prim,Prim.d,1,1);
% Prim_pt = sampling_pt(Prim,Prim.d,1,0);

% Compute the covariance of the data error
CSigma= (gen.std*1.4)^(-2)* gen.i.output.Cov;

% Compute the covariance of the spatial model
covar = kriginginitiaite(gen.covar);
Cz = covar.g(squareform(pdist([Prim.X(:) Prim.Y(:)]*covar.cx)));
Cz=sparse(Cz);
CzZ = Cz * G';
CZ = G * CzZ;

% Combine both covariance an built the Kriging System
CZest = CZ+CSigma;
Czh = Cz(Prim_pt.id,Prim_pt.id);
Czzh = Cz(Prim_pt.id,:);
Czhd = CzZ( Prim_pt.id ,:);
CCa = [ CZest Czhd' ; Czhd Czh ];
CCb = [ CzZ' ; Czzh ];


% Solve the kriging system
Cz0=covar.g(0);
W=zeros(numel(Prim.x)*numel(Prim.y),numel(Sec.d(:))+Prim_pt.n);
S=nan(numel(Prim.y),numel(Prim.x));
CCainv = inv(CCa);
parfor ij=1:numel(Prim.x)*numel(Prim.y)
     W(ij,:) = CCainv * CCb(:,ij);
     S(ij) =  Cz0 - W(ij,:) * CCb(:,ij);
end
% save(['result/' fieldname '_cond'],'W','S','Prim_pt','G','Nscore','Sec','Prim')
% load(['result/' fieldname '_cond'],'W','S','Prim_pt','G','Nscore','Sec','Prim')
% save(['result/' fieldname '_cond_noHD'],'W','S','Prim_pt','G','Nscore','Sec','Prim')
% load(['result/' fieldname '_cond_noHD'],'W','S','Prim_pt','G','Nscore','Sec','Prim')


% Compute the Kriging map
zh = reshape( W * [Sec.d(:)-Sec.dnorm(:) ; Prim_pt.d], numel(Prim.y), numel(Prim.x));
zhtest = reshape( W * [Test_Sec_d(:)-Sec.dnorm(:) ; Prim_pt.d], numel(Prim.y), numel(Prim.x));

figure(5); clf;   colormap(viridis())
subplot(2,1,1);surf(Prim.x,Prim.y,zh,'EdgeColor','none','facecolor','flat'); caxis([-3 3])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); colorbar('southoutside'); title('Kriging Estimate')
subplot(2,1,2);surf(Prim.x,Prim.y,S,'EdgeColor','none','facecolor','flat'); caxis([0 1])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); colorbar('southoutside'); title('Kriging Estimate Error Variance')
%export_fig -eps 'Krig'

figure(55); clf;   colormap(viridis())
subplot(3,1,1);surf(Prim.x,Prim.y,zh,'EdgeColor','none','facecolor','flat'); caxis([-3 3])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); title('Kriging Estimate')
subplot(3,1,2);surf(Prim.x,Prim.y,zhtest,'EdgeColor','none','facecolor','flat');
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y');  title('Kriging Estimate ')
subplot(3,1,3);surf(Prim.x,Prim.y,zh-zhtest,'EdgeColor','none','facecolor','flat');
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); colorbar('southoutside'); title('Error Kriging Estimate')

figure(51); clf;
S(S<=eps)=eps;
histogram( (Prim.d(:)-zh(:))./ sqrt(S(:)) )

figure(31); clf; colormap(viridis())
subplot(2,1,1);hold on; surface(Sec.X,Sec.Y,reshape(diag(CSigma),numel(Sec.y),numel(Sec.x)),'EdgeColor','none','facecolor','flat');  
axis equal tight; box on; xlabel('x');ylabel('y'); set(gca,'Ydir','reverse'); axis(xy_axis);
subplot(2,1,2);hold on; surface(Sec.X,Sec.Y,reshape(diag(CZ),numel(Sec.y),numel(Sec.x)),'EdgeColor','none','facecolor','flat');  
axis equal tight; box on; xlabel('x');ylabel('y'); set(gca,'Ydir','reverse'); axis(xy_axis);
%export_fig -eps 'Cov'

%% Simulation of the Area-to-point Kriging 
rng('shuffle');

parm.n_real = 30;
covar = kriginginitiaite(gen.covar);
zcs=nan(numel(Prim.y), numel(Prim.x),parm.n_real);
z=nan(numel(Sec.y), numel(Sec.x),parm.n_real);
for i_real=1:parm.n_real
    % zs(:,:,i_real) = fftma_perso(covar, struct('x',Prim.x,'y',Prim.y));
    zs = fftma_perso(gen.covar, grid_gen);
    zhs = W * [G * zs(:) ; zs(Prim_pt.id)./1.3];
    r = zh(:) + (zs(:) - zhs(:));
    zcs(:,:,i_real) = reshape( r, numel(Prim.y), numel(Prim.x));
    z(:,:,i_real)=reshape(G*r(:), numel(Sec.y), numel(Sec.x) );
end

% Figure
figure(6);clf; colormap(viridis())
c_axis=[ -3 3];
subplot(4,1,1);surf(Prim.x, Prim.y, Prim.d,'EdgeColor','none','facecolor','flat'); caxis(c_axis);view(2); axis tight equal; set(gca,'Ydir','reverse'); colorbar;
% hold on; scatter(Prim_pt.x,Prim_pt.y,'filled','r')
subplot(4,1,2);surf(Prim.x, Prim.y, zcs(:,:,1),'EdgeColor','none','facecolor','flat'); caxis(c_axis);view(2); axis tight equal; set(gca,'Ydir','reverse'); 
subplot(4,1,3);surf(Prim.x, Prim.y, mean(zcs,3),'EdgeColor','none','facecolor','flat'); caxis(c_axis);view(2); axis tight equal; set(gca,'Ydir','reverse'); 
subplot(4,1,4);surf(Prim.x, Prim.y, std(zcs,[],3),'EdgeColor','none','facecolor','flat'); view(2); axis tight equal; set(gca,'Ydir','reverse'); colorbar;
% export_fig -eps 'PrimOverview'

figure(61);clf; colormap(viridis())
surf(Prim.x, Prim.y,  mean(zcs,3)-Prim.d,'EdgeColor','none','facecolor','flat');view(2); axis tight equal; set(gca,'Ydir','reverse'); colorbar;


figure(62); clf;
S(S<=eps)=eps;
tmp = (zcs-repmat(zh,1,1,parm.n_real)) ./ sqrt(repmat(S,1,1,parm.n_real));
histogram( tmp(:) ,-3:.1:3); xlim([-3 3])

figure(63); clf; hold on;
scatter(repmat(zh(:),parm.n_real,1),zcs(:),'.k')
plot([min(zh(:)) max(zh(:))], [min(zh(:)) max(zh(:))], 'r')
xlabel('True Electrical conductivity')
ylabel('Simulated Electrical conductivity')

figure(7);clf; colormap(viridis())
c_axis=[ -3 3];
subplot(4,1,1);surf(Sec.x, Sec.y, Sec.d,'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis(xy_axis); set(gca,'Ydir','reverse');  colorbar;
subplot(4,1,2);surf(Sec.x, Sec.y, z(:,:,1),'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2);axis(xy_axis); set(gca,'Ydir','reverse'); 
subplot(4,1,3);surf(Sec.x, Sec.y, mean(z,3),'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis(xy_axis); set(gca,'Ydir','reverse');
subplot(4,1,4);surf(Sec.x, Sec.y, std(z,[],3),'EdgeColor','none','facecolor','flat');  title('d Average of True field');view(2); axis(xy_axis); set(gca,'Ydir','reverse'); colorbar;
% export_fig -eps 'SecOverview'


figure(8);clf;colormap(viridis())
subplot(2,1,1);surface(Sec.X,Sec.Y,Test_Sec_d-Sec.d,'EdgeColor','none','facecolor','flat'); view(2); set(gca,'Ydir','reverse');  axis equal tight; box on; xlabel('x');ylabel('y'); caxis([-.6 .6]);colorbar('southoutside');
subplot(2,1,2);surf(Sec.x, Sec.y, mean(z,3)-Sec.d,'EdgeColor','none','facecolor','flat'); view(2); axis tight equal; set(gca,'Ydir','reverse'); colorbar('southoutside');caxis([-.6 .6])
% export_fig -eps 'GztrueGzsim'

% Compute the Variogram and Histogram of realiaztions
% parm.n_real=500;
vario_x=nan(parm.n_real,numel(Prim.x));
vario_y=nan(parm.n_real,numel(Prim.y));
for i_real=1:parm.n_real
    r = zcs(:,:,i_real);
    %r = (r(:)-mean(r(:)))./std(r(:));
    [vario_x(i_real,:),vario_y(i_real,:)]=variogram_gridded_perso(reshape( r(:), numel(Prim.y), numel(Prim.x)));
end
[vario_prim_x,vario_prim_y]=variogram_gridded_perso(Prim.d);

figure(9);clf;
subplot(2,1,1);  hold on; 
h1=plot(Prim.x(1:2:end)-Prim.x(1),vario_x(:,1:2:end)','b','color',[.5 .5 .5]);
h2=plot(Prim.x-Prim.x(1),vario_prim_x,'-r','linewidth',2);
h3=plot(Prim.x-Prim.x(1),1-covar.g((Prim.x-Prim.x(1))*covar.cx(1)),'--k','linewidth',2);
xlim([0 60]); xlabel('Lag-distance h_x ');ylabel('Variogram \gamma(h_x)'); ylim([0 1.5])
legend([h1(1) h2 h3],'500 realizations','True field','Theorical Model')
subplot(2,1,2); hold on; 
h1=plot(Prim.y(1:2:end)-Prim.y(1),vario_y(:,1:2:end)','b','color',[.5 .5 .5]);
h2=plot(Prim.y-Prim.y(1),vario_prim_y,'-r','linewidth',2);
h3=plot(Prim.y-Prim.y(1),1-covar.g((Prim.y-Prim.y(1))*covar.cx(4)),'--k','linewidth',2);
xlim([0 6]); xlabel('Lag-distance h_y ');ylabel('Variogram \gamma(h_y)')
legend([h1(1) h2 h3],'500 realizations','True field','Theorical Model'); ylim([0 1.2])
% export_fig -eps 'Vario'


figure(10);clf; hold on; colormap(viridis())
for i_real=1:parm.n_real
    r=zcs(:,:,i_real);
    %r = (r(:)-mean(r(:)))./std(r(:));
    [f,xi] = ksdensity(r(:));
    h1=plot(xi,f,'b','color',[.5 .5 .5]);
end
[f,xi] = ksdensity(Prim.d(:));
h4=plot(xi,f,'-r','linewidth',2);
[f,xi] = ksdensity(Prim_pt.d(:));
h2=plot(xi,f,'-g','linewidth',2);
h3=plot(xi,normpdf(xi),'--k','linewidth',2);
xlabel('Lag-distance h ');ylabel('Variogram \gamma(h)')
legend([h1 h2 h3 h4],'500 realizations','True field','Theorical Model','Sampled value (well)')
% export_fig -eps 'Histogram'



%% Forward simulation
% Put the realization in the forward ERT

% parm.n_real=12;
fsim_pseudo=nan(numel(gen.f.output.pseudo),parm.n_real);
fsim_resistance=nan(numel(gen.f.output.resistance),parm.n_real);
rho = 1000./Nscore.inverse(zcs);

% fieldname= 'AsSumbitted_2018-06-27_10-07';

f  = gen.f;

f.content = fileread( ['result/' fieldname '_IO-file/forward/R2.in']);
F = griddedInterpolant({Sec.y,Sec.x},1000./Nscore.inverse(Sec.d),'linear');
f.rho = F({f.grid.y, f.grid.x});

% f.rho = i_true_mod.output.res;


for i_real=1:parm.n_real
    f.filepath          = ['data_gen/IO-file-' num2str(i_real) '/'];
    mkdir(f.filepath)
    %copyfile('R2/R2.exe',[f.filepath 'R2.exe'])
    copyfile(['result/' fieldname '_IO-file/forward/electrodes.dat'],[f.filepath 'electrodes.dat'])
    copyfile(['result\' fieldname '_IO-file\forward\protocol.dat'],[f.filepath 'protocol.dat'])
end

parfor i_real=1:parm.n_real
    f_tmp=f;
    f_tmp.filepath          = ['data_gen/IO-file-' num2str(i_real) '/'];
    f_tmp.rho(f_tmp.grid.inside)= rho(:,:,i_real);
    %f_tmp.rho= gen.f.rho;
    %f_tmp.rho(f_tmp.grid.inside)= gen.f.rho(f_tmp.grid.inside);
    [fsim_resistance(:,i_real), fsim_pseudo(:,i_real)] = Matlat2R2min(f_tmp);
end

sqrt(mean(((fsim_pseudo - gen.i.output.pseudoObs)./(gen.i.output.pseudoObs) ).^2))
sqrt(mean(((fsim_resistance - gen.f.output.resistancewitherror)./(gen.f.output.resistancewitherror) ).^2))

scatter(gen.f.pseudo_x,gen.f.pseudo_y,[],mean(fsim_resistance,2) - gen.f.output.resistancewitherror,'filled')

% Compute the misfit
WRMSE=nan(1,parm.n_real);
for i_real=1:parm.n_real  
    WRMSE(i_real) = sqrt(mean(((fsim_pseudo(:,i_real) - gen.i.output.pseudoObs)./(gen.i.b_wgt*gen.i.output.pseudoObs) ).^2));
end
% save(['result/' fieldname '_sim'],'zcs','fsim_pseudo','fsim_resistance','WRMSE','WRMSE_uncond')


figure(11); hold on;
histogram(WRMSE); histogram(WRMSE_uncond); plot([1 1],[0 100]); set(gca, 'XScale', 'log')
xlabel('WRMSE'); ylabel('Histogram'); legend('Conditional realizations', 'Unconditional realizations', 'True initial field')
% export_fig -eps 'misfit-hist'


figure(12);clf; colormap(viridis());c_axis=[min(gen.f.output.pseudo(:)) max(gen.f.output.pseudo(:))]; clf;
subplot(4,1,1); scatter(gen.f.pseudo_x,gen.f.pseudo_y,[],gen.f.output.pseudo,'filled');set(gca,'Ydir','reverse');caxis(c_axis);  xlim([0 100]); ylim([0 16]); colorbar('southoutside');
subplot(4,1,2); scatter(gen.f.pseudo_x,gen.f.pseudo_y,[],mean(fsim_pseudo,2),'filled');set(gca,'Ydir','reverse');caxis(c_axis); colorbar('southoutside');xlim([0 100]); ylim([0 16])
subplot(4,1,3); scatter(gen.f.pseudo_x,gen.f.pseudo_y,[],std(fsim_pseudo,[],2)./mean(fsim_pseudo,2),'filled');set(gca,'Ydir','reverse'); colorbar('southoutside');xlim([0 100]); ylim([0 16])
subplot(4,1,4); scatter(gen.f.pseudo_x,gen.f.pseudo_y,[],(mean(fsim_pseudo,2)-gen.f.output.pseudo)./gen.f.output.pseudo,'filled');set(gca,'Ydir','reverse'); colorbar('southoutside');xlim([0 100]); ylim([0 16])
% export_fig -eps 'pseudo-sec-err'

figure(23); clf; hold on; axis equal tight;
for i_real=1:parm.n_real  
    scatter(fsim_pseudo(:,i_real),gen.f.output.pseudo,'.k');
end
scatter(gen.i.output.pseudoCalc,gen.f.output.pseudo,'.r');
scatter(mean(fsim_pseudo,2),gen.f.output.pseudo,'.g');
x=[floor(min(fsim_pseudo(:))) ceil(max(fsim_pseudo(:)))];
plot(x,x,'-r'); 
plot(x,x-x*gen.i.b_wgt,'--r'); 
plot(x,x+x*gen.i.b_wgt,'--r'); 
xlabel('Apparent resistivity measured from simulated fields');
ylabel('Apparent resistivity measured from true fields');
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
% export_fig -eps 'pseudo-sec-err2'



%% Invertion on fine scale grid

i               = gen.i;
i.grid          = gen.f.grid;
i.grid.nx       = numel(i.grid.x);
i.grid.ny       = numel(i.grid.y);
i.elec_spacing  = gen.f.elec_spacing;
i.elec_id       = gen.f.elec_id;

i_fine_scale = Matlat2R2(i,gen.elec);

% save(['result/' fieldname '_i_fine_scale'],'i_fine_scale','-v7.3')
% load(['result/' fieldname '_i_fine_scale'])


%% Inversion base on true model

i               = gen.i;
i.grid          = gen.f.grid;
i.grid.nx       = numel(i.grid.x);
i.grid.ny       = numel(i.grid.y);
i.elec_spacing  = gen.f.elec_spacing;
i.elec_id       = gen.f.elec_id;

fid = fopen([gen.f.filepath 'R2_forward.dat'],'r');
A = textscan(fid,'%f %f %f %f %f %f %f');fclose(fid);
A{end-1}(2:end) = (1+.02*randn(numel(gen.f.output.resistance),1)).*gen.f.output.resistancewitherror;
fid=fopen([gen.f.filepath 'R2_forward.dat'],'w');
A2=[A{:}];
fprintf(fid,'%d\n',A2(1,1));
for u=2:size(A2,1)
    fprintf(fid,'%d %d %d %d %d %.9e %.5f\n',A2(u,:));
end
fclose(fid);

i.max_iterations = 1;
i.tolerance = 1.36;

i.rho               = mean(1000./sigma_true(:))*ones( numel(i.grid.y), numel(i.grid.x));
i.rho(i.grid.inside)= 1000./sigma_true; % f({grid_Rho.y,grid_Rho.x});

i_true_mod = Matlat2R2(i,gen.elec);

% save(['result/' fieldname '_i_true_mod'],'i_true_mod','-v7.3')
% load(['result/' fieldname '_i_true_mod'])


 
%% forward response based on the fine-scale inversion

i =  i_fine_scale;
i =  i_true_mod;

f=gen.f;
f.rho = 1000./flipud(i.output.res);
f = Matlat2R2(f,gen.elec); % write file and run forward modeling
sqrt( mean (((f.output.resistance - gen.f.output.resistancewitherror) ./ (gen.f.output.resistancewitherror)).^2 ))

Nscore.forward = @(x) (log(x./43)/1.4-gen.mu)./gen.std;
Nscore.inverse = @(x) 43.*exp(1.4*(x.*gen.std+gen.mu));
% Sec_d = Nscore.forward(1000./flipud(i.output.res));
Prim.d_full = Nscore.forward(1000./gen.f.rho);
Prim.d = reshape(Prim.d_full(gen.f.grid.inside),grid_gen.ny,grid_gen.nx);
% Prim.x = grid_gen.x; Prim.y = grid_gen.y; 
% [Prim.X, Prim.Y] = meshgrid(Prim.x, Prim.y);


ztrue = reshape(i.output.Res * Prim.d_full(:),i.grid.ny,i.grid.nx);
f.rho = 1000./Nscore.inverse(ztrue);
f = Matlat2R2(f,gen.elec); % write file and run forward modeling
sqrt( mean (((f.output.resistance - gen.f.output.resistancewitherror) ./ (gen.f.output.resistancewitherror)).^2 ))


Jrfs = reshape(sum(abs(i.output.J),2),size(Prim.d_full));
Jrtm = reshape(sum(abs(i_true_mod.output.J),2),size(Prim.d_full));

figure;
subplot(3,1,1); surf(gen.f.grid.x, gen.f.grid.y, log(Jrfs),'EdgeColor','none','facecolor','flat'); view(2); axis tight equal; set(gca,'Ydir','reverse');
subplot(3,1,2); surf(gen.f.grid.x, gen.f.grid.y, log(Jrtm),'EdgeColor','none','facecolor','flat'); view(2); axis tight equal; set(gca,'Ydir','reverse'); 
subplot(3,1,3); surf(gen.f.grid.x, gen.f.grid.y, log(Jrtm)-log(Jrfs),'EdgeColor','none','facecolor','flat'); view(2); axis tight equal; set(gca,'Ydir','reverse'); colorbar;
figure;
subplot(2,1,1); scatter(gen.f.pseudo_x,gen.f.pseudo_y,[],mean(i_fine_scale.output.J,1),'filled'); set(gca,'ydir','reverse')
subplot(2,1,2); scatter(gen.f.pseudo_x,gen.f.pseudo_y,[],mean(i_true_mod.output.J,1),'filled'); set(gca,'ydir','reverse')


%% Forward of upscaled 

filepath='.\data_gen\IO-file-inversion-scaleoftrue\';

f=gen.f;
f.grid              = gen.i.grid;    
f.elec_id           = gen.i.elec_id;

% f(Z^{est})
f.rho               = 1000./Sigma.d;
f.filepath          = ['data_gen/IO-file-forward-tomogram/'];
mkdir(f.filepath)
f                   = Matlat2R2(f,gen.elec); % write file and run forward modeling
sqrt( mean (((f.output.resistance - gen.f.output.resistancewitherror) ./ (gen.f.output.resistancewitherror)).^2 ))
sqrt(mean(((f.output.pseudo - gen.i.output.pseudoObs)./(gen.i.output.pseudoObs) ).^2))


% f( RUz^{true} )
f.rho = 1000./Nscore.inverse(Test_Sec_d);
f = Matlat2R2(f,gen.elec); % write file and run forward modeling
sqrt( mean (((f.output.resistance - gen.f.output.resistancewitherror) ./ (gen.f.output.resistancewitherror)).^2 ))
sqrt(mean(((f.output.pseudo - gen.i.output.pseudoObs)./(gen.i.output.pseudoObs) ).^2))


% f(Uz^{true})
U = zeros(numel(Sec.d), numel(Prim.d));
for i=1:numel(Sec.d)
    Res = zeros(numel(Sec.y),numel(Sec.x));
    Res(i)=1;
    F = griddedInterpolant({Sec.y,Sec.x},Res,'linear');
    res_t = F({Prim.y,Prim.x});
    U(i,:) = res_t(:) ./sum(res_t(:));
end
U(isnan(U))=0;
Uztrue = reshape(U*Prim.d(:),numel(Sec.y),numel(Sec.x));
f.rho               = 1000./Sigma.d;
f.rho(f.grid.inside) = 1000./Nscore.inverse(Uztrue(f.grid.inside));
f = Matlat2R2(f,gen.elec); % write file and run forward modeling
sqrt( mean (((f.output.resistance - gen.f.output.resistancewitherror) ./ (gen.f.output.resistancewitherror)).^2 ))
sqrt(mean(((f.output.pseudo - gen.i.output.pseudoObs)./(gen.i.output.pseudoObs) ).^2))




%% Forward true field

f=gen.f;
F = griddedInterpolant({Sigma.y,Sigma.x},1000./Sigma.d,'linear');
f.rho = F({f.grid.y, f.grid.x});
f.rho(f.grid.inside)= 1000./sigma_true;
f                   = Matlat2R2(f,gen.elec); % write file and run forward modeling
sqrt( mean (((f.output.resistance - gen.f.output.resistancewitherror) ./ (gen.f.output.resistancewitherror)).^2 ))
sqrt(mean(((f.output.pseudo - gen.i.output.pseudoObs)./(gen.i.output.pseudoObs) ).^2))



%% Figure for Synthetic schema
   
figure(199);clf; n=5;
subplot(n,1,1);imagesc(grid_gen.x, grid_gen.y, phi_true); axis tight; colormap(viridis());daspect([2 1 1])
subplot(n,1,2);imagesc(grid_gen.x, grid_gen.y, sigma_true); axis tight equal; colormap(viridis()); daspect([2 1 1])
subplot(n,1,3);scatter(gen.i.pseudo_x,gen.i.pseudo_y,[], gen.i.output.pseudo,'filled'); colormap(viridis());set(gca,'Ydir','reverse'); xlim([0 100]);ylim([0 20]);daspect([2 1 1])
subplot(n,1,4); surface(Sec.x, Sec.y, Sigma.d,'EdgeColor','none','facecolor','flat'); axis tight; colormap(viridis());set(gca,'Ydir','reverse');daspect([2 1 1])
subplot(n,1,5); imagesc(log(abs(Sigma.res))); axis equal tight; colormap(viridis());
