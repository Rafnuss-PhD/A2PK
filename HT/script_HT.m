%% Creation of the Synthetic Data
% This section gather all possible way to create the data. |gen| struct
% store the parameter and |data_generation.m| compute everything

clc; clear all; addpath('../functions','R2'); dbstop if error 

% Grid size
gen.xmax = 60; %total length in unit [m]
gen.ymax = 30; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.nx = gen.xmax*2+1;
gen.ny = gen.ymax*2+1;

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 0;          % number of well or number of point
gen.covar(1).model      = 'exponential';
gen.covar(1).range0     = [5 20]; 
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
gen.mu                  = 0.25;%0.27; % parameter of the first field. 
gen.std                 = 0.05;%.05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.f ={};
gen.Rho.f.res_matrix    = 0;
gen.Rho.elec.spacing_y  = 1; % in unit [m] | adapted to fit the fine grid
gen.Rho.elec.bufzone_y  = 5; % number of electrod to skip 
gen.Rho.elec.x_t        = 30; % in unit [m] | adapted to fit the fine grid
gen.Rho.elec.x_r        = [10 50]; % in unit [m] | adapted to fit the fine grid
gen.Rho.elec.config_max = 6000; % number of configuration of electrode maximal 
gen.Rho.i.grid.nx       = gen.nx; % | adapted to fit the fine grid
gen.Rho.i.grid.ny       = gen.ny; % log-spaced grid  | adapted to fit the fine grid
gen.Rho.i.res_matrix    = 3; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none

% Other parameter
gen.plotit              = true;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = '60x40';
gen.seed                = 8;

% Run the function
fieldname = data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);




%% Reproducing Theim equation
% Uniform K, 1 pumping well and 3 measuring well separated by 10 m each.
load('data_gen/data/FOR-60x20_2018-02-27_09-34.mat')

% Thiem equation
% 
% $$h_{2}-h_{1}={\frac {Q}{2\pi b K}}\ln \left({\frac {r_{1}}{r_{2}}}\right)$$
% 


uex = unique(gen.Rho.elec.X(gen.Rho.elec.data(:,1)));
id = bsxfun(@eq,gen.Rho.elec.X(gen.Rho.elec.data(:,1))',uex );
h = id* gen.Rho.f.output.resistance./sum(id,2); 
dh = bsxfun(@minus,h',h);

Q=1;
b=gen.ymax;
K=mean(K_true(:));
r=abs(gen.Rho.elec.x_r-gen.Rho.elec.x_t);
dr = bsxfun(@rdivide,r',r);

disp(-Q/(2*pi*b*K).*log(dr))
disp(dh)
disp(-(dh+Q/(2*pi*b*K).*log(dr))./(Q/(2*pi*b*K).*log(dr))*100)


% Point-source injection
%
% $$ Q =4\pi^2 K \frac{\partial h}{\partial r}$$
% $$h_{2}-h_{1}={\frac {Q}{4\pi K}}\ln \left({\frac {1}{r_{2}} - \frac {1}{r_{1}}}\right)$$


[~,id] = min((gen.Rho.elec.Y(gen.Rho.elec.data(:,3))-30).^2);
id = gen.Rho.elec.Y(gen.Rho.elec.data(:,3))==gen.Rho.elec.Y(gen.Rho.elec.data(id,3));
h = gen.Rho.f.output.resistance(id);

x=gen.Rho.elec.X(gen.Rho.elec.data(id,1));
y=gen.Rho.elec.Y(gen.Rho.elec.data(id,1));
x0=gen.Rho.elec.X(gen.Rho.elec.data(id,3));
y0=gen.Rho.elec.Y(gen.Rho.elec.data(id,3));

figure(2); hold on;
scatter(x, y,[],log(h),'filled')
plot(x0,y0,'xk')

Q=1;
K=mean(K_true(:));
r=sqrt( (x-x0).^2 + (y-y0).^2 );
R= 5; H=mean(h(r==R));
ht = -Q/(4*pi*K).*(-1./r+1/R) + H;

figure; hold on;
scatter3(x,y,h,'.')
scatter3(x,y,ht,'.')


%% Inversion
load('data_gen/data/GEN-60x40_2018-03-01_12-45.mat')
addpath('../functions','R2');

% Normal Score based on known distribution of Prim and Sec
Nscore.forward = @(x) ( (log10(x)+4.97)./6.66 - gen.mu)./gen.std;
Nscore.inverse = @(x) 6.66*10^(x.*gen.std+gen.mu)-4.97;
Sec=K; 
Sec.d = Nscore.forward(Sec.d);
Prim.d = Nscore.forward(K_true);
Prim.x = grid_gen.x; Prim.y = grid_gen.y; Prim.X = grid_gen.X; Prim.Y = grid_gen.Y;


% Figure of intial data
figure(1); clf;c_axis=[ min(Prim.d(:)) max(Prim.d(:)) ];
subplot(2,1,1);imagesc(grid_gen.x, grid_gen.y, Prim.d); caxis(c_axis); title('True field'); axis tight equal; 
subplot(2,1,2); hold on;
imagesc(Sec.x, Sec.y, Sec.d); title('Inverted field'); box on;colormap(viridis()); axis tight equal;set(gca,'Ydir','reverse')
plot(gen.Rho.elec.x_t,gen.Rho.elec.y,'or')
plot(gen.Rho.elec.x_r(1),gen.Rho.elec.y,'xb')
plot(gen.Rho.elec.x_r(2),gen.Rho.elec.y,'xb')

% export_fig -eps 'Prim_and_sec'



% Built the matrix G which link the true variable Prim.d to the measured coarse scale d
G = zeros(numel(Sec.d), numel(Prim.d));
for i=1:numel(Sec.d)
    Res = reshape(Sec.res(i,:)+Sec.res_out(i)/numel(Sec.res(i,:)),numel(Sec.y),numel(Sec.x));
    f = griddedInterpolant({Sec.y,Sec.x},Res,'linear');
    res_t = f({Prim.y,Prim.x});
    G(i,:) = res_t(:) ./sum(res_t(:));
end

Test_Sec_d = reshape(G * Prim.d(:), numel(Sec.y), numel(Sec.x));
figure(4); clf; colormap(viridis())
subplot(3,1,1);
surface(Sec.X,Sec.Y,Sec.d,'EdgeColor','none','facecolor','flat'); 
view(2); set(gca,'Ydir','reverse'); caxis([-1.5 1.5]); axis equal tight; box on; xlabel('x');ylabel('y'); colorbar; title('Inverted field Z^{est}')
subplot(3,1,2);surface(Sec.X,Sec.Y,Test_Sec_d,'EdgeColor','none','facecolor','flat'); 
view(2); set(gca,'Ydir','reverse');  caxis([-1.5 1.5]);axis equal tight; box on; xlabel('x');ylabel('y'); colorbar; title('G-transform of the true field Gz^{true}')
subplot(3,1,3);surface(Sec.X,Sec.Y,Test_Sec_d-Sec.d,'EdgeColor','none','facecolor','flat'); 
view(2); set(gca,'Ydir','reverse');  axis equal tight; box on; xlabel('x');ylabel('y');colorbar;  title('Error Gz^{true}-Z^{est}')
% export_fig -eps 'SecvsTestSecD'


% Generate a sampling
Prim_pt = sampling_pt(Prim,Prim.d,2,0);


% Compute the covariance of the data error
Cm = inv(sqrtm(full(gen.Rho.i.output.R(gen.Rho.i.output.inside,gen.Rho.i.output.inside))));
Cmt=(eye(size(Sec.res))-Sec.res)*Cm;


% Compute the covariance of the spatial model
covar = kriginginitiaite(gen.covar);
Cz = covar.g(squareform(pdist([Prim.X(:) Prim.Y(:)]*covar.cx)));
Cz=sparse(Cz);
Czd = Cz * G';
Cd = G * Czd;


% Combine both covariance an built the Kriging System
Cd2 = Cd+Cmt;
Czh = Cz(Prim_pt.id,Prim_pt.id);
Czzh = Cz(Prim_pt.id,:);
Czhd = Czd( Prim_pt.id ,:);
CCa = [ Cd2 Czhd' ; Czhd Czh ];
CCb = [ Czd' ; Czzh' ];


% Solve the kriging system
W=zeros(numel(Prim.x)*numel(Prim.y),numel(Sec.d(:))+Prim_pt.n);
parfor ij=1:numel(Prim.x)*numel(Prim.y)
     W(ij,:) = CCa \ CCb(:,ij);
end
% save(['ERT/result/' fieldname '_cond'],'W','Prim_pt','G','Nscore','Sec','Prim')

% Compute the Kriging map
zh = reshape( W * [Sec.d(:) ; Prim_pt.d], ny, nx);
%zhtest = reshape( W * [Test_Sec_d(:) ; Prim_pt.d], ny, nx);

figure(5); clf;   colormap(viridis())
surf(Prim.x,Prim.y,zh,'EdgeColor','none','facecolor','flat'); caxis([-3 3])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); colorbar('southoutside'); title('Kriging Estimate')
export_fig -eps 'Krig'


%% Simulation of the Area-to-point Kriging 
rng('shuffle');

parm.n_real = 500;
zcs=nan(ny, nx,parm.n_real);
z=nan(numel(Sec.y), numel(Sec.x),parm.n_real);
for i_real=1:parm.n_real
    zs = fftma_perso(covar, struct('x',Prim.x,'y',Prim.y));
    %zs = fftma_perso(gen.covar, grid_gen);
    zhs = W * [G * zs(:) ; zs(Prim_pt.id)./1.3];
    r = zh(:) + (zs(:) - zhs(:));
    zcs(:,:,i_real) = reshape( r, ny, nx);
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
export_fig -eps 'PrimOverview'

figure(7);clf; colormap(viridis())
c_axis=[ -3 3];
subplot(4,1,1);surf(Sec.x, Sec.y, Sec.d,'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis tight equal; set(gca,'Ydir','reverse');  colorbar;
subplot(4,1,2);surf(Sec.x, Sec.y, z(:,:,1),'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis tight equal; set(gca,'Ydir','reverse'); 
subplot(4,1,3);surf(Sec.x, Sec.y, mean(z,3),'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis tight equal; set(gca,'Ydir','reverse');
subplot(4,1,4);surf(Sec.x, Sec.y, std(z,[],3),'EdgeColor','none','facecolor','flat');  title('d Average of True field');view(2); axis tight equal; set(gca,'Ydir','reverse'); colorbar;
export_fig -eps 'SecOverview'

figure(8);clf;colormap(viridis())
subplot(2,1,1);surface(Sec.X,Sec.Y,Test_Sec_d-Sec.d,'EdgeColor','none','facecolor','flat'); view(2); set(gca,'Ydir','reverse');  axis equal tight; box on; xlabel('x');ylabel('y'); caxis([-.6 .6]);colorbar('southoutside');
subplot(2,1,2);surf(Sec.x, Sec.y, mean(z,3)-Sec.d,'EdgeColor','none','facecolor','flat'); view(2); axis tight equal; set(gca,'Ydir','reverse'); colorbar('southoutside');caxis([-.6 .6])
export_fig -eps 'GztrueGzsim'

% Compute the Variogram and Histogram of realiaztions
parm.n_real=500;
vario_x=nan(parm.n_real,nx);
vario_y=nan(parm.n_real,ny);
for i_real=1:parm.n_real
    r = zcs(:,:,i_real);
    %r = (r(:)-mean(r(:)))./std(r(:));
    [vario_x(i_real,:),vario_y(i_real,:)]=variogram_gridded_perso(reshape( r(:), ny, nx));
end
[vario_prim_x,vario_prim_y]=variogram_gridded_perso(Prim.d);

figure(9);clf;
subplot(2,1,1);  hold on; 
h1=plot(Prim.x(1:2:end),vario_x(:,1:2:end)','b','color',[.5 .5 .5]);
h2=plot(Prim.x,vario_prim_x,'-r','linewidth',2);
h3=plot(Prim.x,1-covar.g(Prim.x*covar.cx(1)),'--k','linewidth',2);
xlim([0 60]); xlabel('Lag-distance h_x ');ylabel('Variogram \gamma(h_x)')
legend([h1(1) h2 h3],'500 realizations','True field','Theorical Model')
subplot(2,1,2); hold on; 
h1=plot(Prim.y,vario_y','b','color',[.5 .5 .5]);
h2=plot(Prim.y,vario_prim_y,'-r','linewidth',2);
h3=plot(Prim.y,1-covar.g(Prim.y*covar.cx(4)),'--k','linewidth',2);
xlim([0 6]); xlabel('Lag-distance h_y ');ylabel('Variogram \gamma(h_y)')
legend([h1(1) h2 h3],'500 realizations','True field','Theorical Model')
export_fig -eps 'Vario'


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

parm.n_real=500;
fsim_pseudo=nan(numel(gen.Rho.f.output.pseudo),parm.n_real);
fsim_resistance=nan(numel(gen.Rho.f.output.resistance),parm.n_real);
rho = 1000./Nscore.inverse(zcs);
for i_real=1:parm.n_real
    f={};
    f.res_matrix        = gen.Rho.f.res_matrix;
    f.grid              = gen.Rho.f.grid;    
    f.header            = 'Forward';  % title of up to 80 characters
    f.job_type          = 0;
    f.filepath          = ['data_gen/IO-file-' num2str(i_real) '/'];
    f.readonly          = 0;
    f.alpha_aniso       = gen.Rho.f.alpha_aniso;
    f.elec_spacing      = gen.Rho.f.elec_spacing;
    f.elec_id           = gen.Rho.f.elec_id;
    f.rho               = rho(:,:,i_real);
    f.num_regions       = gen.Rho.f.num_regions;
    f.rho_min           = gen.Rho.f.rho_min;
    f.rho_avg           = gen.Rho.f.rho_avg;
    f.rho_max           = gen.Rho.f.rho_max;
    
    mkdir(f.filepath)
    f                   = Matlat2R2(f,gen.Rho.elec); % write file and run forward modeling
    fsim_pseudo(:,i_real) = f.output.pseudo;
    fsim_resistance(:,i_real) = f.output.resistance;
end


% Compute the misfit
fsim_misfit=nan(numel(gen.Rho.f.output.resistancewitherror),parm.n_real);
err=nan(1,parm.n_real);
for i_real=1:parm.n_real  
    fsim_misfit(:,i_real) = (fsim_resistance(:,i_real) - gen.Rho.f.output.resistancewitherror) ./ (gen.Rho.i.a_wgt + gen.Rho.i.b_wgt*gen.Rho.f.output.resistancewitherror);
    err(i_real) = sqrt(mean(fsim_misfit(:,i_real).^2));
end
% export_fig -eps 'Histogram_of_misfit'
% save(['result/' fieldname '_sim'],'zcs','fsim_pseudo','fsim_resistance')


figure(11);
histogram(err);
xlabel('Misfit'); ylabel('Histogram')
% export_fig -eps 'misfit-hist'




figure(12);clf; colormap(viridis());c_axis=[min(gen.Rho.f.output.pseudo(:)) max(gen.Rho.f.output.pseudo(:))]; clf;
subplot(3,1,1); scatter(gen.Rho.f.pseudo_x,gen.Rho.f.pseudo_y,[],gen.Rho.f.output.pseudo,'filled');set(gca,'Ydir','reverse');caxis(c_axis);  xlim([0 100]); ylim([0 16]); colorbar('southoutside');
subplot(3,1,2); scatter(gen.Rho.f.pseudo_x,gen.Rho.f.pseudo_y,[],mean(fsim_pseudo,2),'filled');set(gca,'Ydir','reverse');caxis(c_axis); colorbar('southoutside');xlim([0 100]); ylim([0 16])
subplot(3,1,3); scatter(gen.Rho.f.pseudo_x,gen.Rho.f.pseudo_y,[],std(fsim_pseudo,[],2)./mean(fsim_pseudo,2),'filled');set(gca,'Ydir','reverse'); colorbar('southoutside');xlim([0 100]); ylim([0 16])
% export_fig -eps 'pseudo-sec'

figure(13);clf;colormap(viridis()); c_axis=[min(gen.Rho.f.output.pseudo(:)) max(gen.Rho.f.output.pseudo(:))]; clf;
subplot(3,1,1); scatter(gen.Rho.f.pseudo_x,gen.Rho.f.pseudo_y,[],(mean(fsim_pseudo,2)-gen.Rho.f.output.pseudo)./gen.Rho.f.output.pseudo,'filled');set(gca,'Ydir','reverse'); colorbar('southoutside');xlim([0 100]); ylim([0 16])
caxis([-.1 .1])
subplot(3,1,2); scatter(gen.Rho.f.pseudo_x,gen.Rho.f.pseudo_y,[],(gen.Rho.i.output.pseudo-gen.Rho.f.output.pseudo)./gen.Rho.f.output.pseudo,'filled');set(gca,'Ydir','reverse'); colorbar('southoutside');xlim([0 100]); ylim([0 16])
caxis([-.05 .05])
% export_fig -eps 'pseudo-sec-err'

figure(23); clf; hold on; axis equal tight;
for i_real=1:parm.n_real  
    scatter(fsim_pseudo(:,i_real),gen.Rho.f.output.pseudo,'.k');
end
scatter(gen.Rho.i.output.pseudo,gen.Rho.f.output.pseudo,'.r');
scatter(mean(fsim_pseudo,2),gen.Rho.f.output.pseudo,'.g');
x=[floor(min(fsim_pseudo(:))) ceil(max(fsim_pseudo(:)))];
plot(x,x,'-r'); 
plot(x,x-x*3*gen.Rho.i.b_wgt,'--r'); 
plot(x,x+x*3*gen.Rho.i.b_wgt,'--r'); 
xlabel('Apparent resistivity measured from simulated fields');
ylabel('Apparent resistivity measured from true fields');
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
% export_fig -eps 'pseudo-sec-err'

   

%% Figure for Synthetic schema
   
figure(199);clf; n=5;
subplot(n,1,1);imagesc(grid_gen.x, grid_gen.y, phi_true); axis tight; colormap(viridis());daspect([2 1 1])
subplot(n,1,2);imagesc(grid_gen.x, grid_gen.y, sigma_true); axis tight equal; colormap(viridis()); daspect([2 1 1])
subplot(n,1,3);scatter(gen.Rho.i.pseudo_x,gen.Rho.i.pseudo_y,[], gen.Rho.i.output.pseudo,'filled'); colormap(viridis());set(gca,'Ydir','reverse'); xlim([0 100]);ylim([0 20]);daspect([2 1 1])
subplot(n,1,4); surface(Sec.x, Sec.y, Sigma.d,'EdgeColor','none','facecolor','flat'); axis tight; colormap(viridis());set(gca,'Ydir','reverse');daspect([2 1 1])
subplot(n,1,5); imagesc(log(abs(Sigma.res))); axis equal tight; colormap(viridis());
