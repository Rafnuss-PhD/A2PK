%% DATA CREATION
clc; clear all; addpath(genpath('./.')); dbstop if error 
% This section gather all possible way to create the data. |gen| struct
% store the parameter and |data_generation.m| compute everything

% Grid size
gen.xmax = 100; %total length in unit [m]
gen.ymax = 20; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.nx = 800;
gen.ny = 40;

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
gen.covar(1).range0     = [5 50];
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
gen.mu                  = -1.579;%0.27; % parameter of the first field. 
gen.std                 = 0.22361;%.05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.f ={};
gen.Rho.f.res_matrix    = 0;
gen.Rho.elec.spacing    = 8; % in grid spacing unit of the forward grid
gen.Rho.elec.config_max = 6000; % number of configuration of electrode maximal 
gen.Rho.i.grid.nx       = 200;
gen.Rho.i.grid.ny       = 20; % log-spaced grid.
gen.Rho.i.res_matrix    = 3; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = '800x40';
gen.seed                = 123;

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);


%%
clear all; addpath(genpath('./.')); dbstop if error 
load('result-A2PK/GEN-800x40_2017-11-02_11-32'); %alpha = 27339.756;

%addpath('C:\Users\Raphael\Documents\MATLAB\ColorBrewer_ Attractive and Distinctive Colormaps\DrosteEffect-BrewerMap-aee8a46\')
%colormap(brewermap([],'YlOrBr'))
addpath('C:\Users\Raphael\Documents\MATLAB\Colormaps\Colormaps (5)\Colormaps\')
colormap(viridis())

[kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
parm.nscore=1;
Nscore = nscore(kern, parm, 0);
Sec=Sigma; 
Sec.d = reshape(Nscore.forward(Sec.d(:)) ,numel(Sec.y),numel(Sec.x));
Prim.d = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);
Prim.x = grid_gen.x; Prim.y = grid_gen.y; Prim.X = grid_gen.X; Prim.Y = grid_gen.Y;


figure(10); clf;c_axis=[ min(Prim.d(:)) max(Prim.d(:)) ];
subplot(2,1,1);imagesc(grid_gen.x, grid_gen.y, Prim.d); caxis(c_axis); title('zt True field');
subplot(2,1,2); imagesc(Sec.x, Sec.y, Sec.d); caxis(c_axis); title('Inverted field'); axis tight equal; box on;colormap(viridis())


figure(101); clf; c_axis_n=log10([ min(sigma_true(:)) max(sigma_true(:)) ]);
subplot(2,1,1);surface(grid_gen.x, grid_gen.y, log10(sigma_true),'EdgeColor','none','facecolor','flat'); 
caxis(c_axis_n); title('True Electrical Conductivity \sigma^{true}');axis equal tight; box on; xlabel('x');ylabel('y'); set(gca,'Ydir','reverse');
subplot(2,1,2); surface(Sec.x, Sec.y, log10(Sigma.d),'EdgeColor','none','facecolor','flat'); 
caxis(c_axis_n); title('Inverted Electrical Conductivity U \sigma^{est}'); axis equal tight; box on; xlabel('x');ylabel('y');set(gca,'Ydir','reverse'); colorbar('southoutside')
colormap(viridis())

% Built the matrix G which link the true variable Prim.d to the measured coarse scale d
G = zeros(numel(Sec.d), numel(Prim.d));
for i=1:numel(Sec.d)
    Res = reshape(Sigma.res(i,:)+Sigma.res_out(i)/numel(Sigma.res(i,:)),numel(Sec.y),numel(Sec.x));
    f = griddedInterpolant({Sec.y,Sec.x},Res,'linear');
    res_t = f({Prim.y,Prim.x});
    G(i,:) = res_t(:) ./sum(res_t(:));
end



i=[905 1838 3855]; th=.05;
figure(1); clf; colormap(viridis())
subplot(3,numel(i),[1 numel(i)]);hold on; surface(Sec.X,Sec.Y,reshape(log(diag(Sigma.res)),numel(Sec.y),numel(Sec.x)),'EdgeColor','none','facecolor','flat');  scatter3(Sec.X(i),Sec.Y(i),100*ones(numel(i),1),'sr','filled')
view(2); axis tight equal; set(gca,'Ydir','reverse'); box on; axis equal tight; xlabel('x');ylabel('y'); title('Diagonal of the resolution matrix R'); colorbar;
for ii=1:numel(i)
    subplot(3,numel(i),numel(i)+ii);hold on; surface(Sec.X,Sec.Y,reshape(Sigma.res(i(ii),:),numel(Sec.y),numel(Sec.x)),'EdgeColor','none','facecolor','flat');  scatter3(Sec.X(i(ii)),Sec.Y(i(ii)),100,'sr','filled')
    view(2); axis tight equal; set(gca,'Ydir','reverse'); box on; axis equal tight; xlabel('x');ylabel('y'); title('Resolution of the red dot R(i,:)')
    subplot(3,numel(i),2*numel(i)+ii);hold on; surface(Prim.X,Prim.Y,reshape(G(i(ii),:), numel(Prim.y), numel(Prim.x)),'EdgeColor','none','facecolor','flat'); scatter3(Sec.X(i(ii)),Sec.Y(i(ii)),100,'sr','filled')
    view(2); axis tight equal; set(gca,'Ydir','reverse'); box on;  axis equal tight; colorbar('southoutside'); xlabel('x');ylabel('y')
end


% Compute coarse scale d
Test_Sec_d = reshape(G * Prim.d(:), numel(Sec.y), numel(Sec.x));
figure(2); clf; colormap(viridis())
subplot(2,1,1);surface(Sec.X,Sec.Y,Sec.d,'EdgeColor','none','facecolor','flat'); view(2); set(gca,'Ydir','reverse'); caxis([-3 3]); axis equal tight; box on; xlabel('x');ylabel('y'); title('Inverted Electrical Conductivity \Sigma^{est}')
subplot(2,1,2);surface(Sec.X,Sec.Y,Test_Sec_d,'EdgeColor','none','facecolor','flat'); view(2); set(gca,'Ydir','reverse');  caxis([-3 3]);axis equal tight; box on; xlabel('x');ylabel('y'); colorbar('southoutside'); title('G-transform of the True Electrical Conductivity Gz^{true}')


% Simulate sampling
Prim_pt = sampling_pt(Prim,Prim.d,1,1);


% Compute the covariance of the data error

Cm = inv(sqrtm(full(gen.Rho.i.output.R(gen.Rho.i.output.inside,gen.Rho.i.output.inside))));
Cmt=(eye(size(Sigma.res))-Sigma.res)*Cm;


% Subsample D
% k=floor(trace(Sigma.res));
% pt=net(haltonset(2),k);
% [~,idx] = min(bsxfun(@minus,Sec.X(:),pt(:,1)'.*Sec.x(end)).^2 + bsxfun(@minus,Sec.Y(:),pt(:,2)'.*Sec.y(end)).^2);
% idx=unique(idx);
% Sec.d(:) = Sec.d(idx)';
% G=G(idx,:);



% Simulation 
covar = kriginginitiaite(gen.covar);
[X, Y] = meshgrid(Prim.x, Prim.y); X=X(:);Y=Y(:);
nx=numel(Prim.x); ny=numel(Prim.y);


Cz = covar.g(squareform(pdist([X Y]*covar.cx)));
%disp('Cz computed')
Cz=sparse(Cz);
%disp('Cz sparse computed')

Czd = Cz * G';
disp('Czd computed')

Cd = G * Czd;
disp('Cd computed')


Cd2 = Cd+Cmt;
Czh = covar.g(squareform(pdist([Prim_pt.x(:) Prim_pt.y(:)]*covar.cx)));% Czh =covar.g(0);
Czzh = covar.g(pdist2([X Y]*covar.cx,[Prim_pt.x(:) Prim_pt.y(:)]*covar.cx));
Czhd = Czd( Prim_pt.id ,:);

CCa = [ Cd2 Czhd' ; Czhd Czh ];
CCb = [ Czd' ; Czzh' ];
disp('Covariances computed')

% warning('off','all')
W=zeros(nx*ny,numel(Sec.d(:))+Prim_pt.n);
parfor ij=1:nx*ny
     W(ij,:) = CCa \ CCb(:,ij);
end
disp('W computed')
save('result-A2PK/GEN-400x40_2017-10-30_17-53_cond','W','Prim_pt')

zh = reshape( W * [Sec.d(:) ; Prim_pt.d], ny, nx);
%zhtest = reshape( W * [Test_Sec_d(:) ; Prim_pt.d], ny, nx);
figure(8); clf;   colormap(viridis())
surf(Prim.x,Prim.y,zh,'EdgeColor','none','facecolor','flat'); caxis([-3 3])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); colorbar('southoutside'); title('Kriging Estimate')

%subplot(2,1,2);surf(Prim.x,Prim.y,zhtest,'EdgeColor','none'); view(2); axis tight equal; set(gca,'Ydir','reverse'); box on; colorbar;



parm.n_real = 500;

rng('shuffle');
zcs=nan(ny, nx,parm.n_real);
z=nan(numel(Sec.y), numel(Sec.x),parm.n_real);
for i_real=1:parm.n_real
    zs = fftma_perso(covar, struct('x',Prim.x,'y',Prim.y));
    zhs = W * [G * zs(:) ; zs(Prim_pt.id)];
    r = zh(:) + (zs(:) - zhs(:));
    zcs(:,:,i_real) = reshape( (r-mean(r))./std(r), ny, nx);
    z(:,:,i_real)=reshape(G*r(:), numel(Sec.y), numel(Sec.x) );
end


figure(10);clf; colormap(viridis())
c_axis=[ -3 3];
subplot(8,1,1);surf(Prim.x, Prim.y, Prim.d,'EdgeColor','none','facecolor','flat'); caxis(c_axis); title('True Electrical Conductivity');view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
% hold on; scatter(Prim_pt.x,Prim_pt.y,'filled','r')
subplot(8,1,2);surf(Prim.x, Prim.y, zcs(:,:,1),'EdgeColor','none','facecolor','flat'); caxis(c_axis);view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
subplot(8,1,3);surf(Prim.x, Prim.y, mean(zcs,3),'EdgeColor','none','facecolor','flat'); caxis(c_axis);view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
subplot(8,1,4);surf(Prim.x, Prim.y, std(zcs,[],3),'EdgeColor','none','facecolor','flat'); view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
subplot(8,1,5);surf(Sec.x, Sec.y, Sec.d,'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
subplot(8,1,6);surf(Sec.x, Sec.y, z(:,:,1),'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
subplot(8,1,7);surf(Sec.x, Sec.y, mean(z,3),'EdgeColor','none','facecolor','flat'); caxis(c_axis); view(2); axis tight equal; set(gca,'Ydir','reverse'); box on
subplot(8,1,8);surf(Sec.x, Sec.y, std(z,[],3),'EdgeColor','none','facecolor','flat');  title('d Average of True field');view(2); axis tight equal; set(gca,'Ydir','reverse'); box on




figure(11);clf; hold on; colormap(viridis())
for i_real=1:parm.n_real
   vario=variogram_gridded_perso(zcs(:,:,i_real));
   h1=plot(Prim.x,vario,'b','color',[.5 .5 .5]);
end
vario=variogram_gridded_perso(Prim.d);
h2=plot(Prim.x,vario,'-r','linewidth',2);
h3=plot(Prim.x,1-covar.g(Prim.x*covar.cx(1)),'--k','linewidth',2);
xlim([0 70]); xlabel('Lag-distance h ');ylabel('Variogram \gamma(h)')
legend([h1 h2 h3],'500 realizations','True field','Theorical Model')



%% Forward simulation

for i_real=1:parm.n_real
    z1=zcs(:,:,i_real);
    s1 = reshape(Nscore.inverse(z1(:)) ,numel(Prim.y),numel(Prim.x));
    

f = gen.Rho.f;
f.res_matrix = gen.Rho.f.res_matrix;
f.grid.x            = gen.Rho.f.x;
f.grid.y            = gen.Rho.f.y;
f.grid.x_n          = gen.Rho.f.x_n;
f.grid.y_n          = gen.Rho.f.y_n;


% Forward
f.header            = 'Forward';  % title of up to 80 characters
f.job_type          = 0;
f.filepath          = filepath;
f.readonly          = 0;
f.alpha_aniso       = gen.covar.range0(2)/gen.covar.range0(1);

% Rho value
% f                  = griddedInterpolant({grid.y,grid.x},rho_true,'nearest','nearest');
f.rho               = rho_true; % f({grid_Rho.y,grid_Rho.x});
% f.dmin.filename       = 'gtrue.dat';
f.num_regions       = 1+numel(f.rho);
f.rho_min           = min(rho_true(:));
f.rho_avg           = mean(rho_true(:));
f.rho_max           = max(rho_true(:))*2;
f                   = Matlat2R2(f,gen.Rho.elec); % write file and run forward modeling

    
    
    
    
    Matlat2R2(gen.Rho.grid, dmin, gen.Rho.elec);
end

