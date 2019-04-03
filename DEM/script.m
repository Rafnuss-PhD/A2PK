%% Load and define variable
load('stuffForMathieu')

xfs = 1:size(trainingFine,1);
yfs = 1:size(trainingFine,2);
[Xfs, Yfs] = meshgrid(xfs, yfs);
dx=size(trainingFine,1)/size(trainingCoarse,1);
dy=size(trainingFine,2)/size(trainingCoarse,2);

xcs = linspace(1.5,size(trainingFine,1)-.5,size(trainingCoarse,1));
ycs = linspace(1.5,size(trainingFine,2)-.5,size(trainingCoarse,2));
[Xcs, Ycs] = meshgrid(xcs,ycs);

G = zeros(numel(trainingCoarse), numel(trainingFine));
for ij=1:numel(trainingFine)
     [~,id_z]=min((Xfs(ij)-Xcs(:)).^2 + (Yfs(ij)-Ycs(:)).^2);
     G(id_z,ij)=1;%/(Z.dx*Z.dy);
end
for ij = 1:numel(trainingCoarse)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end

% Check upscaling in intial space (non-transformed)
figure;
subplot(1,3,1); imagesc(trainingCoarse); caxis([-25 25])
tmp = reshape(G*trainingFine(:),size(trainingCoarse,1),[]);
subplot(1,3,2); imagesc(tmp); caxis([-25 25])
subplot(1,3,3); imagesc((tmp-trainingCoarse)/std(trainingCoarse(:))); caxis([-.1 .15])


%% Compute the transform
[cdf,xi] = ksdensity(trainingFine(:),'NumPoints',1000,'Function','cdf');

T_F = griddedInterpolant(xi,cdf,'pchip','pchip');
Tinv_F = griddedInterpolant(cdf,xi,'pchip','pchip');
Phi_fs_inv = @(y) Tinv_F(normcdf(y)); % back-transform a value in normal space by taking the normcdf.
Phi_fs = @(y) norminv( T_F(y) );

trainingFineTrans = Phi_fs(trainingFine);


[cdf1,xi1] = ksdensity(G*trainingFineTrans(:),'NumPoints',1000,'Function','cdf');
[cdf2,xi2] = ksdensity(trainingCoarse(:),'NumPoints',100,'Function','cdf');

T_F1 = griddedInterpolant(xi1,cdf1,'pchip','pchip');
Tinv_F1 = griddedInterpolant(cdf1,xi1,'pchip','pchip');
T_F2 = griddedInterpolant(xi2,cdf2,'pchip','pchip');
Tinv_F2 = griddedInterpolant(cdf2,xi2,'pchip','pchip');

Phi_cs = @(y) Tinv_F1(T_F2(y)); % back-transform a value in normal space by taking the normcdf.
Phi_cs_inv = @(y) Tinv_F2( T_F1(y) );

trainingCoarseTrans = Phi_cs(trainingCoarse);

% Check transformation
figure;
hold on;
plot(trainingCoarse(:), G*trainingFineTrans(:),'.k')
plot(trainingCoarse(:), Phi_cs(trainingCoarse(:)),'.r')

figure; hold on;
plot(G*trainingFineTrans(:), trainingCoarse(:),'.k')
plot(G*trainingFineTrans(:), Phi_cs_inv(G*trainingFineTrans(:)),'.r')

% Check upscaling
plot(G*trainingFineTrans(:),trainingCoarseTrans(:),'.k')

figure;
subplot(1,3,1); imagesc(trainingCoarseTrans); caxis([-4 4])
tmp = reshape(G*trainingFineTrans(:),size(trainingCoarse,1),[]);
subplot(1,3,2); imagesc(tmp); caxis([-4 4])
subplot(1,3,3); imagesc((tmp-trainingCoarseTrans)/std(trainingCoarseTrans(:)));  caxis([-.1 .15])

%% Learn covariance
addpath('../functions')
addpath('../FastGaussianSimulation/')
[gamma_x, gamma_y] = variogram_gridded_perso(trainingFineTrans);

y = mean([gamma_x(1:51) gamma_y(1:51)],2);
fun1 = @(x) fun(x,xi,y);
x0=[1 1 1];
xmin = fmincon(fun1, x0,[],[],[],[],[0.001 .9 0.1],[100 1.1 10],[]);

covar.model = 'k-bessel';
covar.alpha = xmin(1);
covar.var = xmin(2);
covar.range = [xmin(3) xmin(3)];
covar = covarIni(covar);

covar2.model = 'spherical';
covar2.var = 1;
covar2.range = [7 7];
covar2 = covarIni(covar2);

xi=0:50;
figure; hold on;
plot(xi,gamma_x(1:51));
plot(xi,gamma_y(1:51));
plot(xi,1-covar.g(xi*covar.cx(1)),'linewidth',2);
plot(xi,1-covar2.g(xi*covar2.cx(1)),'linewidth',2);



%% Simulation A2PK
addpath('..')

xfs = 1:size(targetFine,1);
yfs = 1:size(targetFine,2);
[Xfs, Yfs] = meshgrid(xfs, yfs);
xcs = linspace(1.5,size(targetFine,1)-.5,size(targetCoarse,1));
ycs = linspace(1.5,size(targetFine,2)-.5,size(targetCoarse,2));
[Xcs, Ycs] = meshgrid(xcs,ycs);

G = zeros(numel(targetCoarse), numel(targetFine));
for ij=1:numel(targetFine)
     [~,id_z]=min((Xfs(ij)-Xcs(:)).^2 + (Yfs(ij)-Ycs(:)).^2);
     G(id_z,ij)=1;%/(Z.dx*Z.dy);
end
for ij = 1:numel(targetCoarse)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end

targetCoarseTrans = Phi_cs(targetCoarse);
n_real = 5;

[targetfineTransSim,targetfineTransEst] = A2PK(xfs,yfs,[],targetCoarseTrans,G,covar,n_real);

targetfineEst = Phi_fs_inv(targetfineTransEst);
targetfineSim = Phi_fs_inv(targetfineTransSim);



%% Comparing

figure; colormap gray
subplot(2,2,1); imagesc(targetFine); caxis([-25 25]); axis square;
subplot(2,2,2); imagesc(targetfineSim(:,:,1)); caxis([-25 25]); axis square;
subplot(2,2,3); imagesc(targetfineSim(:,:,2)); caxis([-25 25]); axis square;
subplot(2,2,4); imagesc(targetfineSim(:,:,3)); caxis([-25 25]); axis square;

figure;
imagesc((targetfineEst-targetFine)./targetFine); axis square;

[gamma_x, gamma_y] = variogram_gridded_perso(targetfineSim(:,:,1));
[gamma_x2, gamma_y2] = variogram_gridded_perso(targetFine);
[gamma_x3, gamma_y3] = variogram_gridded_perso(trainingFine);

xi=0:50;
figure; hold on;
plot(xi,gamma_x(1:51));
plot(xi,gamma_y(1:51));
plot(xi,gamma_x2(1:51));
plot(xi,gamma_y2(1:51));
plot(xi,gamma_x3(1:51));
plot(xi,gamma_y3(1:51));
% plot(xi,1-covar.g(xi*covar.cx(1)),'linewidth',2);


