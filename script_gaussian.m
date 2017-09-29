clc; clear all;addpath(genpath('./.')); dbstop if error 

% Define the true case
covar.model='exponential'; covar.range0=[5 5]; covar.c0=1;covar.azimuth=0;
covar= kriginginitiaite(covar);
Prim.x=1:200;
Prim.y=1:2;
[Prim.X,Prim.Y] = meshgrid(Prim.x, Prim.y);
Prim.nxy =numel(Prim.X);
Prim.nx=numel(Prim.x);
Prim.ny=numel(Prim.y);
Prim.d = fftma_perso(covar, Prim);
Sec.dx=2;
Sec.dy=2;
Sec.x = Sec.dx/2:Sec.dx:Prim.x(end);
Sec.y = Sec.dy/2:Sec.dy:Prim.y(end);
[Sec.X,Sec.Y] = meshgrid(Sec.x, Sec.y);

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
