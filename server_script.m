clear all; addpath(genpath('./.')); dbstop if error 
load('result-A2PK/GEN-120x20_2017-10-02_17-25.mat'); 

[kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
parm.nscore=1;
Nscore = nscore(kern, parm, 0);
Sec.x=Sigma.x_raw; Sec.y=Sigma.y_raw; [Sec.X , Sec.Y] = meshgrid(Sigma.x_raw, Sigma.y_raw); 
Sec.d = reshape(Nscore.forward(Sigma.d_raw(:)) ,numel(Sec.y),numel(Sec.x));
Prim.d = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);
Prim.x = grid_gen.x; Prim.y = grid_gen.y; Prim.X = grid_gen.X; Prim.Y = grid_gen.Y;

Sigma.res = Sigma.res +  repmat(Sigma.res_out./numel(Sigma.d_raw),1,numel(Sigma.d_raw));

% Built the matrix G which link the true variable Prim.d to the measured coarse scale d
G = zeros(numel(Sec.d), numel(Prim.d));
for i=1:numel(Sec.d)
    Res = reshape(Sigma.res(i,:)+Sigma.res_out(i)/numel(Sigma.res(i,:)),numel(Sec.y),numel(Sec.x));
    f = griddedInterpolant({Sec.y,Sec.x},Res,'nearest');
    res_t = f({Prim.y,Prim.x});
    G(i,:) = res_t(:) ./sum(res_t(:));
end

i=10;
figure(1); clf; 
subplot(2,1,1);surface(Sec.X,Sec.Y,reshape(Sigma.res(i,:),numel(Sec.y),numel(Sec.x)),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; hold on; scatter3(Sec.X(i),Sec.Y(i),100,'filled','r')
subplot(2,1,2);surface(Prim.X,Prim.Y,reshape(G(i,:), numel(Prim.y), numel(Prim.x)),'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; hold on; scatter3(Sec.X(i),Sec.Y(i),100,'filled','r')


% Compute coarse scale d
Test_Sec_d = reshape(G * Prim.d(:), numel(Sec.y), numel(Sec.x));
figure(4); clf;  c_axis=[ min(Prim.d(:)) max(Prim.d(:))];
subplot(3,1,1);surf(Prim.x,Prim.y,Prim.d,'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis(c_axis); 
subplot(3,1,2);surf(Sec.x,Sec.y,Sec.d,'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis(c_axis); 
subplot(3,1,3);surf(Sec.x,Sec.y,Test_Sec_d,'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; caxis(c_axis); 



% Simulate sampling
prim = sampling_pt(Prim,Prim.d,1,2);



% Compute the covariance of the data error
inside=false( numel(gen.Rho.i.yy)-1, numel(gen.Rho.i.xx)-1);
dx = (numel(gen.Rho.i.xx)-1-numel(gen.Rho.grid.x))/2 +1;
inside( 1:numel(gen.Rho.grid.y), dx:end-dx+1) = true;
Cm = inv(sqrtm(full(gen.Rho.i.output.R(inside,inside))));
Cmt=(eye(size(Sigma.res))-Sigma.res)*Cm;


% Subsample D
% k=floor(trace(Sigma.res));
% pt=net(haltonset(2),k);
% [~,idx] = min(bsxfun(@minus,Sec.X(:),pt(:,1)'.*Sec.x(end)).^2 + bsxfun(@minus,Sec.Y(:),pt(:,2)'.*Sec.y(end)).^2);
% idx=unique(idx);
% D = Sec.d(idx)';
% G2=G(idx,:);
D=Sec.d(:);
G2=G;

% Simulation 
covar = kriginginitiaite(gen.covar);
parm.seed_U='shuffle';
parm.n_real = 10;
wradius = 2;
[X, Y] = meshgrid(Prim.x, Prim.y); X=X(:);Y=Y(:);
nx=numel(Prim.x); ny=numel(Prim.y);



%Cz = covar.g(squareform(pdist([X Y]*covar.cx)));
%disp('Cz computed')

%Cz=sparse(Cz);
%disp('Cz sparse computed')

Czd = Cz * G2';
disp('Czd computed')

Cd = G2 * Czd;
disp('Cd computed')


Cd2 = Cd+Cmt;


Czh = covar.g(squareform(pdist([prim.x(:) prim.y(:)]*covar.cx)));% Czh =covar.g(0);
Czzh = covar.g(pdist2([X Y]*covar.cx,[prim.x(:) prim.y(:)]*covar.cx));

Czhd = Czd( prim.id ,:);

CCa = [ Cd2 Czhd' ; Czhd Czh ];
CCb = [ Czd' ; Czzh' ];

disp('Covariances computed')

% warning('off','all')
W=zeros(nx*ny,numel(D)+prim.n);
parfor ij=1:nx*ny
     W(ij,:) = CCa \ CCb(:,ij);
end

disp('W computed')


zh = reshape( W * [D(:) ; prim.d], ny, nx);

zhtest = reshape( W * [Test_Sec_d(:) ; prim.d], ny, nx);
figure(8); clf;  
subplot(2,1,1);surf(Prim.x,Prim.y,zh,'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; colorbar;
subplot(2,1,2);surf(Prim.x,Prim.y,zhtest,'EdgeColor','none'); view(2); axis tight; set(gca,'Ydir','reverse'); box on; colorbar;





rng(parm.seed_U);
zcs=nan(ny, nx,parm.n_real);
for i_real=1:parm.n_real
    zs = fftma_perso(covar, struct('x',Prim.x,'y',Prim.y));
    zhs = W * [G2 * zs(:) ; prim.d];
    r = zhtest(:) + (zs(:) - zhs(:));
    zcs(:,:,i_real) = reshape( (r-mean(r))./std(r), ny, nx);
end


figure(10);clf
c_axis=[ min(Prim.d(:)) max(Prim.d(:))];
subplot(4,1,1);surf(Prim.x, Prim.y, Prim.d,'EdgeColor','none'); caxis(c_axis); title('Prim.d True field');;view(2); axis tight; set(gca,'Ydir','reverse'); box on
subplot(4,1,2);surf(Sec.x, Sec.y, Sec.d,'EdgeColor','none'); caxis(c_axis); title('d Average of True field');view(2); axis tight; set(gca,'Ydir','reverse'); box on
subplot(4,1,3);surf(Prim.x, Prim.y, zcs(:,:,1),'EdgeColor','none'); caxis(c_axis);view(2); axis tight; set(gca,'Ydir','reverse'); box on
z=nan(numel(Sec.y), numel(Sec.x),parm.n_real);
for i_real=1:parm.n_real
    r=zcs(:,:,i_real);
   z(:,:,i_real)=reshape(G*r(:), numel(Sec.y), numel(Sec.x) );
end
subplot(4,1,4);surf(Sec.x, Sec.y, mean(z,3),'EdgeColor','none'); caxis(c_axis); title('d Average of True field');view(2); axis tight; set(gca,'Ydir','reverse'); box on
% hold on; scatter(Sec.X(idx),Sec.Y(idx),'filled','r')

figure(101);clf
c_axis=[ min(Prim.d(:)) max(Prim.d(:))];
subplot(4,1,1);surf(Prim.x, Prim.y, Prim.d,'EdgeColor','none'); caxis(c_axis); title('Prim.d True field');;view(2); axis tight; set(gca,'Ydir','reverse'); box on
subplot(4,1,2);surf(Sec.x, Sec.y, Sec.d,'EdgeColor','none'); caxis(c_axis); title('d Average of True field');view(2); axis tight; set(gca,'Ydir','reverse'); box on
subplot(4,1,3);surf(Prim.x, Prim.y, mean(zcs,3),'EdgeColor','none'); caxis(c_axis);view(2); axis tight; set(gca,'Ydir','reverse'); box on
subplot(4,1,4);surf(Prim.x, Prim.y, std(zcs,[],3),'EdgeColor','none'); caxis(c_axis); title('d Average of True field');view(2); axis tight; set(gca,'Ydir','reverse'); box on
% hold on; scatter(Sec.X(idx),Sec.Y(idx),'filled','r')


figure(11);clf; hold on;
for i_real=1:parm.n_real
   z=variogram_gridded_perso(zcs(:,:,i_real));
   plot(Prim.x,z,'b')
end
z=variogram_gridded_perso(Prim.d);
plot(Prim.x,z,'-r')
plot(Prim.x,1-covar.g(Prim.x*covar.cx(1)),'--k');
xlim([0 40])

figure(12); clf; hold on;

