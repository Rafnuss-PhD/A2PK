function [zcs] = ESA2P(y,x,prim,D,G,parm,Cmt)
% Kriging Error Simulation of Areal-to-point
% Paramter settings
if ~isfield(parm, 'seed_U'),        parm.seed_U         = 'shuffle'; end
if ~isfield(parm, 'n_real'),        parm.n_real         = 1; end
if ~isfield(parm.k, 'wradius'),        parm.k.wradius         = 2; end

covar = kriginginitiaite(parm.k.covar);
[X, Y] = meshgrid(x, y); X=X(:);Y=Y(:);
nx=numel(x); ny=numel(y); nd=numel(D);


if nx*ny>100000
    Cz=zeros(nx*ny,nx*ny);
    wradius = parm.k.wradius;
    parfor ixy=1:nx*ny
        u=zeros(1,nx*ny);
        id = X-X(ixy)<covar.range(1)*wradius  & Y-Y(ixy)<covar.range(2)*wradius;
        u(id) = covar.g(pdist2([X(id) Y(id)]*covar.cx,[X(ixy) Y(ixy)]*covar.cx));
        Cz(ixy,:) = u;
    end
    Cz=sparse(Cz);
else
    Cz = covar.g(squareform(pdist([X Y]*covar.cx)));
end
disp('Cz computed')
Czd = Cz * G';
Cd = G * Czd;
Cd = Cd;

Czh = covar.g(squareform(pdist([prim.y(:) prim.x(:)]*covar.cx)));
Czzh = covar.g(pdist2([Y X]*covar.cx,[prim.y(:) prim.x(:)]*covar.cx));

Czhd = Czd( prim.id ,:);

CCa = [ Cd Czhd' ; Czhd Czh ];
CCb = [ Czd' ; Czzh' ];

disp('Covariances computed')

W=zeros(nx*ny,nd+prim.n);
parfor ij=1:nx*ny
     W(ij,:) = CCa \ CCb(:,ij);
end
zh = reshape( W * [D(:) ; prim.d], ny, nx);

disp('W computed')

rng(parm.seed_U);
zcs=nan(ny, nx,parm.n_real);
parfor i_real=1:parm.n_real
    zs = fftma_perso(covar, struct('x',x,'y',y));
    zhs = reshape( W * [G * zs(:) ; prim.d], ny, nx);
    zcs(:,:,i_real) = zh + (zs - zhs);
end

%dcs = reshape(G * zcs(:), numel(Sec.y), numel(Sec.x));
