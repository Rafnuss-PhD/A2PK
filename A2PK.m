%% Area-to-Point kriging A2PK
% *Area-to-Point kriging A2PK* generates stocastic Gaussian realization
% |*z*| constrained to a variable |*Z*| which is lineraly related
%  to |*z*| by |*G*|
%
% $$\mathbf{Z = Gz}$$
%
% 
% The argument of the function are
% * x
% * y
% * hd
% * Z
% * G
% * covar: covariance structure defined as

% * n_real: number of realization desired
% 


function [zcs] = A2PK(x,y,hd,Z,G,covar,n_real)

%% Checkin input argument
%TO BE DONE

%% Inlude the dependancy
addpath('../functions');

%% Convert the covariance
covar = kriginginitiaite(covar);
[X, Y] = meshgrid(x, y); X=X(:);Y=Y(:);
nx=numel(x); ny=numel(y); nZ=numel(Z);

%% Calcul of the covariance and cross covariance.
% Because the covariance matrix can be quite large, an iterative loop might
% be faster
if nx*ny>100000
    Czz=zeros(nx*ny,nx*ny);
    wradius = parm.k.wradius;
    parfor ixy=1:nx*ny
        u=zeros(1,nx*ny);
        id = X-X(ixy)<covar.range(1)*wradius  & Y-Y(ixy)<covar.range(2)*wradius;
        u(id) = covar.g(pdist2([X(id) Y(id)]*covar.cx,[X(ixy) Y(ixy)]*covar.cx));
        Czz(ixy,:) = u;
    end
    Czz=sparse(Czz);
else
    Czz = covar.g(squareform(pdist([X Y]*covar.cx)));
end

CzZ = Czz * G';
CZZ = G * CzZ;

Chd = Czz(hd.id,hd.id);
Czhd = Czz(:,hd.id);
CZhd = CzZ(hd.id,:);

CCa = [ CZZ CZhd' ; CZhd Chd ];
CCb = [ CzZ' ; Czhd' ];


%% Compute the kriging weights and kriging map
% note that this can be very long... 
W=zeros(nx*ny,nZ+hd.n);
parfor ij=1:nx*ny
    W(ij,:) = CCa \ CCb(:,ij);
end
zh = reshape( W * [Z(:) ; hd.d], ny, nx);


%% Create Simulation
zcs=nan(ny, nx,n_real);
parfor i_real=1:n_real
    zs = fftma_perso(covar, struct('x',x,'y',y));
    zhs = reshape( W * [G * zs(:) ; hd.d], ny, nx);
    zcs(:,:,i_real) = zh + (zs - zhs);
end
