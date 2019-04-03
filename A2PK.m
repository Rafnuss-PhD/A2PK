%% Area-to-Point kriging A2PK
% *Area-to-Point kriging A2PK* generates stocastic Gaussian realization
% |*z*| constrained to a variable |*Z*| which is lineraly related
%  to |*z*| by |*G*|:
%
% $$\mathbf{Z = Gz}$$
%
% 
% The argument of the function are
% * |x|: vector of coordinates along the first axes
% * |y|: vector of coordinates along the second axes
% * |hd|: Hard data, see |functions/sampling_pt.m|
% * |Z|: Coarse scale 
% * |G|: Linear link
% * |covar|: covariance structure defined as in
% |FastGaussianSimulation/covarIni.m|
% * |n_real|: number of realization desired
% 
% The output of the function are
% * |zcs|: Conditional realizations cells of size n_real
% * |zh|: Krigging estimation
% * |S|: Variance of kriging estimation 
%
% *Script*
% 
% *Exemples*: Available in the folder |examples| with unconditional and
% conditional Gaussian simulation and a case study of electrical tomography


function [zcs,zh,S] = A2PK(x,y,hd,Z,G,covar,n_real)

%% Checkin input argument
validateattributes(x,{'numeric'},{'vector'})
validateattributes(y,{'numeric'},{'vector'})
if isempty(hd)
    hd.id=[];
    hd.n=0;
    hd.d=[];
end
validateattributes(hd,{'struct'},{})
% validateattributes(hd.id,{'numeric'},{'vector','integer'})
% validateattributes(hd.d,{'numeric'},{'vector'})
validateattributes(hd.n,{'numeric'},{'integer','nonnegative','scalar'})
validateattributes(Z,{'numeric'},{'2d'})
validateattributes(G,{'numeric'},{'2d'})
validateattributes(n_real,{'numeric'},{'integer','positive','scalar'})

assert(all(size(G)==[numel(Z) numel(x)*numel(y)]),'The size of the upscaling matrix G need to be equal to [numel(Z) numel(x)*numel(y)]')

%% Inlude the dependancy
addpath('./FastGaussianSimulation');
addpath('./functions');

%% Convert the covariance
covar = covarIni(covar);
[X, Y] = meshgrid(x, y); X=X(:); Y=Y(:);
nx=numel(x); ny=numel(y); nZ=numel(Z);

%% Calcul of the covariance and cross covariance.
% Because the covariance matrix can be quite large, an iterative loop might
% be faster
if nx*ny>100000
    Czz=zeros(nx*ny,nx*ny);
    wradius = parm.k.wradius;
    for ixy=1:nx*ny
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


%% Compute the kriging weights W, variance S and kriging map zh

W = (CCa \ CCb)';

zh = reshape( W * [Z(:) ; hd.d], ny, nx);
S =  reshape(covar.g(0) - diag(W * CCb), ny, nx);

%% Create Simulation
zcs=nan(ny, nx,n_real);
for i_real=1:n_real
    zs = FGS(struct('s',[numel(x) numel(y)]), covar);
    zhs = reshape( W * [G * zs{1}(:) ; zs{1}(hd.id)'], ny, nx);
    zcs(:,:,i_real) = zh + (zs{1} - zhs);
end
