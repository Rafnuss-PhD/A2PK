clear all; addpath(genpath('./.')); dbstop if error 
load('result-A2PK/GEN-test_2017-09-26_17-22.mat'); 

covar=gen.covar;

x=grid_gen.x; 
y = grid_gen.y;
[X, Y] = meshgrid(x, y);
X=X(:); Y=Y(:);
nx=numel(x); ny=numel(y);

Cz=zeros(nx*ny,nx*ny);
parpool(20)
parfor ixy=1:nx*ny
    u=zeros(1,nx*ny)
    id = X-X(ixy)<covar.range(1)*parm.k.wradius  & Y-Y(ixy)<covar.range(2)*parm.k.wradius;
    u(id) = covar.g(pdist2([X(id) Y(id)]*covar.cx,[X(ixy) Y(ixy)]*covar.cx));
    Cz(ixy,:) = u;
end

Cz=sparse(Cz);
save('Cz','Cz');

