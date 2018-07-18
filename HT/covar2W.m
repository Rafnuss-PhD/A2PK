function W = covar2W(covar,Prim,Sec,G,Prim_pt,Cmt)

% Compute the covariance of the spatial model
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
% parpool(48)
parfor ij=1:numel(Prim.x)*numel(Prim.y)
     W(ij,:) = CCa \ CCb(:,ij);
end

end