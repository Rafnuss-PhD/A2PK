function elec=config_elec(elec)
%% CONFIGURATION compute the different configuration possible for a total
% number of electrode. There is no unit, you can think of it as an index
% for the matrix.
%
% INPUT:
% * method      defining the method of configuration 'wenner', 'pole-pole',
%               'dipole-dipole',...
% * elec.n      total number of electrode [no unit]
% * elec.depth_max   max depth. Above this depth, all configuration will be
%               removed
% * elec.config_max  Total number of configuration. If the code generate more
%               than elec.config_max, configuration will be remove by how close it
%               is from its neighbourhood configuration up to elec.config_max
%               configurations.
% * elec.plotit      Bolean (1,0) decide for ploting or not
%
% OUPUT:
% * data        in coelec.numn order:
%               1: +ve emmiteur electrode
%               2: -ve emmiteur electrode
%               3: +ve receiveur electrode
%               4: -ve receiveur electrode
% * pos         average position of the measure for an homogenous media.
%               This is compute in x as the middle of the electrodes and z
%               as the investigation depth
% * elec.depth_n    finaelec.n total number of configuration
% * k



assert(elec.n>4 && mod(elec.n,1)==0,'Invalid elec.n value')
assert(mod(elec.depth_max,1)==0,'Invalid elec.depth_max value')
assert(mod(elec.config_max,1)==0,'Invalid elec.config_max value')

assert(elec.n<=300,'R2 cannot have more than 300 unique electrodes site')
assert(elec.config_max<=6000,'R2 cannot have more than 6000 measurement')

%%
% * Defining for the seelec.nected method
%  A: +ve emmiteur electrode
%  B: -ve emmiteur electrode
%  M: +ve receiveur electrode
%  N: -ve receiveur electrode
switch elec.method
    case 'wenner'                   % A <--a--> M <--a--> N <--a--> B          k = 2 pi a
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a;
        N_fx    = @(u,n,a) u+a*2;
        B_fx    = @(u,n,a) u+a*3;
        good    = @(u,n,a,nn) u+(a+1)*3<=elec.n && nn<=1;
        ze_fx   = @(u,n,a) 0.173*a*3;
        xa_fx   = @(u,n,a) u+a*1.5;
        k_fx    = @(n,a) 2*pi*a;
    case 'pole-pole'                % A <--a--> M                               k = 2 pi a
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a;
        N_fx    = @(u,n,a) NaN;
        B_fx    = @(u,n,a) NaN;
        good    = @(u,n,a,nn) u+(a+1)<=elec.n && nn<=1;
        ze_fx   = @(u,n,a) 0.35*a;
        xa_fx   = @(u,n,a) u+a*0.5;
        k_fx    = @(n,a) 2*pi*a;
    case 'wenner-schelec.numberger' % A <--na--> M <--a--> N <--na--> B         k = pi n(n+1)a
        zec=[.173 .186 .189 .190];
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a*n;
        N_fx    = @(u,n,a) u+a*n+a;
        B_fx    = @(u,n,a) u+a*n+a+a*n;
        good    = @(u,n,a,nn) (u+(a+1)*n+(a+1)+(a+1)*n)<=nn;
        ze_fx   = @(u,n,a) zec(min(n,length(zec)))*(a*n+a+a*n);
        xa_fx   = @(u,n,a) u+n*a+a*0.5;
        k_fx    = @(n,a) pi*n*(n+1)*a;
    case 'dipole-dipole'            % B <--a--> A <--na--> M <--a--> N          k = pi n(n+1)(n+2)a
        zec=[.139 .174 .192 .203 .211 .216 .22 .224 .255];
        B_fx    = @(u,n,a) u;
        A_fx    = @(u,n,a) u+a;
        M_fx    = @(u,n,a) u+a+n*a;
        N_fx    = @(u,n,a) u+a+a*n+a;
        good    = @(u,n,a,nn) (u+(a+1)+(a+1)*n+(a+1))<=nn;
        ze_fx   = @(u,n,a) zec(min(n,length(zec)))*(a+a*n+a);
        xa_fx   = @(u,n,a) u+a+n*a*0.5;
        k_fx    = @(n,a)   pi*n*(n+1)*(n+2)*a;
    otherwise
        error('unvalid variabel method')
end


%%
% Computing all possibel configuration
config_n_pos = ceil((elec.n^3)/3/2);
data        = nan(config_n_pos,4);
pos         = nan(config_n_pos,2);
k           = nan(config_n_pos,1);
i=0;
u=0;
a=0;
good_u=true;
while good_u
    u=u+1;
    n=0;
    if ~good(u,n,a,elec.n) || ze_fx(u,n+1,a)> elec.depth_max
        good_u=false;
    else
        good_n=true;
    end
    while good_n
        n=n+1;
        a=0;
        if ~good(u,n,a,elec.n) || ze_fx(u,n,a+1)> elec.depth_max
            good_n=false;
        else
            good_a=true;
        end
        while good_a
            a=a+1;
            i=i+1;
            data(i,:)   = [A_fx(u,n,a) B_fx(u,n,a) M_fx(u,n,a) N_fx(u,n,a)];
            pos(i,:)    = [xa_fx(u,n,a) ze_fx(u,n,a)];
            k(i)        = k_fx(n,a);
            assert(any(data(i,:)<=elec.n))
            if ~good(u,n,a,elec.n) || ze_fx(u,n,a+1)> elec.depth_max
                good_a=false;
            end
        end
    end
end
data=data(1:i,:); pos=pos(1:i,:);k=k(1:i);
elec.config_n=size(data,1);

elec.pos_ini = pos;
elec.k_ini = k;

%%
% * Removing depth
idx=pos(:,2)<elec.depth_max;
if any(idx==0) && elec.plotit
    disp(['We removed data below the max depth ', num2str(elec.depth_max), 'm which correspond to ', num2str(elec.depth_n-sum(idx)),' point(s)' ])
end
data=data(idx,:); pos=pos(idx,:);k=k(idx);
elec.config_n=size(data,1);


%%
% * Removing duplicate
C = unique(pos,'rows');idx=[];
for i=1:size(C,1)
    u = find(pos(:,1)==C(i,1) & pos(:,2)==C(i,2));
    if numel(u)>1
        if k(u(1))<k(u(2))
            idx=[idx; u(1)];
        else
            idx=[idx; u(2)];
        end
    end
end
data(idx,:)=[];
pos(idx,:)=[];
k(idx)=[];

%%
% * Removing for very large dataset
% this look at the 5% of the max depth. It take the value of k which is
% greater than 10% of this data (which are locatated between 95% and 100% of maxdepth)
tmp=sort(k(0.95*elec.depth_max<pos(:,2)));
k_max=tmp(ceil(end/10));
idx=k>k_max;
data(idx,:)=[];
pos(idx,:)=[];
k(idx)=[];
elec.config_n=size(data,1);


%%
% * Datasample 10'000 points
if elec.config_n>20000
    [~,idx]=datasample(pos,10000,'Replace',false,'weight',max(0.000001,1/k));
    data=data(idx,:); pos=pos(idx,:);k=k(idx);
    elec.config_n=size(data,1);
end
%%
% * with different method....
% 1. k-mean, 2.iterative closest neighboor, 3. iterative avg closest
% neighboor, 4.voronoi, 5. random


if elec.config_max<elec.config_n
    switch elec.selection
        case 1 % k-mean
            idx = kmeans(pos,elec.config_max);
            data_n = nan(elec.config_max,4);
            pos_n  = nan(elec.config_max,2);
            k_n    = nan(elec.config_max);
            for i=1:elec.config_max
                u=find(idx==i);
                [~, idx2]=min(k(u));
                data_n(i,:)=data(u(idx2),:); pos_n(i,:)=pos(u(idx2),:); k_n(i)=k(u(idx2));
            end
            data=data_n;
            pos=pos_n;
            k=k_n;
            
        case 2 % iterative removal of the closest neighboohood
            D = pdist2(pos,pos);
            D(D==0)=Inf;
            while size(D,1)>elec.config_max
                [D_min,idx]=min(D);
                [~,idx2]=min(D_min);
                idx1=idx(idx2);
                if k(idx1)<k(idx2)
                    idx=idx1;
                else
                    idx=idx2;
                end
                D(idx,:)=[];
                D(:,idx)=[];
                data(idx,:)=[];
                pos(idx,:)=[];
                k(idx)=[];
            end
        case 3 % iterative removal of the averaged closest point
            D = pdist2(pos,pos);
            while size(D,1)>elec.config_max
                D_avg=mean(D);
                [~,idx]=min(D_avg);
                D(idx,:)=[];
                D(:,idx)=[];
                data(idx,:)=[];
                pos(idx,:)=[];
                k(idx)=[];
            end
        case 4 % voronoi
            while size(pos,1)>elec.config_max
                [v , c] = voronoin(pos); 
                tess_area=zeros(size(c,1),1);
                for i = 1:size(c,1)
                    ind = c{i}';
                    tess_area(i,1) = polyarea( v(ind,1) , v(ind,2) );
                end
                [~,idx]=min(tess_area);
                if elec.plotit
                    clf;hold on; voronoi(pos(:,1),pos(:,2));
                    scatter(pos(idx,1),pos(idx,2));drawnow
                end
                data(idx,:)=[];
                pos(idx,:)=[];
                k(idx)=[];
            end
        case 5 % random
            D=pdist2(pos,pos);
            tic
            dist_best=0;
            dist_c =0;
            %i=1;
            while toc<5;
                [~,idx]=datasample(pos,elec.config_max,'Replace',false);
                for u=1:elec.config_max
                    for uu= u:elec.config_max
                        dist_c = dist_c + D(u,uu);
                    end
                end
                if  dist_c > dist_best
                    idx_best = idx;
                    dist_best = dist_c;
                end
                dist_c =0;
%                 be(i)=dist_best;
%                 i=i+1;
            end
            data=data(idx_best,:);
            pos=pos(idx_best,:);
            k=k(idx_best);
            
    end
    elec.config_n=elec.config_max;
end


%% SORT
%
[elec.data,idx]=sortrows(data);
elec.pos=pos(idx,:);
elec.k=k(idx);


elec.pseudo_x       = elec.pos(:,1)*(elec.x(2)-elec.x(1));
elec.pseudo_y       = elec.pos(:,2)*(elec.x(2)-elec.x(1));
end