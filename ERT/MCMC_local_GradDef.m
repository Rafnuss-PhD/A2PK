
%% Load data
fieldname = 'GEN-AsSumbitted_2018-07-09_16-34';
addpath('../functions','R2');
load(['result/' fieldname]);
load(['result/' fieldname '_cond']);
%addpath('C:\Users\rnussba1\Documents\MATLAB\Colormaps\')

%% Set-up forward simulator
f=gen.f;

f.content = fileread( ['result/' fieldname '_IO-file/forward/R2.in']);
F = griddedInterpolant({Sec.y,Sec.x},1000./Nscore.inverse(Sec.d),'linear');
f.rho = F({f.grid.y, f.grid.x});

f.filepath          = 'data_gen/IO-file-prop/';
mkdir(f.filepath)
copyfile('R2/R2.exe',[f.filepath 'R2.exe'])
copyfile(['result/' fieldname '_IO-file/forward/electrodes.dat'],[f.filepath 'electrodes.dat'])

gen_f_output_resistancewitherror=gen.f.output.resistancewitherror;
gen_i_b_wgt = gen.i.b_wgt;


%% Conditioning
covar = kriginginitiaite(gen.covar);
Czzh = covar.g(pdist2([Prim.X(Prim_pt.id ) Prim.Y(Prim_pt.id)]*covar.cx,[Prim.X(:) Prim.Y(:)]*covar.cx));
Czh = Czzh(:,Prim_pt.id);

W=zeros(numel(Prim.x)*numel(Prim.y),Prim_pt.n);
Czhinv = inv(Czh);
for ij=1:numel(Prim.x)*numel(Prim.y)
     W(ij,:) = Czhinv * Czzh(:,ij);
end

zh = reshape( W * Prim_pt.d, numel(Prim.y), numel(Prim.x));


%% Simple parrallelization
theta= 0.05; %0.05 -> 0.002
n = 6*60*24*3; %chain length
ns = 1000;
np=10;

covar = kriginginitiaite(gen.covar);

accrate=false(1,n);
l=nan(1,n);

i_s=0;
m = nan(numel(Prim.y),numel(Prim.x),ceil(n/ns));
m_pres = fftma_perso(covar, struct('x',Prim.x,'y',Prim.y));
f.rho(f.grid.inside) = 1000./Nscore.inverse(m_pres);
fsim_resistance = Matlat2R2min(f);
l_pres = sum (((fsim_resistance - gen_f_output_resistancewitherror) ./ ( gen_i_b_wgt*gen_f_output_resistancewitherror)).^2 );


i_prop=0;


f.filepath = 'data_gen/IO-file/';
mkdir(f.filepath)
copyfile('R2/R2.exe',[f.filepath 'R2.exe'])
copyfile(['result/' fieldname '_IO-file/forward/electrodes.dat'],[f.filepath 'electrodes.dat'])
copyfile(['result\' fieldname '_IO-file\forward\protocol.dat'],[f.filepath 'protocol.dat'])


f_tmp=f;

while i_prop<n-10
    
    % New proposal
    i_prop = i_prop+1;
    
    m_prop = m_pres*cos(theta) + fftma_perso(covar, struct('x',Prim.x,'y',Prim.y))*sin(theta);

    % Compute proposal
    f_tmp.rho(f_tmp.grid.inside) = 1000./Nscore.inverse(m_prop_par(:,:,i_prop_par));
    fsim_resistance = Matlat2R2min(f_tmp);
    l(i_prop) = sum (((fsim_resistance - gen_f_output_resistancewitherror) ./ ( gen_i_b_wgt*gen_f_output_resistancewitherror)).^2 );


    % Update the step \theta
    if mod(i_prop,20)==0 &&  sqrt(l_pres/numel(gen.f.output.resistancewitherror))>2
        dtheta = 0.001*(mean(accrate(i_prop-19:i_prop))-.3)*10;
        disp(['Change of theta: ' num2str(theta) '->' num2str(theta + dtheta) ' | '  num2str(100*mean(accrate(i_prop-19:i_prop))) '%)'])
        if theta + dtheta>0 && theta + dtheta<1
            theta = theta + dtheta;
        end
    end

    % Save the present model if  not burning zone.
    if mod(i_prop,ns)==0 && sqrt(l_pres/numel(gen.f.output.resistancewitherror))>2
        i_s=i_s+1;
        m(:,:,i_s)=m_pres;
    end

    % Test for acceptance.
    if min(exp(-1/2*(l(i_prop)-l_pres)),1) > rand()
        accrate(i_prop)=true;
        m_pres = m_prop_par(:,:,i_prop_par);
        l_pres = l(i_prop);
        
        if l_pres<l_best
            l_best=l_pres;
            m_best=m_pres;
        end
        
        disp([num2str(i_prop) '| Accepted (' num2str(100*mean(accrate(2:i_prop))) '%)'])
        
    else
        disp([num2str(i_prop) '| Rejected (' num2str(100*mean(accrate(2:i_prop))) '%)'])
    end
end



%% Complex parrallelization
% Initiate Monte-carlo param
n = 6*60*24*9*3; %chain length

numcores = feature('numcores');
numcores=10;

i_prop=0;

covar = kriginginitiaite(gen.covar);

accrate=false(1,n);
l=nan(1,n);

nt=20; % change theta  every nt proposition
i_t=1;
theta=nan(1,round(n/nt));
theta(1)= 0.05; %0.05 -> 0.002

ns = 1000; % save every ns proposition
i_s=0;
s_prop=nan(1,ceil(n/ns));
m = nan(numel(Prim.y),numel(Prim.x),ceil(n/ns));
m_pres = fftma_perso(covar, struct('x',Prim.x,'y',Prim.y));
f.rho(f.grid.inside) = 1000./Nscore.inverse(m_pres);
fsim_resistance = Matlat2R2min(f);
l_pres = sum (((fsim_resistance - gen_f_output_resistancewitherror) ./ ( gen_i_b_wgt*gen_f_output_resistancewitherror)).^2 );
m_best = nan(numel(Prim.y),numel(Prim.x));
l_best = l_pres;

l_true = sum (((gen.f.output.resistance - gen.f.output.resistancewitherror) ./ (gen.i.a_wgt + gen.i.b_wgt*gen.f.output.resistancewitherror)).^2 );

accrate=[accrate false(1,n-n_old)];
l=[l nan(1,n-n_old)];
s_prop=[s_prop nan(1,ceil((n-n_old)/ns))];
m_old=m;
m = nan(numel(Prim.y),numel(Prim.x),ceil(n/ns));
m(:,:,1:size(m_old,3))=m_old;


% parrallele structure
nt=41; % size of the tree. need to be odd
acc_rate=.3; % Target acceptance rate. 

tmp=[0 repelem(1:nt,1,2)];
adjacency=zeros(nt);
for i=2:nt
    adjacency(tmp(i),i) = 1;
end
G = digraph(adjacency);%,cellstr(num2str((1:tree.size)')));
G.Nodes.acc=[NaN repmat([true false],1,(nt-1)/2)]';
G.Nodes.prob = ones(nt,1);
for i=2:2:nt-1
    G.Nodes.prob(i)=acc_rate*G.Nodes.prob(predecessors(G,i));
    G.Nodes.prob(i+1)=(1-acc_rate)*G.Nodes.prob(predecessors(G,i));
end

figure(1);plot(G,'NodeLabel',G.Nodes.acc)
figure(1);plot(G,'NodeLabel',G.Nodes.prob)

[~,i_max] = maxk(G.Nodes.prob,numcores+1);
H = subgraph(G,i_max);
figure(2); plot(H,'NodeLabel',H.Nodes.acc)

% Initiate model for the parallele structure

for i_prop_par=1:numcores+1
    f.filepath = ['data_gen/IO-file-' num2str(i_prop_par) '/'];
    mkdir(f.filepath)
    copyfile('R2/R2.exe',[f.filepath 'R2.exe'])
    delete([f.filepath 'R2.in'])
    copyfile(['result/' fieldname '_IO-file/forward/R2.in'],[f.filepath 'R2.in'])
    copyfile(['result/' fieldname '_IO-file/forward/electrodes.dat'],[f.filepath 'electrodes.dat'])
    copyfile(['result\' fieldname '_IO-file\forward\protocol.dat'],[f.filepath 'protocol.dat'])
end


tic
while i_prop<n-10
    
    m_prop_par=nan(numel(Prim.y),numel(Prim.x),numcores+1);
    m_prop_par(:,:,1) = m_pres; %*cos(theta) + fftma_perso(covar, struct('x',Prim.x,'y',Prim.y))*sin(theta);
    for i=2:numcores+1
        
        X=randn();
        L=bwdist(padarray(1,[250,250]);
        
        filter2(X,L));
        
        
        
        m_prop_par(:,:,i) = m_prop_par(:,:,predecessors(H,i))*cos(theta(i_t)) + fftma_perso(covar, struct('x',Prim.x,'y',Prim.y))*sin(theta(i_t));
        
    end
    
    l_par = nan(numcores+1,1);
    parfor i_prop_par = 2:numcores+1
        f_tmp=f;
        f_tmp.filepath = ['data_gen/IO-file-' num2str(i_prop_par) '/'];
        m_prop_par_uncond = m_prop_par(:,:,i_prop_par);
        m_prop_par_cond = m_prop_par_uncond + zh - reshape( W * m_prop_par_uncond(Prim_pt.id), numel(Prim.y), numel(Prim.x));
        f_tmp.rho(f_tmp.grid.inside) = 1000./Nscore.inverse(m_prop_par_cond);
        fsim_resistance = Matlat2R2min(f_tmp);
        l_par(i_prop_par) = sum (((fsim_resistance - gen_f_output_resistancewitherror) ./ ( gen_i_b_wgt*gen_f_output_resistancewitherror)).^2 );
    end

    assert(all(~isnan(l_par(2:end))))
    
    i_prop_par=2;
    while ~isempty(i_prop_par)
        
        % New proposal
        i_prop = i_prop+1;
        
        % Update the step \theta. every 20 propositions and before burning
        % zone (WRMSE>2)
        if mod(i_prop,nt)==0 &&  sqrt(l_pres/numel(gen.f.output.resistancewitherror))>2
            % last change of thetat i_prop=8692
            % change at i_prop=9352 to theta=0.003
            % change at i_prop=49467 to theta=0.002
            % change at i_prop=68285 to theta=0.004
            % change at i_prop=72106 to theta=0.003
            dtheta = 0.001*(mean(accrate(i_prop-19:i_prop))-.3)*10;
            i_t=i_t+1;
            disp(['Change of theta: ' num2str(theta(i_t-1)) '->' num2str(theta(i_t-1) + dtheta) ' | '  num2str(100*mean(accrate(i_prop-19:i_prop))) '%)'])
            theta(i_t) = theta(i_t-1) + dtheta;
            assert(~isnan(theta(i_t)),'theta is nan')
        end
        
        % Save the present model if not burning zone.
        if mod(i_prop,ns)==0 && sqrt(l_pres/numel(gen.f.output.resistancewitherror))<6
            i_s=i_s+1;
            s_prop(i_s)=i_prop;
            m(:,:,i_s)=m_pres;
        end
        
        % Save the sum of square error
        l(i_prop)=l_par(i_prop_par);
        

        % Test for acceptance.
        succ = successors(H,i_prop_par);
        if min(exp(-1/2*(l(i_prop)-l_pres)),1) > rand()
            accrate(i_prop)=true;
            m_pres = m_prop_par(:,:,i_prop_par);
            l_pres = l(i_prop);
            disp([num2str(i_prop) '| Accepted (' num2str(100*mean(accrate(2:i_prop))) '%)'])
            i_prop_par = succ(H.Nodes.acc(succ)==1);
            
            if l_pres<l_best
                l_best=l_pres;
                m_best=m_pres;
            end
        else
            disp([num2str(i_prop) '| Rejected (' num2str(100*mean(accrate(2:i_prop))) '%)'])
            
            i_prop_par = succ(H.Nodes.acc(succ)==0);
        end
    end
    toc
end

% save('result/MCMC','accrate','l','m','n','l_best','m_best','m_pres','l_pres','s_prop','theta','-v7.3');

%% Show result
figure(1); plot(1:nt:nt*numel(theta),theta)

figure(2); clf;  
subplot(2,1,1); imagesc(Prim.d); caxis([-3 3])
subplot(2,1,2); imagesc(m_pres); caxis([-3 3])
subplot(2,1,2); imagesc(m(:,:,i_s)); caxis([-3 3])

figure(25); clf;  
subplot(5,1,1); imagesc(Prim.d); caxis([-3 3])
subplot(5,1,2); imagesc(mean(m(:,:,100:10:i_s),3)); caxis([-3 3])
subplot(5,1,3); imagesc(var(m(:,:,100:10:i_s),[],3));
subplot(5,1,4); imagesc(m(:,:,i_s)); caxis([-3 3])
subplot(5,1,5); imagesc(m_best); caxis([-3 3])

% figure; plot(l(accrate))
figure(3); clf; hold on;
plot(sqrt(l/numel(gen.f.output.resistancewitherror))); set(gca,'YScale','log');% set(gca,'XScale','log')
plot([0 i_prop],sqrt(l_true/numel(gen.f.output.resistancewitherror)).*[1 1],'-r')

post = find(accrate);

maccmean = mean(m(:,:,accrate),3);
maccstd = std(m(:,:,accrate),[],3);

figure; 
subplot(4,1,1); imagesc(Prim.d); caxis([-3 3])
subplot(4,1,2); imagesc(maccmean); caxis([-3 3])

subplot(4,1,3); imagesc(maccstd); caxis([-3 3])


figure(5); clf;   colormap(viridis())
subplot(7,1,1);surf(Prim.x,Prim.y,Prim.d,'EdgeColor','none','facecolor','flat'); caxis([-3 3])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); %colorbar('southoutside'); title('Kriging Estimate')
subplot(7,1,2);surf(Prim.x,Prim.y,zh,'EdgeColor','none','facecolor','flat'); caxis([-3 3])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); %colorbar('southoutside'); title('Kriging Estimate')
subplot(7,1,3);surf(Prim.x,Prim.y,S,'EdgeColor','none','facecolor','flat'); caxis([0 1])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y');% colorbar('southoutside'); title('Kriging Estimate')
subplot(7,1,4);surf(Prim.x,Prim.y,mean(zcs,3),'EdgeColor','none','facecolor','flat'); caxis([-3 3])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); %colorbar('southoutside'); title('Kriging Estimate')
subplot(7,1,5);surf(Prim.x,Prim.y,std(zcs,[],3),'EdgeColor','none','facecolor','flat'); caxis([0 1])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y'); %colorbar('southoutside'); title('Kriging Estimate')
subplot(7,1,6);surf(Prim.x,Prim.y, maccmean,'EdgeColor','none','facecolor','flat'); caxis([0 1])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y');% colorbar('southoutside'); title('Kriging Estimate')
subplot(7,1,7);surf(Prim.x,Prim.y, maccstd,'EdgeColor','none','facecolor','flat'); caxis([0 1])
view(2); axis equal tight; set(gca,'Ydir','reverse'); xlabel('x');ylabel('y');% colorbar('southoutside'); title('Kriging Estimate')
