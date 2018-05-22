
%% Load data
fieldname = 'GEN-600x40_2017-12-21_15-44';
addpath('../functions','R2');
load(['result/' fieldname]);
addpath('C:\Users\rnussba1\Documents\MATLAB\Colormaps\')

%% Set-up fprward simulator
f={};
f.res_matrix        = gen.Rho.f.res_matrix;
f.grid              = gen.Rho.f.grid;    
f.header            = 'Forward';  % title of up to 80 characters
f.job_type          = 0;
f.filepath          = 'data_gen/IO-file-prop/';
f.readonly          = 0;
f.alpha_aniso       = gen.Rho.f.alpha_aniso;
f.elec_spacing      = gen.Rho.f.elec_spacing;
f.elec_id           = gen.Rho.f.elec_id;
f.num_regions       = gen.Rho.f.num_regions;
f.rho_min           = gen.Rho.f.rho_min;
f.rho_avg           = gen.Rho.f.rho_avg;
f.rho_max           = gen.Rho.f.rho_max;
mkdir(f.filepath)

%%
theta= 0.05; %0.05 -> 0.002
n = 6*60*24*3; %chain length
ns = 1000;

covar = kriginginitiaite(gen.covar);

accrate=false(1,n);
l=nan(1,n);

i_s=0;
m = nan(numel(Prim.y),numel(Prim.x),ceil(n/ns));
m_pres = fftma_perso(covar, struct('x',Prim.x,'y',Prim.y));
f.rho = 1000./Nscore.inverse(m_pres);
f = Matlat2R2(f,gen.Rho.elec);
l_pres = sum (((f.output.resistance - gen.Rho.f.output.resistancewitherror) ./ (gen.Rho.i.a_wgt + gen.Rho.i.b_wgt*gen.Rho.f.output.resistancewitherror)).^2 );


l=[l false(1,5000)];
m = cat(3, m, nan(numel(Prim.y),numel(Prim.x),5000));
n=n+5000;

for i_prop=1:n
    
    m_prop = m_pres*cos(theta) + fftma_perso(covar, struct('x',Prim.x,'y',Prim.y))*sin(theta);
    
    % sim
    f.rho = 1000./Nscore.inverse(m_prop);
    f = Matlat2R2(f,gen.Rho.elec);
    l(i_prop) = sum (((f.output.resistance - gen.Rho.f.output.resistancewitherror) ./ (gen.Rho.i.a_wgt + gen.Rho.i.b_wgt*gen.Rho.f.output.resistancewitherror)).^2 );

    if min(exp(-1/2*(l(i_prop)-l_pres)),1) > rand()
        accrate(i_prop)=true;
        m_pres = m_prop;
        l_pres = l(i_prop);
        %plot(i_prop,l(i_prop),'or')
    else
        %plot(i_prop,l(i_prop),'ok')
    end
    
%     if mod(i_prop,10)==0
%         dtheta = 0.001*(mean(accrate(i_prop-9:i_prop))-.3)*10;
%         disp(['Change of theta: ' num2str(theta) '->' num2str(theta + dtheta) ' | '  num2str(100*mean(accrate(i_prop-9:i_prop))) '%)'])
%         if theta + dtheta>0 && theta + dtheta<1
%             theta = theta + dtheta;
%         end
%     end

    if mod(i_prop,ns)==0
        i_s=i_s+1;
        m(:,:,i_s)=m_pres;
    end

    disp([num2str(i_prop) '| Accepted: ' num2str(accrate(i_prop)) ' (' num2str(100*mean(accrate(2:i_prop))) '%)'])
end



figure; plot(l(accrate))

figure; 
subplot(2,1,1); imagesc(Prim.d); caxis([-3 3])
subplot(2,1,2); imagesc(m(:,:,i_prop)); caxis([-3 3])

save('result/MCMC','accrate','l','m','n','-v7.3');

%% Show result

l_true = sum (((gen.Rho.f.output.resistance - gen.Rho.f.output.resistancewitherror) ./ (gen.Rho.i.a_wgt + gen.Rho.i.b_wgt*gen.Rho.f.output.resistancewitherror)).^2 );

imagesc(m(:,:,i_pres))

figure(3); clf; hold on;
plot(l); set(gca,'YScale','log')
plot([0 numel(accrate)],[l_true l_true],'-r')

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
