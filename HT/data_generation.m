function [filename, grid_gen, K_true, K, gen] = data_generation(gen)
% DATA_GENERATION is basically creating all the data required for a simulation.
% INPUT:
% GRID
%   * xmax:         length of the grid_gen of x [unit]
%   * ymax:         length of the grid_gen of y [unit]
%   * sx:          scale level. number of cell is 2^gen.scale.x(i)+1
%   * sy:          scale level. number of cell is 2^gen.scale.y(i)+1
% STRUCTURE
%   * method:       method of generation: 'fromRho', ''
%   * samp:         Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
%   * samp_n:       Number of well or number of point
%   * covar:
%       * modele    covariance structure
%       * c
%   * mu            parameter of the first field.
%   * std
% RHO R2
%   * Rho.grid_gen.nx           = 300;
%   * Rho.grid_gen.ny           = 60; % log-spaced grid_gen.
%   * Rho.elec.spacing      = 2; % in grid_gen spacing unit.
%   * Rho.elec.config_max   = 6000; % number of configuration of electrode maximal
%   * Rho.method          = 'R2';
% OTHER
%   * plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
%   * saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
%   * name                = 'Small_range';
%   * seed                = 123456;
% OUTPUT:
%       - K_true:      	Hydraulic conductivity true field, matrix (grid_gen.nx x grid_gen.ny) (data or generated)
%       - rho_true:      	Electrical conductivity true field, matrix (grid_gen.nx x grid_gen.ny)  (data or from K_true)
%       - K:          	Hydraulic conductivity at some point, structure: location (K.x, K.y) of data (K.d) (sampled from K_true)
%       - g:          	Electrical conductivity at some point, structure: location (g.x, g.y) of data (g.d) (sampled from rho_true)
%       - G:          	Electrical conductivity measured grid_gen, matrix (G.nx x G.ny) of data (G.d) (ERT inverse)
%
% Author: Raphael Nussbaumer

%% * *INPUT CEHCKING*
assert(isfield(gen, 'xmax'))
assert(isfield(gen, 'ymax'))
assert(isfield(gen, 'nx'))
assert(isfield(gen, 'ny'))
if ~isfield(gen, 'method'); gen.method = 'Random'; end
if ~isfield(gen, 'samp'); gen.samp = 2; end
if ~isfield(gen, 'samp_n'); gen.samp_n = round(1/100 * gen.nx*gen.ny); end
if ~isfield(gen, 'mu'); gen.mu = 0; end
if ~isfield(gen, 'std'); gen.std = 1; end
% Other
if ~isfield(gen, 'plotit'); gen.plotit = 0; end
if ~isfield(gen, 'saveit'); gen.saveit = 1; end
if ~isfield(gen, 'forwardonly'); gen.saveit = 0; end
if ~isfield(gen, 'name'); gen.name = ''; end
if ~isfield(gen, 'seed'); gen.seed = 'default'; end
if ~isfield(gen, 'plot'); gen.plot = true; end
tic

rng(gen.seed)

%% * 2. *construction of the grid_gen*
grid_gen.nx     = gen.nx;
grid_gen.ny     = gen.ny;
grid_gen.nxy    = grid_gen.nx*grid_gen.ny; % total number of cells

grid_gen.dx=gen.xmax/(grid_gen.nx-1);
grid_gen.dy=gen.ymax/(grid_gen.ny-1);

grid_gen.x=linspace(0, gen.xmax, grid_gen.nx); % coordinate of cells center
grid_gen.y=linspace(0, gen.ymax, grid_gen.ny);
grid_gen.xy=1:grid_gen.nxy;

[grid_gen.X, grid_gen.Y] = meshgrid(grid_gen.x,grid_gen.y); % matrix coordinate




%% * 3. *Generate field*
f_Heinz         = @(phi,a,b) 10.^(a *phi - b);

phi_true        = gen.mu + gen.std*fftma_perso(gen.covar, grid_gen);
K_true          = f_Heinz(phi_true,6.66,4.97);
R_true          = 1./K_true;


%% Sampling
K = sampling_pt(grid_gen,K_true,gen.samp,gen.samp_n);

% Plot
if gen.plotit
    K_true_t = (log(K_true) - mean(log(K_true(:)))) ./ std(log(K_true(:)));
    K_dt = (log(K.d) - mean(log(K_true(:)))) ./ std(log(K_true(:)));
    
    
    figure(1);clf; subplot(2,1,1); hold on; title('Electrical Conductivity [mS/m]');xlabel('x [m]'); ylabel('y [m]')
    imagesc(grid_gen.x,grid_gen.y,K_true_t);colorbar;  scatter(K.x,K.y,K.d); legend({'Sampled location'}); axis tight;
    subplot(2,1,2); hold on; title('Histogram'); xlabel('Electrical Conductivity [mS/m]');
    ksdensity(K_true_t(:)); ksdensity(K_dt(:)); legend({'True','Sampled'}); box on;
    
    [gamma_x, gamma_y] = variogram_gridded_perso(K_true_t);
    figure(2); clf; subplot(2,1,1); hold on; title('Horizontal (x) Variogram')
    plot(grid_gen.x(1:end/2),gamma_x(1:end/2)./std(K_true_t(:))^2);
    % plot([gen.covar.modele(1,2) gen.covar.modele(1,2)],[0 1])
    plot(grid_gen.x(1:end/2),1-gen.covar.g(grid_gen.x(1:end/2)/gen.covar.range(2)),'linewidth',2)
    subplot(2,1,2); hold on; title('Vertical (y) Variogram')
    plot(grid_gen.y(1:end/2),gamma_y(1:end/2)./std(K_true_t(:))^2);
    % plot([gen.covar.modele(1,3) gen.covar.modele(1,3)],[0 1])
    plot(grid_gen.y(1:end/2),1-gen.covar.g(grid_gen.y(1:end/2)/gen.covar.range(1)),'linewidth',2)
    
    keyboard;close all; % try different initial data if wanted
end

%%
filepath       = 'data_gen/IO-file/';
delete([filepath '*']);

% Forward Grid
f                   = gen.Rho.f;
f.grid.x            = grid_gen.x;
f.grid.y            = grid_gen.y;
f.grid.nx           = grid_gen.nx;
f.grid.ny           = grid_gen.ny;
cell2vertex         = @(x) [x(1)-(x(2)-x(1))/2  x(1:end-1)+diff(x)/2  x(end)+(x(end)-x(end-1))/2];
f.grid.x_n          = cell2vertex(f.grid.x);
f.grid.y_n          = cell2vertex(f.grid.y);

% Electrod config
elec                = gen.Rho.elec;
[~,min_spacing]     = min(abs(f.grid.y-elec.spacing_y));
elec.spacing_y      = f.grid.y(min_spacing);
f.elec_spacing_y    = min_spacing-1;
f.elec_y_id         = f.elec_spacing_y*(elec.bufzone_y)+1 : f.elec_spacing_y : numel(f.grid.y_n)-elec.bufzone_y*f.elec_spacing_y;
elec.y              = f.grid.y_n(f.elec_y_id);
[~,f.elec_x_t_id]   = min(abs(bsxfun(@minus,f.grid.x_n,elec.x_t(:))),[],2);
[~,f.elec_x_r_id]   = min(abs(bsxfun(@minus,f.grid.x_n,elec.x_r(:))),[],2);
elec.x_t            = f.grid.x_n(f.elec_x_t_id);
elec.x_r            = f.grid.x_n(f.elec_x_r_id)';
[elec.X, elec.Y] = meshgrid([elec.x_t; elec.x_r],elec.y);
[f.elec_X_id, f.elec_Y_id] = meshgrid([f.elec_x_t_id; f.elec_x_r_id],f.elec_y_id);

% elec.config_max   = 3000;
elec.n              = (numel(f.elec_x_t_id)+numel(f.elec_x_r_id))*numel(f.elec_y_id);
elec                = config_elec(elec);



% Inverse Grid
i                   = gen.Rho.i;
i.grid.x_n          = f.grid.x_n(1: round((f.grid.nx+1)/(i.grid.nx+1)) :end);
i.grid.y_n          = f.grid.y_n(1: round((f.grid.ny+1)/(i.grid.ny+1)) :end); % cell center
assert(all(ismember(elec.x_t,i.grid.x_n )) & all(ismember(elec.x_r,i.grid.x_n )),'The grid of the inverse is not possible because the electrode have not corresponding position.')
assert(all(ismember(elec.y,i.grid.y_n )) ,'The grid of the inverse is not possible because the electrode have no corresponding position.')
i.grid.x            = i.grid.x_n(1:end-1)+diff(i.grid.x_n)/2;
i.grid.y            = i.grid.y_n(1:end-1)+diff(i.grid.y_n)/2;
i.grid.nx           = numel(i.grid.x);
i.grid.ny           = numel(i.grid.y);

i.elec_y_id         = find(sum(bsxfun(@eq,i.grid.y_n',elec.y),2))';
i.elec_x_t_id       = find(sum(bsxfun(@eq,i.grid.x_n',elec.x_t),2));
i.elec_x_r_id       = find(sum(bsxfun(@eq,i.grid.x_n',elec.x_r'),2));
[i.elec_X_id, i.elec_Y_id] = meshgrid([i.elec_x_t_id; i.elec_x_r_id],i.elec_y_id);

% Forward
f.header            = 'Forward';  % title of up to 80 characters
f.job_type          = 0;
f.filepath          = filepath;
f.readonly          = 0;
f.alpha_aniso       = gen.covar.range0(2)/gen.covar.range0(1);

% Rho value
% f                  = griddedInterpolant({grid.y,grid.x},rho_true,'nearest','nearest');
f.rho               = R_true; % f({grid_Rho.y,grid_Rho.x});
% f.filename       = 'gtrue.dat';
f.rho_avg           = mean(R_true(:));
f.rho_min           = min(R_true(:))/100;
f.rho_max           = max(R_true(:))*100;
f                   = Matlat2R2(f,elec); % write file and run forward modeling


if gen.forwardonly 
    gen.Rho.f  = f;
    gen.Rho.elec        = elec;
    filename = ['data_gen/data/FOR-', gen.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM'), '.mat'];
    if gen.saveit
    save(filename, 'K_true', 'grid_gen', 'gen') %'sigma',
    end
    keyboard
end

% Add some error to the observation
i.a_wgt = 0;%0.01;
i.b_wgt = 0.05;
% var(R) = (a_wgt*a_wgt) + (b_wgt*b_wgt) * (R*R)
f.output.resistancewitherror = i.a_wgt.*randn(numel(f.output.resistance),1) + (1+i.b_wgt*randn(numel(f.output.resistance),1)).*f.output.resistance;
%f.output.resistancewitherror(f.output.resistancewitherror>0) = -f.output.resistancewitherror(f.output.resistancewitherror>0);
%f.output.resistancewitherror(f.output.resistancewitherror<-10) = -10;


% fid = fopen([f.filepath 'R2_forward.dat'],'r');
% A = textscan(fid,'%f %f %f %f %f %f %f');fclose(fid);
% A{end-1}(2:end) = f.output.resistancewitherror;
% fid=fopen([f.filepath 'R2_forward.dat'],'w');
% A2=[A{:}];
% fprintf(fid,'%d\n',A2(1,1));
% for u=2:size(A2,1)
%     fprintf(fid,'%d %d %d %d %d %f %f\n',A2(u,:));
% end
% fclose(fid);


if 0==1
    figure(4); clf; hold on;
    [X,Y] = meshgrid(f.grid.x_n,f.grid.y_n);
    mesh(X,Y,0*X,'EdgeColor','b','facecolor','none')
    [X,Y] = meshgrid(i.grid.x_n,i.grid.y_n);
    mesh(X,Y,0*X,'EdgeColor','r','facecolor','none')
    [X,Y] = meshgrid(f.grid.x,f.grid.y);
    scatter(X(:),Y(:),'b')
    [X,Y] = meshgrid(i.grid.x,i.grid.y);
    scatter(X(:),Y(:),'r')
    scatter(f.grid.x_n(f.elec_X_id(:)),f.grid.y_n(f.elec_Y_id(:)),'xb')
    scatter(i.grid.x_n(i.elec_X_id(:)),i.grid.y_n(i.elec_Y_id(:)),'xr')
    view(2); axis tight; set(gca,'Ydir','reverse');
end


% Inverse
i.header            = 'Inverse';  % title of up to 80 characters
i.job_type          = 1;
i.filepath          = filepath;
i.readonly          = 0;
i.alpha_aniso       = f.alpha_aniso;
i.num_regions       = 1;
i.tolerance         = 1;
% i.rho               = f.rho;
i.rho_min           = f.rho_min;
i.rho_max           = f.rho_max;
i.rho_avg           = f.rho_avg;
i                   = Matlat2R2(i,elec);

%% Ouput Sigma
K.d             = 1./flipud(i.output.res);

if i.res_matrix ==  1 && any(~isnan(i.output.sen(:)))
    K.sen       = 1./flipud(i.output.sen);
elseif i.res_matrix ==  2 && any(~isnan(i.output.rad(:)))
    K.rad       = flipud(i.output.rad);
elseif i.res_matrix ==  3 && any(~isnan(i.output.Res(:)))
    K.res=i.output.Res;
    K.res(~i.output.inside,:)=[]; K.res(:,~i.output.inside(:))=[];
    K.res_out=i.output.Res;
    K.res_out(~i.output.inside,:)=[]; K.res_out(:,i.output.inside(:))=[];
    K.res_out = sum(K.res_out,2);
end
rmpath data_gen/R2
K.x = i.grid.x;
K.y = i.grid.y;
[K.X,K.Y] = meshgrid(K.x, K.y);

gen.Rho.i           = i;
gen.Rho.f           = f;
gen.Rho.elec        = elec;

%% * 4.*SAVING*

if gen.saveit
    filename = ['data_gen/data/GEN-', gen.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM'), '.mat'];
    save(filename, 'K_true','K','grid_gen', 'gen') %'sigma',
end

fprintf('  ->  finish in %g sec\n', toc)
end