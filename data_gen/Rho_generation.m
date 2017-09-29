function [Sigma, gen] = Rho_generation(grid, rho_true, gen)
% MEAS_Rho_GRID simulate the measurement of low resolution of the grid from the g
% field (rho_true). Two method are implemented (data or ERT simulation). In
% the ERT option, ERT is simulated with ERT2D package. We still need to
% implement the inverse modeling and noise to get the low resolution G
% field
% INPUT:
%       - grid:     grid of the matrix to generate (see Import_dat for more gen)
%       - rho_true:   variogram (see mGstat or gStat)
%       - method:   choice of method : 1. borehole | 2. random generated
%       - plotit:   1 or 0 to disply a plot or not
%
% OUTPUT:
%       - Rho.d:      grid measurement (G)
%       - Rho.std     std error associated to G
%
% Author: Raphael Nussbaumer
% date : January 2014
% need to do : add assert() for input, flexible number of input var

switch gen.Rho.method
    case 'Paolo' % read data
        
        load('data.mat')
        [X,Y]=meshgrid(x,y);
        F = scatteredInterpolant(X(:),Y(:),sigma_obs(:),'linear','nearest');
        Rho.d = F(grid.X,grid.Y);
        F = scatteredInterpolant(X(:),Y(:),sigma_obs_err(:),'linear','nearest');
        Rho.err = F(grid.X,grid.Y);
        assert(~any(isnan(Rho.d(:))),'error');
        assert(~any(isnan(Rho.err(:))),'error');
        clear sigma_obs sigma_obs_err sigma_true
        
        Rho.std = Rho.err/100.*Rho.d;
        
        if plotit
            figure;hold on
            subplot(3,1,1); imagesc(grid.x,grid.y,rho_true);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}')
            subplot(3,1,2); imagesc(grid.x,grid.y,Rho.d);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('G')
            subplot(3,1,3); imagesc(grid.x,grid.y,Rho.std);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('Rho_{std}')
        end
        
    case 'noise' % add simple noise
        
        noise=std(rho_true(:))*randn(size(rho_true));
        Rho_noised=rho_true+1*noise;
        
        grid_new.nx=20;
        grid_new.ny=5;
        grid_new.x=linspace(grid.x(1),grid.x(end),grid_new.nx);
        grid_new.y=linspace(grid.y(1),grid.y(end),grid_new.ny);
        [grid_new.X, grid_new.Y]=meshgrid(grid_new.x,grid_new.y);
        dist= @(x,y) min(sqrt(sum(bsxfun(@minus,y,x).^2,2)));
        Rho.d=nan(size(rho_true));
        for i=1:grid.nx
            for j=1:grid.ny
                [~,idx]=dist([grid.x(i) grid.y(j)],[grid_new.X(:) grid_new.Y(:)]);
                Rho.d(j,i)=idx;
            end
        end
        for i=1:grid_new.nx*grid_new.ny
            Rho.d(Rho.d==i)=mean(rho_true(Rho.d==i));
        end
        
        Rho.std=1.5*ones(size(Rho.d));
        
        
        
        if plotit
            figure;hold on
            subplot(2,2,1); imagesc(grid.x,grid.y,rho_true);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}')
            subplot(2,2,2); imagesc(grid.x,grid.y,Rho_noised);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+noise')
            subplot(2,2,3); imagesc(grid.x(1:ratio_x:end), grid.y(1:ratio_y:end), Rho_upscaled);shading flat;  xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+noise upscaled')
            subplot(2,2,4); imagesc(grid.x,grid.y,Rho.d); shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+scale upscaled and then downscaled in a grid')
        end
        
    case 'RESINV3D' %% NOT WORKING

        
        % conductivity
        s=rho_true';
        % grid size
        para.dx=grid.dx*ones(grid.nx,1);
        para.dz=grid.dy*ones(grid.ny,1);
        
        % electrode number, position, configuration
        elec.dx=.5; % electrode spacing in unit
        elec.n=ceil(grid.x(end)/elec.dx);
        elec.x=(-elec.n/2:elec.n/2)*elec.dx;
        
        [data,pos,n_config]=configuration('dipole-dipole',elec.n,grid.y(end)/elec.dx,1000,0);
        
        para.srcloc=[elec.x(data(:,1))' grid.dy*ones(n_config,1) elec.x(data(:,2))' grid.dy*ones(n_config,1)];
        para.recloc=[elec.x(data(:,3))' grid.dy*ones(n_config,1) elec.x(data(:,4))' grid.dy*ones(n_config,1)];
        para.srcnum=1:n_config;
        
        
        tic;
        [Para] = get_2_5Dpara(para.srcloc,para.dx,para.dz,[],0,para.recloc,para.srcnum);%
        toc;
        
        [dobs,U] = dcfw2_5D(s,Para);
        
        % save('Forward_calc.mat','dobs','U','Para','para','n_config','s')
        
        %         figure;
        %         i=1; % visualized the ith
        %         imagesc(grid.x,grid.y,reshape(U(:,i),length(para.dx),length(para.dz))')
        %         figure; scatter(pos(:,1),pos(:,2),[],dobs)
        %         xlabel('Electrode position');ylabel('depth')
        %         set(gca, 'YDir', 'reverse');
        
        stop
        script_RESINVM3D(dobs,para,n_config,mean2(s))
        
        rmpath FW2_5D RESINVM3D
        
    case 'R2' % R2
        grid_Rho = gen.Rho.grid;
        elec = gen.Rho.elec;
        dmin = gen.Rho.dmin;
        
        dmin.readonly = 0;
        
        addpath data_gen/R2
        
        %where to write stuff
        dmin.filepath       = 'data_gen/data/IO-file/';
        dmin.res_matrix     = gen.Rho.dmin.res_matrix; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
        delete([dmin.filepath '*']);
        
        % grid_Rho.nx           = 300;
        % grid_Rho.ny           = 30;
        grid_Rho.nxy          = grid_Rho.ny*grid_Rho.nx;
        grid_Rho.x            = linspace(grid.x(1),grid.x(end),grid_Rho.nx); % cell center
        grid_Rho.y            = logspace(log10(grid.y(1)+5),log10(grid.y(end)+5),grid_Rho.ny)-5; % cell center
        grid_Rho.x_n          = [grid_Rho.x(1)-(grid_Rho.x(2)-grid_Rho.x(1))/2  grid_Rho.x(1:end-1)+diff(grid_Rho.x)/2  grid_Rho.x(end)+(grid_Rho.x(end)-grid_Rho.x(end-1))/2]; % this are nodes and not element (center)
        grid_Rho.y_n          = [grid_Rho.y(1)-(grid_Rho.y(2)-grid_Rho.y(1))/2  grid_Rho.y(1:end-1)+diff(grid_Rho.y)/2  grid_Rho.y(end)+(grid_Rho.y(end)-grid_Rho.y(end-1))/2]; % this are nodes and not element (center)
        assert((grid_Rho.ny+16)*(grid_Rho.nx+32)<=30000,'the number of cell cannot be greater than 30''0000')
        gen.Rho.grid          = grid_Rho;
        
        % Electrod config
        % elec.spacing      = 2; % in grid_Rho.dx unit
        elec.x              = grid_Rho.x_n(1:elec.spacing:end);
        elec.n              = numel(grid_Rho.x_n(1:elec.spacing:end)); % max is 300
        % elec.config_max   = 3000; % max is 6000
        elec.method         = 'dipole-dipole';
        elec.depth_max      = ceil(elec.n*grid_Rho.y(end)/grid_Rho.x(end));
        elec.selection      = 5; %1: k-mean, 2:iterative removal of the closest neighboohood 3:iterative removal of the averaged closest point 4:voronoi 5:random
        elec                = config_elec(elec); % create the data configuration.
        gen.Rho.elec        = elec; % put the elec struct in the dmin
        
        
        % Forward
        dmin.header         = 'Forward';  % title of up to 80 characters
        dmin.job_type       = 0;
        % Interpolation
        f                   = griddedInterpolant({grid.y,grid.x},rho_true,'nearest','nearest');
        dmin.rho_true       = f({grid_Rho.y,grid_Rho.x});
        dmin.filename       = 'gtrue.dat';
        dmin.num_regions    = 1+numel(dmin.rho_true);
        dmin.rho_min        = min(rho_true(:));
        dmin.rho_avg        = mean(rho_true(:));
        dmin.rho_max        = max(rho_true(:));
        gen.Rho.f           = Matlat2R2(grid_Rho,dmin,elec); % write file and run forward modeling
        
        % Inverse
        dmin.header         = 'Inverse';  % title of up to 80 characters
        dmin.job_type       = 1;
        dmin.num_regions    = 1;
        % dmin.tolerance    = 0.1;
        gen.Rho.i           = Matlat2R2(grid_Rho,dmin,elec);
        
        
        % Interpolation
        Rho.d_raw           = flipud(gen.Rho.i.output.res);
        f                   = griddedInterpolant({grid_Rho.y,grid_Rho.x},Rho.d_raw,'nearest','nearest');
        Rho.d               = f({grid.y,grid.x});
        Sigma.d             = 1000./Rho.d;
        Sigma.d_raw         = 1000./Rho.d_raw;
        
        if dmin.res_matrix ==  1 && any(~isnan(gen.Rho.i.output.sen(:)))
            Rho.sen_raw         = flipud(gen.Rho.i.output.sen);
            Sigma.sen_raw       = 1000./Rho.sen_raw;
            f                   = griddedInterpolant({grid_Rho.y,grid_Rho.x},Sigma.sen_raw,'nearest','nearest');
            % f                   = griddedInterpolant({gen.Rho.grid.y,gen.Rho.grid.x},Sigma.sen_raw,'nearest','nearest');
            Sigma.sen           = f({grid.y,grid.x});
            % Sigma.sen           = f({grid_gen.y,grid_gen.x});
            
        elseif dmin.res_matrix ==  2 && any(~isnan(gen.Rho.i.output.rad(:)))
            Sigma.rad_raw       = flipud(gen.Rho.i.output.rad);
            f                   = griddedInterpolant({grid_Rho.y,grid_Rho.x},Sigma.rad_raw,'nearest','nearest');
            Sigma.rad           = f({grid.y,grid.x});
        elseif dmin.res_matrix ==  3 && any(~isnan(gen.Rho.i.output.Res(:)))
            inside=false( numel(gen.Rho.i.yy)-1, numel(gen.Rho.i.xx)-1);
            dx = (numel(gen.Rho.i.xx)-1-numel(grid_Rho.x))/2 +1;
            inside( 1:numel(grid_Rho.y), dx:end-dx+1) = true;
            
            Sigma.Res=gen.Rho.i.output.Res;
            Sigma.Res(~inside,:)=[]; Sigma.Res(:,~inside(:))=[];
            Sigma.Res_out=gen.Rho.i.output.Res;
            Sigma.Res_out(~inside,:)=[]; Sigma.Res_out(:,inside(:))=[];
            Sigma.Res_out = sum(Sigma.Res_out,2);
            

            Sigma.res_raw       = flipud();
            f                   = griddedInterpolant({grid_Rho.y,grid_Rho.x},Sigma.res_raw,'nearest','nearest');
            Sigma.res           = f({grid.y,grid.x});
        
        else
            Sigma.sen           = ones(size(grid.X));
        end
        
        rmpath data_gen/R2

end

Sigma.x_raw = grid_Rho.x;
Sigma.y_raw = grid_Rho.y;
Sigma.dx_raw = [diff(Sigma.x_raw(1:2)) diff(Sigma.x_raw(1:end-1)+.5*diff(Sigma.x_raw)) diff(Sigma.x_raw(end-1:end))];
Sigma.dy_raw = [diff(Sigma.y_raw(1:2)) diff(Sigma.y_raw(1:end-1)+.5*diff(Sigma.y_raw)) diff(Sigma.y_raw(end-1:end))];


% grid_Rho.x = Sigma.x_raw;
% grid_Rho.y = Sigma.y_raw;
% grid=grid{end};

Sigma.x = grid.x;
Sigma.y = grid.y;
Sigma.xy = grid.xy;



end
