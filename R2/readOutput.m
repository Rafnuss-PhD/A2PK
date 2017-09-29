function output=readOutput(d,grid_G)

if d.job_type == 1 % inverse solution
    
    dataset='f001';
    
    if exist([d.filepath dataset '_res.dat'], 'file')
        data=dlmread([d.filepath dataset '_res.dat']);
    elseif exist([d.filepath dataset '.' sprintf('%03d',d.max_iterations) '_res.dat'], 'file')
        data=dlmread([d.filepath dataset '.' sprintf('%03d',d.max_iterations) '.dat']);
    else 
        error('no file...')
    end
    
    % read resistivity result
    output.res=flipud(reshape(data(:,3),grid_G.ny,grid_G.nx));

    % read error result
    data            = dlmread([d.filepath dataset '_err.dat'],'',1,0);
    output.err      = data(:,1);
    output.pseudo   = data(:,3);
    output.wight    = data(:,5);
    
    if d.res_matrix==1 % 1-'sensitivity' matrix
        try
            data=dlmread([d.filepath dataset '_sen.dat']);
            output.sen=flipud(reshape(data(:,3),grid_G.ny,grid_G.nx));
        catch
            output.sen=NaN;
        end
    elseif d.res_matrix==2 % 2-true resolution matrix - diagonal element
        try
            data=dlmread([d.filepath dataset '_rad.dat']);
            output.rad=flipud(reshape(data(:,3),grid_G.ny,grid_G.nx));
        catch
            output.rad=NaN;
        end
    elseif d.res_matrix==3 % 2-true resolution matrix - Jacobian Matrix + roughness matrix
        try

            data=dlmread([d.filepath dataset '_sen.dat']);
            output.sen=flipud(reshape(data(:,3),grid_G.ny,grid_G.nx));
            
            n_parm = (d.numnp_y-1)*(d.numnp_x-1);
            n_obs = numel(d.pseudo_x);
            
            % Jacobian Matrix
            data=dlmread([d.filepath dataset '_J.dat']);
            J = data(2:end,:)';
            output.J=reshape( J(:), n_parm,n_obs);
            %imagesc(output.J)
            
            % Roughness matrix m-m_ref
            data=dlmread([d.filepath dataset '_Rindex.dat']);
            Rind = data(2:end,:);
            data=dlmread([d.filepath dataset '_R.dat']);
            R_raw = data(2:end,:);
            
            output.R=sparse(n_parm,n_parm);
            for i=1:size(Rind,1)
                for j=1:size(Rind,2)
                    if Rind(i,j)>0
                        output.R(Rind(i,1), Rind(i,j)) = R_raw(i,j);
                    end
                end
            end
            
            % Data Weight
            data=importdata([d.filepath dataset '_err.dat']);
            output.Wd=sparse(diag(data.data(:,5)));
             
            % Compute the sensitivity matrix
            sens = output.J*(output.Wd*output.Wd')*output.J';
            
            % Compute the Resolution Matrix
            alpha = input('Find the last alpha value in the .out file');
            output.Res = (sens + alpha * output.R)^(-1) * sens;

            
        catch
            output.sen=NaN;
            output.J=NaN;
            output.R=NaN;
            output.Res=NaN;
        end
    end
    
else
    data=dlmread([d.filepath 'forward_model.dat']);
    % output.x=unique(data(:,1));
    % output.y=-unique(data(:,2));  
    output.re=flipud(reshape(data(:,3),grid_G.ny,grid_G.nx));
    
    data=dlmread([d.filepath 'R2_forward.dat'],'',1,0);
    assert(size(data,2)==7)
    output.pseudo=data(:,7);
end


% Interpolation
if numel(d.pseudo_x) == numel(output.pseudo)
    f=scatteredInterpolant(d.pseudo_y,d.pseudo_x,output.pseudo,'nearest','none');
    output.pseudo_interp = f({grid_G.y,grid_G.x});
    
    if d.job_type == 1
        f=scatteredInterpolant(d.pseudo_y,d.pseudo_x,output.err,'nearest','none');
        output.err_interp = f({grid_G.y,grid_G.x});
    end
else
    output.pseudo_interp = nan(numel(grid_G.x),numel(grid_G.y));
end


end
