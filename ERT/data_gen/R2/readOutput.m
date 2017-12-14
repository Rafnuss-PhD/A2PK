function output=readOutput(d)

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
    output.res=flipud(reshape(data(:,3),d.grid.ny,d.grid.nx));

    % read error result
    data            = dlmread([d.filepath dataset '_err.dat'],'',1,0);
    output.err      = data(:,1);
    output.pseudo   = data(:,3);
    output.wight    = data(:,5);
    
    if d.res_matrix==1 % 1-'sensitivity' matrix
        try
            data=dlmread([d.filepath dataset '_sen.dat']);
            output.sen=flipud(reshape(data(:,3),d.grid.ny,d.grid.nx));
        catch
            output.sen=NaN;
        end
    elseif d.res_matrix==2 % 2-true resolution matrix - diagonal element
        try
            data=dlmread([d.filepath dataset '_rad.dat']);
            output.rad=flipud(reshape(data(:,3),d.grid.ny,d.grid.nx));
        catch
            output.rad=NaN;
        end
    elseif d.res_matrix==3 % 2-true resolution matrix - Jacobian Matrix + roughness matrix
        try

            data=dlmread([d.filepath dataset '_sen.dat']);
            output.sen=flipud(reshape(data(:,3),d.grid.ny,d.grid.nx));
            
            n_parm = (d.numnp_y-1)*(d.numnp_x-1);
            n_obs = numel(d.pseudo_x);
            
            % Jacobian Matrix
            data=dlmread([d.filepath dataset '_J.dat']);
            J = data(2:end,:)';
            if numel(J)== n_parm * n_obs
                output.J=reshape( J, n_parm, n_obs);
            else
                output.J=reshape( J(J(:)~=0), n_parm, n_obs);
            end
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
            diagWd=zeros(n_obs,1);
            diagWd(sum(abs(output.J))~=0) = data.data(:,5);
            output.Wd=sparse(diag(diagWd));

             
            % Compute the sensitivity matrix
            sens = output.J*(output.Wd*output.Wd')*output.J';
            
            % Compute the Resolution Matrix
            filetext = fileread([d.filepath 'R2.out']);
            beg = strfind(filetext,'Alpha:');
            ed = strfind(filetext(beg:end),'RMS Misfit:');  
            alpha = str2double(filetext(beg(end)+6:beg(end)+ed(1)-3));
            
            output.Res = (sens + alpha * output.R)^(-1) * sens;

            % Compute the zone inside, that is removing the buffer zone
            output.inside=false( d.numnp_y-1, d.numnp_x-1);
            dx = (d.numnp_x-1-d.grid.nx)/2 +1;
            output.inside( 1:d.grid.ny, dx:end-dx+1) = true;
            
        catch
            warning('Not reading the jacobien, ressolution... matrices')
            keyboard;
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
    % output.re=flipud(reshape(data(:,3),d.grid.ny,d.grid.nx));
    
    data=dlmread([d.filepath 'R2_forward.dat'],'',1,0);
    assert(size(data,2)==7)
    output.pseudo=data(:,7);
    output.resistance=data(:,6);
    
    % Interpolation
    if numel(d.pseudo_x) == numel(output.pseudo)
        f=scatteredInterpolant(d.pseudo_y,d.pseudo_x,output.pseudo,'nearest','none');
        output.pseudo_interp = f({d.grid.y,d.grid.x});
        
%         if d.job_type == 1
%             f=scatteredInterpolant(d.pseudo_y,d.pseudo_x,output.err,'nearest','none');
%             output.err_interp = f({d.grid.y,d.grid.x});
%         end
    else
        output.pseudo_interp = nan(numel(d.grid.x),numel(d.grid.y));
    end

end




end
