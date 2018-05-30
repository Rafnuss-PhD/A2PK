function createR2in(d)

fileID = fopen([d.filepath 'R2.in'],'w');

%% BASIC GENERAL SETTING
assert(ischar(d.header),'header need to be a string')
assert(length(d.header)<80,'header need to be smaller than 80 characteres')
fprintf(fileID,'%s\n\n\n',d.header);

assert(d.job_type==0 || d.job_type==1,'job_type need to be ''0'' or ''1''')
assert(d.mesh_type==3 || d.mesh_type==4 || d.mesh_type==5,'mesh_type need to be ''3'', ''4'' or ''5''')
assert(d.flux_type==2.0 || d.flux_type==3,'flux_type need to be ''2.0'' or ''3.0''')
assert(d.singular_type==0 || d.singular_type==1,'singular_type need to be ''0'' or ''1''')
assert(d.res_matrix==0 ||d.res_matrix==1 || d.res_matrix==2 || d.res_matrix==3,'flux_type need to be ''1'', ''2'' or ''3''')
fprintf(fileID,'%d  %d  %1.1f  %d  %d    << job_type, mesh_type, flux_type, singular_type, res_matrix\n\n',d.job_type, d.mesh_type, d.flux_type, d.singular_type, d.res_matrix);

%% MESH SETTING
if d.mesh_type==3 %Trigular Mesh
    fprintf(fileID,'%s    << scale\n\n',d.scale);
else
    assert(mod(d.numnp_x,1)==0 && mod(d.numnp_y,1)==0,'numnp_x and numnp_y need to be integer')
    fprintf(fileID,'%d  %d    << numnp_x, numnp_y\n\n',d.numnp_x,d.numnp_y);
    assert(d.numnp_x==length(d.xx))
    fprintf(fileID,'%5.3f ', d.xx); fprintf(fileID,'    << xx\n\n');
    if d.mesh_type==4
        assert(d.numnp_x==length(d.topog))
        fprintf(fileID,'%5.3f ', d.topog); fprintf(fileID,'    << topog\n\n');
    end
    assert(d.numnp_y==length(d.yy))
    fprintf(fileID,'%5.3f ', d.yy); fprintf(fileID,'    << yy\n\n');
end
fprintf(fileID,'\n');

%% RESISTIVITY
fprintf(fileID,'%d    << num_regions    << elem_1, elem_2, value\n',d.num_regions);
if d.num_regions==0
    %fprintf(fileID,'    << file_name\n');
    fprintf(fileID,'%s\n\n',d.file_name);
else
    fprintf(fileID,'%5d  %5d  %5d\n',[d.elem_1;d.elem_2;d.value]);
    fprintf(fileID,'\n\n');
end

%% INVERSE SOLUTION
if d.job_type==1
    if d.mesh_type==4 || d.mesh_type==5
        fprintf(fileID,'%d  %d    << patch_size_x, patch_size_y\n',d.patch_size_x,d.patch_size_y);
        if d.patch_size_x==0 &&  d.patch_size_y==0
            fprintf(fileID,'%d  %d    << num_param_x, num_param_y\n',d.num_param_x,d.num_param_y);
            fprintf(fileID,'%d  %d    << npxstart, npx\n',d.npxstart,d.npx);
            fprintf(fileID,'%d  %d    << npystart, npy\n',d.npystart,d.npy);
        end
    end
    fprintf(fileID,'%d  %d   << inverse_type, target_decrease \n',d.inverse_type, d.target_decrease);
    if d.inverse_type==3
        fprintf(fileID,'%d    << qual_ratio\n',d.inverse_type);
        fprintf(fileID,'%d  %d    << rho_min, rho_max\n',d.rho_min,d.rho_max);
    else
        fprintf(fileID,'%d  %d    << data_type, reg_mode\n',d.data_type,d.reg_mode);
        
        if d.reg_mode==0 || d.reg_mode==2
            fprintf(fileID,'%d  %d  %d  %d    << tolerance, max_iterations, error_mod, alpha_aniso\n',d.tolerance, d.max_iterations, d.error_mod, d.alpha_aniso);
        else
            fprintf(fileID,'%f  %d  %d  %f  %f    << tolerance, max_iterations, error_mod, alpha_aniso, alpha_s\n',d.tolerance, d.max_iterations, d.error_mod, d.alpha_aniso, d.alpha_s);
        end
        fprintf(fileID,'%f  %f  %f  %f   << a_wgt, b_wgt, rho_min, rho_max\n', d.a_wgt, d.b_wgt, d.rho_min, d.rho_max);
    end
    fprintf(fileID,'\n\n');
end

%% OUPUT REGION
assert(d.num_xy_poly==0 || (d.num_xy_poly==numel(d.x_poly) && d.num_xy_poly==numel(d.y_poly)) ,'x_poly and y_poly need to be equal to num_xy_poly')
fprintf(fileID,'%d    << num_xy_poly    << x_poly, y_poly\n',d.num_xy_poly);
fprintf(fileID,'%8.3f  %8.3f\n',[d.x_poly;d.y_poly]);
fprintf(fileID,'\n\n');

%% ELECTRODE
fprintf(fileID,'%d    << num_electrodes',d.num_electrodes);
if d.job_type==1 && d.inverse_type==3
    fprintf(fileID,'%d  %d    << j, node\n',d.j_e,d.node);
else
    assert(d.num_electrodes==numel(d.j_e) && d.num_electrodes==numel(d.column) && d.num_electrodes==numel(d.row),'j, column and row need to be equal to num_electrodes')
    fprintf(fileID,'    << j, column, row\n');
    fprintf(fileID,'%3d  %3d  %3d\n', [d.j_e;d.column;d.row]);
end
fprintf(fileID,'\n\n');

%% CLOSE
fclose(fileID);
end
