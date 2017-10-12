function createProtocoldat(d)

%% OPEN FILE
fileID = fopen([d.filepath 'protocol.dat'],'w');

%%
fprintf(fileID,'%d    << num_ind_meas',d.num_ind_meas);

if d.job_type==1
    if d.inverse_type==3
        fprintf(fileID,'%d  %d    << j, node\n',[d.j_e,d.node]);
    end
else
    assert(d.num_ind_meas==numel(d.j_p) && d.num_ind_meas==size(d.elec,1),'j_p and elec need to be equal to num_ind_meas')
    fprintf(fileID,'    << j, elec\n');
    fprintf(fileID,'%3d      %3d  %3d  %3d  %3d\n', [d.j_p; d.elec']);
end

fclose(fileID);
end