%% Create R2.in and protocol.dat
%
% Comment are from the Readme Manual (2.7a)
%
% INPUT (generated with this script):
% * R2.in           : geometry informaton 
% * protocol.dat    : index, 4 electrodes index 
% ( * mesh.dat      : for triangulare meshing)
%
% OUPUT:
% * R2.out            : main log exectution
% * electrodes.dat    : x,y-coordinate of electrodes
% * electrodes.vtk    : idem in vtk format
% FORWARD OUTPUT
% * R2_forward.dat    : similare to protocol.dat + calculated resistances +
% calculated apparent resistivities
% * forward_model.dat : x, y, resis, log10(resis)
% INVERSE

function [resistance, pseudo] = Matlat2R2min(d)
% d is either the inverser (i) or forward (f) structure

% RESISITIVITY
d.value =  d.rho(:)' ;

% PROTOCOL
% createProtocoldat(d)

% CREATE THE FILES
% createR2in(d)
content_start=strfind(d.content,'<< elem_1, elem_2, value')+24;
content_end=strfind(d.content,'<< num_xy_poly')-8;
content = [d.content(1:content_start ) sprintf('%5d  %5d  %5d\n',[d.elem_1;d.elem_2;d.value]) d.content(content_end:end)];
fId = fopen( [d.filepath 'R2.in'], 'w' ) ;
fwrite( fId, content);
fclose( fId );

% RUN .EXE
pwd_temp = pwd;
cd(d.filepath);
[~,~]=system('R2.exe');
cd(pwd_temp);

% OUPUT
data=dlmread([d.filepath 'R2_forward.dat'],'',1,5);
resistance=data(:,1);
pseudo=data(:,2);
end




