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


%assert(elec.n<=300,'R2 cannot have more than 300 unique electrodes site')
%assert(elec.config_max<=6000,'R2 cannot have more than 6000 measurement')

%%
% * Defining for the seelec.nected method
%  A: +ve emmiteur electrode (C+)
%  B: -ve emmiteur electrode
%  M: +ve receiveur electrode (P+)
%  N: -ve receiveur electrode


emmiteur = find(ismember(elec.X,elec.x_t));
receiver = find(ismember(elec.X,elec.x_r));

% elec.data =nan(numel(emmiteur)/2*numel(receiver),4);
% i=1;
% for i_e1=1:numel(emmiteur)
%     for i_e2=i_e1+1:numel(emmiteur)
%         for i_r=1:numel(receiver)/2
%             elec.data(i,:) = [ receiver(i_r) receiver(i_r+numel(receiver)/2) emmiteur(i_e1) emmiteur(i_e2)];
%             i=i+1;
%         end
%     end
% end

[E, R] = meshgrid(emmiteur,receiver);

ref = ones(numel(E),1);

% [M N A B]
elec.data = [R(:) (elec.n+1)*ref E(:) (elec.n+2)*ref];

elec.pseudo_x = mean(elec.X(elec.data(:,[1 3])),2);
elec.pseudo_y = mean(elec.Y(elec.data(:,[1 3])),2);

% figure; hold on;
% for i=1:size(elec.data,1)
%     plot([elec.X(elec.data(i,3)) elec.X(elec.data(i,1))],[elec.Y(elec.data(i,3)) elec.Y(elec.data(i,1))])
% end

end
