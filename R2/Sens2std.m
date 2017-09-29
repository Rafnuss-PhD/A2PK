function Sigma =Sens2std(Sigma,sigma,grid)

Sigma_sen_log       = log10(Sigma.sen);
% normalized sensitivity. then because we have the side effect
Sigma_sen_norm      = 1.4*(Sigma_sen_log-min(Sigma_sen_log(:)))/(max(Sigma_sen_log(:)-min(Sigma_sen_log(:))));

% 
% Sigma_d_well        = Sigma.d(ismember(sigma.y,grid.Y)&ismember(sigma.x,grid.X));
% err                 = sigma.d-Sigma_d_well;
% predictor           = Sigma_sen_norm(ismember(sigma.y,grid.Y)&ismember(sigma.x,grid.X)) .* Sigma_d_well;
% 
% [pred_so,id] = sort(predictor);
% err_so = err(id);
% 
% n=floor(numel(predictor)/10);
% data=nan(n,2);
% for i=1:n
%     idx = (i-1)*10+(1:10);
%     data(i,1) = mean( pred_so(idx) );
%     data(i,2) = std( err_so(idx) ) + abs(mean(err_so(idx)));
% end
% 
% figure; hold on;
% plot(predictor,err,'o')
% plot(data(:,1),data(:,2),'d')
% 
% 
% x=[ones(n,1) data(:,1)]\data(:,2);
% 

% Sigma.std  = x(1) + x(2).* Sigma_sen_norm .* Sigma.d;


%Sigma.std           = .015 + .08*Sigma_sen_norm.*Sigma.d;

Sigma.std           = .01 + .3*Sigma_sen_norm.*Sigma.d;
% figure; imagesc(Sigma.std);colorbar
% figure; 
% subplot(4,1,1); imagesc(Sigma.d);colorbar
% subplot(4,1,2); imagesc(sigma_true); colorbar;
% subplot(4,1,3); imagesc(sigma_true-Sigma.d); colorbar;
% subplot(4,1,4); imagesc(Sigma.std); colorbar;

clear Sigma_sen_log Sigma_sen_norm 
end


%% TRUE FIELD
% 
% err                 = sigma_true-Sigma.d;
% predictor           = Sigma_sen_norm .* Sigma.d;
% 
% [pred_so,id] = sort(predictor(:));
% err_so = err(id);
% 
% dn=1000;
% n=floor(numel(predictor)/dn);
% data=nan(n,2);
% for i=1:n
%     idx = (i-1)*dn+(1:dn);
%     data(i,1) = mean( pred_so(idx) );
%     data(i,2) = std( err_so(idx) ) + mean(abs(err_so(idx)));
% end
% 
% figure; hold on;
% plot(predictor(:),err(:),'o')
% plot(data(:,1),data(:,2),'d')
% 
% 
% x=[ones(n,1) data(:,1)]\data(:,2);
% 
% 
% Sigma.std  = x(1) + x(2).* Sigma_sen_norm .* Sigma.d;
% 
