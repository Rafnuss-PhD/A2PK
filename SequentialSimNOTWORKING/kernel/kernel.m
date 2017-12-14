function kern = kernel(Prim, Sec, parm, plotit)
% kernel computes the non-parametric density
% INPUT :
%         - Prim    : Primary variable
%         - Sec     : Secondary variable
%         - parm    : parameter structure
% 
% OUTPUT :
%         - brandwidth   : Primary variable (1D-X.n elts)
%         - density   : Primary variable (1D-X.n elts)
%         - x   : Primary variable (1D-X.n elts)
%         - y   : Primary variable (1D-X.n elts)
%         - n   : Number of ?? (1D-X.n elts)


if isfield(parm, 'dens') % Kernel Given
    kern.dens=parm.dens;
    
    assert(isfield(parm, 'axis_sec'))
    assert(isfield(parm, 'axis_prim'))
    kern.axis_sec = parm.axis_sec(:);
    kern.axis_prim = parm.axis_prim(:);
    
    if(~isfield(parm, 'prior'))
        kern.prior = sum(kern.dens,2);
    else
        kern.prior = parm.prior;
    end
    
    [X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
    kern.XY = [X(:),Y(:)];
    
    
elseif Prim.n==0 % Not possible to built one.
    warning('No hard data !')
    assert(all(parm.p_w(1,:)==0),'The Secondary cannot be used !')
    kern.axis_prim=linspace(-5,5,100)';
    kern.axis_sec=kern.axis_prim;
    kern.dens=ones(100,100);
    [X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
    kern.XY = [X(:),Y(:)];
    
    
else % Build from hard data

    if (~isfield(parm, 'axis_prim') || ~isfield(parm, 'axis_sec') )
        % Scott's rules
        dP = 3.5*std(Prim.d(:))*Prim.n^(-1/3);
        dS = 3.5*std(Sec.d(:))*numel(Sec.d)^(-1/3);
        kern.axis_prim = (min(Prim.d(:))-.2*range(Prim.d(:))):dP:(max(Prim.d(:))+.2*range(Prim.d(:)));
        kern.axis_sec = (min(Sec.d(:))-.2*range(Sec.d(:))):dS:(max(Sec.d(:))+.2*range(Sec.d(:)));
    end
    [X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
    kern.XY = [X(:),Y(:)];

    % Find the colocated value of Prim in Sec.
    if Prim.n==numel(Sec.xy) && all(Prim.y(:)==Sec.y) && all(Prim.x(:)==Sec.x)
        krange.min = [min(Sec.d(:))-.2*range(Sec.d(:)) min(Prim.d(:))-.2*range(Prim.d(:))];
        krange.max = [max(Sec.d(:))+.2*range(Sec.d(:)) max(Prim.d(:))+.2*range(Prim.d(:))];
        data=[Sec.d(:) Prim.d(:)];
    else
        kern.id=nan(Prim.n,1);
        for i=1:Prim.n
            [~,idy] = min((Prim.y(i)-Sec.y).^2);
            [~,idx] = min((Prim.x(i)-Sec.x).^2);
            kern.id(i)= sub2ind(size(Sec.d),idy, idx);
        end
        data=[Sec.d(kern.id) Prim.d];
    end 

    kern.dens = ksdensity(data,kern.XY);
    kern.prior = sum(kern.dens,2);
    
end 




% Plot
if plotit
    figure; hold on;
    imagesc(kern.axis_sec, kern.axis_prim, kern.dens)
    scatter(data(:,1), data(:,2),'xk');
    ylabel('Primary variable'); xlabel('Secondary variable');
    axis tight; colorbar('Northoutside')
    % plot([0 20],[0 20],'--k')
    % set(gca,'xticklabel',{''});  xlabel(''); 
    % set(gca,'yticklabel',{''});  ylabel('');
    keyboard
end


% Check for error
[~,id] = min( (kern.axis_sec-min(Sec.d(:))).^2 );
assert(sum(kern.dens(:,id))>10^(-10) ,'The current joint pdf is not suitable for the whole secondary variable range. add more collocated data.');

end

% OLD METHOD
% krange.min = [min(Sec.d(:))-.2*range(Sec.d(:)) min(Prim.d(:))-.2*range(Prim.d(:))];
% krange.max = [max(Sec.d(:))+.2*range(Sec.d(:)) max(Prim.d(:))+.2*range(Prim.d(:))];
% 
% % Find the colocated value of Prim in Sec.
% if Prim.n==numel(Sec.xy) && all(Prim.y(:)==Sec.y) && all(Prim.x(:)==Sec.x)
%     krange.min = [min(Sec.d(:))-.2*range(Sec.d(:)) min(Prim.d(:))-.2*range(Prim.d(:))];
%     krange.max = [max(Sec.d(:))+.2*range(Sec.d(:)) max(Prim.d(:))+.2*range(Prim.d(:))];
%     data=[Sec.d(:) Prim.d(:)];
% else
%     kern.id=nan(Prim.n,1);
%     for i=1:Prim.n
%         [~,idy] = min((Prim.y(i)-Sec.y).^2);
%         [~,idx] = min((Prim.x(i)-Sec.x).^2);
%         kern.id(i)= sub2ind(size(Sec.d),idy, idx);
%     end
%     data=[Sec.d(kern.id) Prim.d];
% end
% 
% % Number of point in the joint pdf
% n = 2^8;
% 
% 
% % Define the range of the grid used in the density. Range have been
% % computed as input data
% assert( all(max(data) < krange.max) )
% assert( all(min(data) > krange.min) )
% 
% 
% % Bandwidth proposed by Foster and Bowman
% l_sec=std(data(:,1))* Prim.n^(-1/6);
% l_prim=std(data(:,2))* Prim.n^(-1/6);
% 
% % Compute the density
% [~,kern.dens,x,y]=kde2d_pr(data,n,krange.min,krange.max,l_sec,l_prim);
% 
% % recover and re-name the output
% kern.axis_sec = x(1,:)';
% kern.axis_prim = y(:,1); % convert the meshgrid to simple vector
% kern.daxis_prim = (kern.axis_prim(2)-kern.axis_prim(1));
% kern.n = n;
% 
% % Normalizing...
% kern.dens(kern.dens<0) = 0;
% kern.dens = kern.dens./sum(kern.dens(:));
% 
% % Compute the marginal distribution
% kern.prior = hist(Prim.d(:),kern.axis_prim)';
% kern.prior = kern.prior./sum(kern.prior(:));
% 
% end
