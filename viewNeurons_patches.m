function [A_new, C_new] = viewNeurons_patches(A,C,b,f,bkgnd,ctr,options, ind, C2, folder_nm)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.
%   folder_nm: string, the folder to output images neuron by neuron.

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

% allocate outputs
% A_new = A; C_new = C;

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(A, 2);
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save resulted images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end

% obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1));

ind_del = false(size(ind));     % indicator of deleting neurons
% ctr = obj.estCenter();      %neuron's center
gSiz = options.gSiz;        % maximum size of a neuron
if isempty(gSiz)
    gSiz = 2*options.gSig+1;
end
% time
T = size(C, 2);
t = 1:T;
if ~isnan(options.fr)
    t = t/options.fr;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

%% start viewing neurons
figure('position', [100, 100, 1024, 512]);
m=1; maxInd = length(ind);
while and(m>=1, m<=length(ind))
    %% full-frame view
    subplot(3,3,[1 2 4 5]); thr = .95;
%     imagesc(reshape(A(:, ind(m)), options.d1, options.d2));
    imagesc(bkgnd); %colormap(gray);
    hold on;
    A_temp = full(reshape(A(:,ind(m)),options.d1,options.d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind_temp] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        x_cont =  reshape(A_temp,options.d1,options.d2);
        y_cont = [0,0]+A_temp(ind_temp(ff));
        [~,ww] = contour(x_cont,y_cont,'LineColor','k');
        ww.LineWidth = 2;
    end
%     imagesc(reshape(A(:, ind(m)), options.d1, options.d2)); 
    
    axis equal; axis off; hold off;
    if ind_del(m)
        title(sprintf('Neuron %d', ind(m)), 'color', 'r');
    else
        title(sprintf('Neuron %d', ind(m)));
    end
    %% zoomed in with background
    subplot(3,3,3);
%     imagesc(bkgnd);
%     axis equal; axis off;
%     x0 = ctr(ind(m), 2);
%     y0 = ctr(ind(m), 1); 
%     hold on;
%     [~,ww] = contour(x_cont,y_cont,'LineColor','k');
%     ww.LineWidth = 2;
%     xlim(x0+[-gSiz, gSiz]*2);
%     ylim(y0+[-gSiz, gSiz]*2); hold off;
    
    
    %% zoomed-in view
    subplot(3,3,6); 
    imagesc(reshape(A(:, ind(m))/max(A(:,ind(m))), options.d1, options.d2));
    axis equal; axis off;
    x0 = ctr(ind(m), 2);
    y0 = ctr(ind(m), 1);
    xlim(x0+[-gSiz, gSiz]*2);
    ylim(y0+[-gSiz, gSiz]*2);
    
    
    %% temporal components
    subplot(3,3,7:9);cla;hold off;
    if ~isempty(C2)
        plot(t, C2(ind(m), :)*max(A(:, ind(m))), 'linewidth', 2); hold on;
        plot(t, C(ind(m), :)*max(A(:, ind(m))), 'r');
    else
        cell_ts =  C(ind(m), :)*max(A(:, ind(m)));
        plot(t, cell_ts); %hold on;
    end
    xlabel(str_xlabel);
    
    %% save images
    if save_img
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
        m = m+1;
    else
        fprintf('Neuron %d of %d, keep(k, default)/delete(d)/split(s)/trim(t)/gif (g)/delete all(da)/backward(b)/end(e):    ', size(C,1), ind(m));
        temp = input('', 's');
        if temp=='d'
            ind_del(m) = true;
            m = m+1;
        elseif strcmpi(temp, 'b')
            m = m-1;
        elseif strcmpi(temp, 'da')
            ind_del(m:end) = true;
            break;
        elseif strcmpi(temp, 'k')
            ind_del(m) = false;
            m= m+1;
        elseif strcmpi(temp, 's')
            try
                subplot(3,3,6);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                tmpA = A(:, ind(m));
                A(:, end+1) = tmpA.*tmp_ind(:);
                C(end+1, :) = C(ind(m), :);
                A(:, ind(m)) = tmpA.*(1-tmp_ind(:));
                maxInd = maxInd+1;
%                 S(end+1, :) = S(ind(m), :);
%                 C_raw(end+1, :) = C_raw(ind(m), :);
%                 P.kernel_pars(end+1, :) = P.kernel_pars(ind(m), :);
            catch
                sprintf('the neuron was not split\n');
            end
        elseif strcmpi(temp, 't')
            try
                subplot(3,3,6);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                A(:, ind(m)) = A(:, ind(m)).*tmp_ind(:);
            catch
                sprintf('the neuron was not trimmed\n');
            end
            
        elseif strcmpi(temp,'g')
            try
               subplot(3,3,7:9); hold on;
               fprintf('select timepoints to make gif\n');
               [xrange, ~]= ginput(2);
               xrange=sort(round(xrange*options.fr));
               
               
               tt = xrange(2)-xrange(1)+1;
               gif = reshape(A*C(:,xrange(1):xrange(2)) + b*f(:,xrange(1):xrange(2)),...
                   options.d1,options.d2,tt);
%                    gif([1:round(max([x0-gSiz*2 1])), round(min([x0+gSiz*2,options.d1])):end],:) = 0;
%                    gif(:,[1:round(max([y0-gSiz*2 1])), round(min([y0+gSiz*2,options.d2])):end])=0;
                while 1
                   Cmax = full(max(cell_ts));
                   subplot(3,3,7:9); tp1=plot([t(xrange(1)) t(xrange(1))],[0 Cmax],'r');
                   tp2 = plot([t(xrange(2)) t(xrange(2))],[0 Cmax],'r');
                   indLine = plot([t(xrange(1)), t(xrange(1))],...
                                [0 Cmax],'k');
%                    gline = animatedline('MaximumNumPoints',tt);
                   
                   tf = input('play(1)/continue(0)?');
                   if tf ==0
                       delete(tp1); delete(tp2); delete(indLine);
                       break
                   else
                       subplot(3,3,3);
                       gifhandle=imagesc(gif(:,:,1)); hold on; axis equal; axis off;
                       [~,ww] = contour(x_cont,y_cont,'LineColor','k');
                       ww.LineWidth = 2; 
                       xlim(x0+[-gSiz, gSiz]*2);
                       ylim(y0+[-gSiz, gSiz]*2); 
                            
                       drawnow;  hold off;
                       for i = 1:tt
                           subplot(3,3,3);
                            set(gifhandle,'CData',gif(:,:,i)); drawnow;
                            
                            subplot(3,3,7:9); delete(indLine);
%                             addpoints(gline,Cmax,t(xrange(1)+i-1));
%                             dt = (xrange(2)-xrange(1))/tt;
                            indLine = plot([t(xrange(1)+(i-1)), t(xrange(1)+(i -1))],...
                                [0 Cmax],'k'); drawnow;

                            
                       end
                   end
               end
            catch ME
                
                fprintf('gif fail');
                rethrow(ME);
            end
            
        elseif strcmpi(temp, 'e')
            break;
        else
            m = m+1;
        end
    end
end
A_new = A(:,~ind_del);
C_new = C(~ind_del,:);
% if save_img
%     cd(cur_cd);
% else
%     obj.delete(ind(ind_del));
%     obj.Coor = obj.get_contours(0.9);
% end

