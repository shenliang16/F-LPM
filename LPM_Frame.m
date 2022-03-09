
function [Crrsp,time, c1]=LPM_Frame(X, Y, SIFTratio, K1, Rotation_thr, Scale_thr, lambda2)
global Crr_GT_4
%% parameters setting
lambda1   = 1;   
if ~exist('K1','var'); K1  = 14;   end
if ~exist('Rotation_thr','var');   Rotation_thr = pi/2;     end
if ~exist('Scale_thr','var');   Scale_thr = 1.5;     end
if ~exist('lambda2','var'); lambda2   = 0.5;   end
if ~exist('flag','var'); flag   = 0;   end

tau1      = 0.2;   tau2      = 0.2;   
%%
tic; 
Xt=[X(:, 1:2)']; 
Yt=[Y(:, 1:2)']; 
%% Frame-based Neighbor
    % Scaling and Rotation
    Scaling = Y(:,3)./X(:,3);                %     Scaling = Scaling*normal.xscale/normal.yscale; %     
    Theta = wrapToPi(Y(:,4) - X(:,4));   
    if ~isempty(SIFTratio)
        idx = SIFTratio2Idx(SIFTratio, 0.85);
    else
        idx = 1:length(Scaling);
    end
    [neighborX, neighborY] = FindNeighbors2Iter(Xt, Yt, idx, Scaling, K1);%[neighborX, neighborY] = FindNeighbors(Xt, Yt, Scaling, K1, 1);
    %Inliers±ÈÂÊ
%     if flag == -1
%         neighborIndex = zeros(K1+1, size(neighborX,2));
%         Ind_left = find(Scaling>1);     Ind_right = find(Scaling<1);
%         neighborIndex(:, Ind_left) = neighborX(1:K1+1,Ind_left);
%         neighborIndex(:, Ind_right) = neighborY(1:K1+1,Ind_right);
%         [c1.InlierRatio(1), c1.InlierRatio_inlier(1), c1.InlierRatio_outlier(1), c1.InlierNum_all(1), c1.InlierNum_inlier(1), c1.InlierNum_outlier(1)] = NeigInlierRatio3(neighborIndex,  Crr_GT_4, K1);
%     else
        c1 = [];
%     end

%     vec=Yt-Xt; d2=vec(1,:).^2+vec(2,:).^2; [p2, C] = LPM_cosF(neighborX, neighborY, 0.9, vec, d2, tau1, K1-2);
    [p2, ~, ~] = Frame_LPM_cosF(X, Y, neighborX, neighborY,  Scaling, Theta, Rotation_thr, Scale_thr, K1, tau1, 1);

%%  iteration 2
idx = find(p2 < lambda1);    K1 = 8;      Rotation_thr = pi/2;    Scale_thr = 1.4;  
if length(idx)>= K1+4
    [neighborX, neighborY] = FindNeighbors2Iter(Xt, Yt, idx, Scaling, K1);
            %         [neighborX, ~] = vl_kdtreequery(kdtreeX, Xt(:,idx), Xt, 'NumNeighbors', K1+1) ;
            %         [neighborY, ~] = vl_kdtreequery(kdtreeY, Yt(:,idx), Yt, 'NumNeighbors', K1+1) ;
            %         neighborX = idx(neighborX);
            %         neighborY = idx(neighborY);
    [p2, ~, c2] = Frame_LPM_cosF(X, Y, neighborX, neighborY,  Scaling, Theta, Rotation_thr, Scale_thr, K1, tau2, 2);
    idx = find(p2 <= lambda2);
end
time=toc;

Crrsp = idx;
end


function idx = SIFTratio2Idx(SIFTratio, SRthr)
    dist_in_sort=sort(SIFTratio);   
    midInd=min(round(length(dist_in_sort)/2),200);  
    threshold=max(dist_in_sort(midInd),SRthr); 
    idx = find(SIFTratio < threshold);
end



