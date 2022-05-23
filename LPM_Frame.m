
function Crrsp=LPM_Frame(X, Y, SIFTratio, K, Rotation_thr, Scale_thr, lambda, tau)
%% Inputs
% X: Feature set 1
% Y: Feature set 2
% SIFTratio:  NNSR values in the ratio test
% K:  Neighborhood size, k
% Rotation_thr:  rotation threshold, ¦È
% Scale_thr:  zoom threshold, ¦Ç
% lambda:  balancing weight, or decision threshold, ¦Ë
% tau: Topology Consensus threshold, ¦Ó

lambda1   = 1;    
Xt=[X(:, 1:2)']; 
Yt=[Y(:, 1:2)']; 
%% Frame-based Neighbor
    % Scaling and Rotation
    Scaling = Y(:,3)./X(:,3);              
    Rotation = wrapToPi(Y(:,4) - X(:,4));   
    if ~isempty(SIFTratio)
        idx = SIFTratio2Idx(SIFTratio, 0.8);
    else
        idx = 1:length(Scaling);
    end
    [neighborX, neighborY] = FindNeighbors3Iter(Xt, Yt, idx, Scaling, K, Rotation, Rotation_thr, Scale_thr);

    [p2, ~, ~] = Frame_LPM_cosF(X, Y, neighborX, neighborY,  Scaling, Rotation, Rotation_thr, Scale_thr, K, tau, 1);

%%  iteration 2
idx = find(p2 < lambda1);   Rotation_thr = max(Rotation_thr*0.8, pi/8);  Scale_thr = max(Scale_thr*0.8, 1.2);  %Rotation_thr = pi/4;    Scale_thr = 1.4;  
if length(idx)>= K+4
    [neighborX, neighborY] = FindNeighbors3Iter(Xt, Yt, idx, Scaling, K, Rotation, Rotation_thr, Scale_thr);
    [p2, ~, c2] = Frame_LPM_cosF(X, Y, neighborX, neighborY,  Scaling, Rotation, Rotation_thr, Scale_thr, K, tau, 2);
    idx = find(p2 <= lambda);
end

Crrsp = idx;
end


function idx = SIFTratio2Idx(SIFTratio, SRthr)
    dist_in_sort=sort(SIFTratio);   
    midInd=max(min(round(length(dist_in_sort)/4), 200), 20);  
    threshold=max(dist_in_sort(midInd),SRthr); 
    idx = find(SIFTratio < threshold);
end



