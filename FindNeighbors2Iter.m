function [neighborX, neighborY] = FindNeighbors2Iter(Xt, Yt, idx, Scaling, K1)
%% 去除重复点（由多对一的匹配引起），这在大缩放下尤为重要
% Xt_idx = Xt(:, idx);
% Yt_idx = Yt(:, idx);
[~, iaX, icX] = unique(Xt(1, idx));
[~, iaY, icY] = unique(Yt(1, idx));
Xt1 = Xt(:, idx(iaX));
Yt1 = Yt(:, idx(iaY));

%% KNN
    Kmax = min([5*K1, length(iaX)-1, length(iaY)-1]);%Kmax=round(max(min(100,(N-1)/2),K1));%ff=(ff-min(ff))./(max(ff)-min(ff));
    Kmax_1 = Kmax+1;
    kdtreeX = vl_kdtreebuild(Xt1); 
    [neighborX, DisXt] = vl_kdtreequery(kdtreeX, Xt1, Xt, 'NumNeighbors', Kmax_1);
    kdtreeX = vl_kdtreebuild(Yt1); 
    [neighborY, DisYt] = vl_kdtreequery(kdtreeX, Yt1, Yt, 'NumNeighbors', Kmax_1);
    neighborX = iaX(neighborX);
    neighborY = iaY(neighborY);
%     neighborX = idx(iaX(neighborX));
%     neighborY = idx(iaY(neighborY));

%% 比例阈值 和 最小数量
    Amp = 1.2;                      % 邻域放大比例
    num_min = round(Amp*K1);        % 待选邻域的最小数量
    max_num = 5*K1;
    
%% 以左侧Xt为基准
%     Xt_zero_flag = (DisXt(1, :)~=0);
    Ind_left = find(Scaling>1);
    Nei_Y_sel = neighborY(2:end, Ind_left);
    neighborY(K1+2:end, :) = 0;                 % 默认只有K1个邻域，其与删除
    Sel_FlagY = DisYt(2:end, Ind_left) < repmat( Amp * Scaling(Ind_left)'.*DisXt(K1+1, Ind_left),   Kmax, 1);
    Sel_FlagY(1:num_min, :) = true;
    num_Right = sum(Sel_FlagY);
    Nei_Y_sel = uint32(double(Nei_Y_sel).*double(Sel_FlagY));
    neighborY(2:end, Ind_left) = Nei_Y_sel;
    
    
%% 以又侧Yt为基准
    Ind_right = find(Scaling<=1);               % 满足以右侧为基准的序号
    Nei_X_sel = neighborX(2:end, Ind_right);    % 取选择的邻域
    neighborX(K1+2:end, :) = 0;                 % 默认只有K1个邻域，其与删除
    Sel_FlagX = Scaling(Ind_right)'.* DisXt(2:end, Ind_right) < repmat( Amp * DisYt(K1+1, Ind_right),   Kmax, 1);
    Sel_FlagX(1:num_min, :) = true;             % 在对侧，至少选择min_num个邻域
    num_Left = sum(Sel_FlagX);                  % Frame改善后的邻域数量
    Nei_X_sel = uint32(double(Nei_X_sel).*double(Sel_FlagX));
    neighborX(2:end, Ind_right) = Nei_X_sel;
    
%%  小图中重复点恢复
    if length(iaX)/length(idx)<0.8 %(如果重复比例大于1/0.8，则恢复)
        neighborX = recoverSamePoint(icX, Ind_left, neighborX);
    end
    if length(iaY)/length(idx)<0.8 %(如果重复比例大于1/0.8，则恢复)
        neighborY = recoverSamePoint(icY, Ind_right, neighborY);
    end
    
%%  在原集合中的索引
    idx1 = [idx(:); 0];     N = length(idx1);
    neighborX(neighborX ==0) = N;
    neighborY(neighborY ==0) = N;
    neighborX = idx1(neighborX);
    neighborY = idx1(neighborY);
end
   
%% 小图中重复点恢复
function neighborY1 = recoverSamePoint(icY, Ind_right, neighborY)
    % 最大5倍
    [m, n] = size(neighborY);
    neighborY1 = zeros(m*4, n);    neighborY1 = [neighborY; neighborY1];
    
    [Sort_icY, idx_sort_icY] = sort(icY);
%     N = length(icY);
%     diff_Sort_icY = diff(Sort_icY);
%     start_end_idx = find(diff_Sort_icY); 
%     % 每组的开始和结束的idx
%     start_idx_group = [1; start_end_idx]; 
%     end_idx_group =  [start_end_idx; N]; 
%     % 属于那一组
%     Nei_group = icY(nonzeros(neighborY(2:end, 1137))); 
%     idx_sort_icY(start_idx_group(Nei_group))
%     end_idx_group(Nei_group)

    % 按重复分组
    Num_Groups = max(Sort_icY);
    Group = cell(1, Num_Groups);
    for i = 1:Num_Groups %分组，每组中为同一个点
        temp = Sort_icY==i;
        Group{i} = idx_sort_icY(temp);
    end
    
    num_nei_XY1 = size(neighborY1, 1);
    for i = 1:length(Ind_right)
        Ind_i =  Ind_right(i);
        temp = [];
        temp_idx = icY(nonzeros(neighborY(2:end, Ind_i)));  
        for j = 1 : length(temp_idx)
             temp = [temp; Group{temp_idx(j)}];
        end
        num_nei_max = min(length(temp)+1, num_nei_XY1);
        neighborY1(2:num_nei_max, Ind_i) = temp(1:num_nei_max-1);
    end
end
    

    
%     A2 = unique(icY);
%     [n, bins_idx] = histc(icY, A2);
% [SortYt, idx_sortYt] = sort(Yt(1, idx));
% same_flag = diff(SortYt)==0;
% cumsum_same = cumsum(same_flag);
% idx_same_Yt = idx_sortYt(same_flag);