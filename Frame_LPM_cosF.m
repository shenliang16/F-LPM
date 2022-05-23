function [p2, c1, c2] = Frame_LPM_cosF(X, Y, neighborX, neighborY,  Scaling, Theta, Rotation_thr, Scale_thr, K1, tau, iter)
    N = size(X, 1);
    c1 = 1e5*ones(1, N);
    c2 = 1e5*ones(1, N);             
    Low_thr = min(Scale_thr, 1/Scale_thr);      High_thr = max(Scale_thr, 1/Scale_thr); 
    ratio = 1;

%     idx1 = [idx(:); 0];     N = length(idx1);
%     neighborX(neighborX ==0) = N;
%     neighborY(neighborY ==0) = N;
%     neighborX = idx1(neighborX);
%     
%     Crr_GT_4 = [Crr_GT_4; 0];
%     Knn_ind(Knn_ind ==0) = length(Crr_GT_4);
%     neighborX
% 
%     neighborIndex = [neighborX(2:end+1, :); neighborY(2:end+1, :)];
%     [Nei_sort, index_sort] = sort(neighborIndex);
%     temp1 = diff(Nei_sort);
%     index_sort(temp1 == 0) = N+1;
%     neighborIndex(index_sort)

    for i = 1:N
        NeiX_i = nonzeros(neighborX(1:end,i)); 
        NeiY_i = nonzeros(neighborY(1:end,i));   
        % 去掉自己
        NeiX_i(NeiX_i == i) = [];
        NeiY_i(NeiY_i == i) = [];
        %% S1. 计算左右两侧 共同的邻域个数;
        [Idx_common, ~, ~] = intersect(NeiY_i, NeiX_i);

        %% S2. 计算被测点和邻域的旋转、缩放是否相近
        Zoom = Scaling(Idx_common);
        Rot = Theta(Idx_common);
        Zoom_ratio = Zoom./Scaling(i);
        Rot_diff = wrapToPi(Rot - Theta(i));
        Idx_common_RS_idx = Zoom_ratio> Low_thr & Zoom_ratio<High_thr  &  abs(Rot_diff)<Rotation_thr ;
        Idx_common_RS = Idx_common(Idx_common_RS_idx)';
        c1(i) = K1 - length(Idx_common_RS);   % 加1 是因为前面考虑了自身，即neighborY(1:end,i)是从1开始取的  
    
        c2(i) = 0;
        if iter >= 1 && ~isempty(Idx_common_RS)
            Neighbors_X = X(Idx_common_RS, 1:2)'; %     Neighbors_X = Y(NeiX_i, 1:2)'; 
            Neighbors_Y = Y(Idx_common_RS, 1:2)'; %     Neighbors_Y = Y(NeiY_i, 1:2)'; 
            %% S3. 利用Frame信息，进行预变换，并计算c2;
            Ctr_X = X(i, 1:2)';                                 
            vec1 =  Neighbors_X - Ctr_X; 
        %     Rotation = [cos(Theta(i)), sin(Theta(i)); -sin(Theta(i)), cos(Theta(i))];   % Rotation
        %     NeighborsX_Transformed = (Ctr_X + Scaling(i) .* (Rotation * vec1))';                      % Error  
            Ctr_Y = Y(i, 1:2)';                                 
            vec2 =  Neighbors_Y - Ctr_Y; 
                % 去掉与被测点任一端点相同的点，
                Same_Idx = all(vec2 == 0) | all(vec2 == 0);
                Same = sum(Same_Idx);
                % 去掉邻居中端点相同的点，
                [~, Idx_sort] = sort(vec2(1,~Same_Idx));% 按大小，列重排列
                vec2_del_Same = vec2(:, Idx_sort);
                Same = Same + sum(all(diff(vec2_del_Same,[],2)==0));
            

            d1 = sum(vec1.^2, 1);
            d2 = sum(vec2.^2, 1);
            
            Rotation = [cos(Theta(i)), sin(Theta(i)); -sin(Theta(i)), cos(Theta(i))];   %             Error = (DisY(:,:,i) - Scaling(i) .* (Rotation * DisX(:,:,i)))'; 
            RotatedVec1 = Rotation * vec1;
            cos_sita_transfored = sum(RotatedVec1.*vec2, 1) ./ sqrt(d1.*d2);
%                 cos_sita = sum(vec1.*vec2, 1) ./ sqrt(d1.*d2);
%                 cos_sita_transfored = cos(acos(cos_sita) + Theta(i));
        %             if any(abs(cos_sita_transfored2-cos_sita_transfored)>0.001) 
        %                 pause(1)
        %             end
            ratio = Scaling(i) * sqrt(d1 ./ d2);
            c2i = abs(log(cos_sita_transfored.*ratio)) > tau;
            %c2i = abs(cos_sita_transfored.*ratio - 1) > tau;
            c2(i) = sum(c2i) + Same;
        end
        
%         if any(i==[1	6	12	82	85	132	177	185	186	191	198])
%            pause(0.001) 
%         end
    end
    %% 换顺序
    p2 = (c1+c2)/K1;%
end
