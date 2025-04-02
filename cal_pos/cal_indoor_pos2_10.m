function [true_pos2,lat0,lon0,const_i,const_j] = cal_indoor_pos2_10()   % 设置轨迹的起始点
% 初始化pos2数组

const_i = 95;
const_j = 40;
m = 4; 
pos2 = [1,1.5];
% 设置六边形的顶点
vertices = [1,1.5;24,1.5;24,0.3;1,0.3;1,1.5];
%vertices = [10,4;16,4;16,16;10,16;10,4];
%vertices = [const_i-m/2,const_j-m/2;const_i+m/2,const_j-m/2;const_i+m/2,const_j+m/2;const_i-m/2,const_j+m/2;const_i-m/2,const_j-m/2];%正方形
% 设置步长
step = 0.2;

% 计算每条边的方向向量和长度
edges = diff(vertices, 1, 1);
edge_lengths = sqrt(sum(edges.^2, 2));

% 生成六边形轨迹
for i = 1:size(edges, 1)
    % 计算当前边的单位方向向量
    dir_vector = edges(i, :) / norm(edges(i, :));
    
    % 计算当前边的点数
    num_points = floor(edge_lengths(i) / step);
    
    % 生成当前边的点
    for j = 1:num_points
        point = vertices(i, :) + j * step * dir_vector;
        pos2(end+1, :) = point;
    end
end


    
     % 绘制四边形
    figure;
    plot(pos2(:,1), pos2(:,2), '-o');
    axis([0 26 0 2]);
    grid on;
    xlabel('X');
    ylabel('Y');
    title('任意四边形闭合轨迹');

    pos2 = round(pos2,1); % 将 pos3 中的数据保留到小数点后一位
    true_pos2 = pos2;
    true_pos2(:,[1 2]) = true_pos2(:,[2 1]); % 交换第一列和第二列，第一列为纬度数据，第二列为经度数据
    lat0 = true_pos2(1,1);
    lon0 = true_pos2(1,2);
    % 如果需要，可以将 pos2 保存到文件
    % save('pos2.mat', 'true_pos2');
    % 如果需要，可以将 pos2 保存到文件
    % save('pos2.mat', 'true_pos2');
end