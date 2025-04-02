% 磁地图理论参考仿真器
% 根据输入经纬度高度，构建以输入经纬度为中心的大小20*20、精度0.1的磁地图
clear all; clc; close all;

load('./BPIT/fx.mat')
load('./BPIT/xfx.mat')
addpath('./multi_algorithm')

disp("*1*输入初始数据******************************************************")
% 输入经纬度和高度
input_lat = input('请输入纬度：');
input_lon = input('请输入经度：');
alt = input('请输入高度（km）: ');
sizelat = 20;
sizelon = 20;
resolutionlat = 0.1;
resolutionlon = 0.1;

% 如果输入为空，设置默认值
if isempty(input_lat)
    input_lat = 0;
end
if isempty(input_lon)                          
    input_lon = 0;
end
if isempty(alt)
    alt = 0;
end
if isempty(sizelat)
    sizelat =180 ;
end
if isempty(sizelon)
    sizelon = 360;
end
if isempty(resolutionlat)
    resolutionlat = 18;
end
if isempty(resolutionlon)
    resolutionlon= 36;
end

% 计算经度和纬度范围
center_lat = input_lat;
center_lon = input_lon;
top_lat = center_lat + sizelat / 2;
bottom_lat = center_lat - sizelat / 2;
left_lon = center_lon - sizelon / 2;
right_lon = center_lon + sizelon / 2;
disp("*2*生成地磁图数据****************************************************")
% 生成经度和纬度网格
lats = bottom_lat:resolutionlat:top_lat;
lons = left_lon:resolutionlon:right_lon;
len_lats = length(lats);
len_lons = length(lons);

if len_lats > len_lons   
    lons = linspace(lons(1), lons(end), len_lats);
else
    lats = linspace(lats(1), lats(end), len_lons);
end
 len_max=length(lats);

% 固定日期为2019年12月31日
date = datenum(2019, 12, 31);

% 获取用户输入位置的地磁场信息
[bx, by, bz] = igrf(date, input_lat, input_lon, alt);%东、北、天地磁分量，量纲：nT

% 计算地磁场总量
tmi = sqrt(bx^2 + by^2 + bz^2);
    
% 计算磁偏角和磁倾角
magnetic_declination = atan2d(by, bx);
magnetic_inclination = atan2d(bz, sqrt(bx^2 + by^2));  % 转换为角度

% 打印地磁场分量、总量、磁偏角、磁倾角、水平分量、偏航角%%%%%？？？？？？
fprintf('磁场分量 (X轴): %.4fnT\n', bx);
fprintf('磁场分量 (Y轴): %.4fnT\n', by);
fprintf('磁场分量 (Z轴): %.4fnT\n', bz);
fprintf('地磁场总量： %.4fnT\n', tmi);
fprintf('磁航向： %.4f degrees\n', magnetic_declination);
fprintf('磁倾角： %.4f degrees\n', magnetic_inclination);

% 初始化地磁场强度、磁偏角和磁倾角的矩阵
magnetic_mod_xy=ones(length(lats), length(lons));
magnetic_mod = ones(length(lats), length(lons));
magnetic_dec = ones(length(lats), length(lons));
magnetic_inc = ones(length(lats), length(lons));
magnetic_x=ones(length(lats), length(lons));
magnetic_y=ones(length(lats), length(lons));
magnetic_z=ones(length(lats), length(lons));
magnetic_alts=alt*ones(length(lats), length(lons));
magnetic_xy_mod = ones(length(lats), length(lons),length(lats)*length(lons));
rmse_matrix=ones(length(lons),length(lats));
[X, Y] = meshgrid(lats, lons);

% 计算整个区域内的地磁场信息
for i = 1:length(lats)
    for j = 1:length(lons)
        [x, y, z] = igrf(date, lats(i), lons(j), alt);    
        magnetic_mod(i, j) = sqrt(x^2 + y^2 + z^2);
        magnetic_dec(i, j) = atan2d((y), (x));
        magnetic_inc(i, j) = atan2d((z), sqrt((x)^2 + (y)^2));
        magnetic_x(i,j)=x;
        magnetic_y(i,j)=y;
        magnetic_z(i,j)=z;
        magnetic_mod_xy(i,j)=sqrt((x)^2 + (y)^2);
        magnetic_xy_mod(j,i,1) = lons(j); 
        magnetic_xy_mod(j,i,2) = lats(i); 
        magnetic_xy_mod(j,i,3) = magnetic_mod(i,j); 
    
    end
end
%save('magnetic_mod_gauss.mat','magnetic_mod');
for i = 1:length(lats)
    for j = 1:length(lons)
        % 提取周围20x20矩阵大小(2*2经纬度范围)方形范围内的点的模值
        j_min = max(1, j-20);
        j_max = min(len_lons, j+20);
        i_min = max(1, i-20);
        i_max = min(len_lats, i+20);
        local_region = magnetic_mod(i_min:i_max,j_min:j_max);     
        % 计算均值
        mean_value = mean(local_region, 'all');
        % 计算RMSE
        rmse = sqrt(mean((local_region - mean_value).^2, 'all'));   
        % 将RMSE值保存到新矩阵中
        rmse_matrix(j, i) = rmse;
    end
end

% 绘制二维坐标系
imagesc(rmse_matrix);
colorbar; % 显示颜色条
title('每个坐标点的均方根误差值');
xlabel('X 坐标');
ylabel('Y 坐标');

magnetic_map=[X;Y;magnetic_alts;magnetic_x;magnetic_y;magnetic_z;magnetic_mod_xy;magnetic_mod;magnetic_dec;magnetic_inc];
% 使用 griddata 插值生成平滑的网格数据
smooth_resolutionlat = resolutionlat/10;
smooth_resolutionlon = resolutionlon/10;
[lats_smooth, lons_smooth] = meshgrid(bottom_lat:smooth_resolutionlat:top_lat, left_lon:smooth_resolutionlon:right_lon);
magnetic_inc_smooth = griddata(lats, lons, magnetic_inc, lats_smooth, lons_smooth, 'cubic');
magnetic_dec_smooth = griddata(lats, lons, magnetic_dec, lats_smooth, lons_smooth, 'cubic');
magnetic_mod_smooth=griddata(lats, lons, magnetic_mod, lats_smooth, lons_smooth, 'cubic');

lats_smooth1=lats_smooth(1,:);
lons_smooth1=lons_smooth(1,:);

%% 生成干扰源
% 添加时空特性的干扰
disp("*3*添加干扰源****************************************************")
sizeMatrix = size(magnetic_mod,1);
% 计算原磁场矩阵的平均值
averageMagneticField = mean(magnetic_mod(:));
averageMagneticField_dec = mean(magnetic_dec(:));
averageMagneticField_inc = mean(magnetic_inc(:));
averageMagneticField_x = mean(magnetic_x(:));
averageMagneticField_y = mean(magnetic_y(:));
averageMagneticField_z = mean(magnetic_z(:));
averageMagneticField_xy = mean(magnetic_mod_xy(:));
% 定义用于调整干扰强度的因子
amplitudeFactor = 0.001;      % 调整正弦波的振幅因子
%stdDeviationFactor = 0.00001;  % 调整高斯分布的标准差因子
stdDeviationFactor = 0.001;  %  调整高斯分布的标准差因子
% 正弦波的时变特性
frequency = 0.05;         % 正弦波频率
amplitude = amplitudeFactor * averageMagneticField;  % 调整正弦波振幅
amplitude_dec = amplitudeFactor * averageMagneticField_dec;  % 调整正弦波振幅
amplitude_inc = amplitudeFactor * averageMagneticField_inc;  % 调整正弦波振幅
amplitude_x = amplitudeFactor * averageMagneticField_x;  % 调整正弦波振幅
amplitude_y = amplitudeFactor * averageMagneticField_y;  % 调整正弦波振幅
amplitude_z = amplitudeFactor * averageMagneticField_z;  % 调整正弦波振幅
amplitude_xy = amplitudeFactor * averageMagneticField_xy;  % 调整正弦波振幅
temporalSineWave = amplitude * sin(2 * pi * frequency * linspace(0, 1, sizeMatrix));
temporalSineWave_dec = amplitude_dec * sin(2 * pi * frequency * linspace(0, 1, sizeMatrix));
temporalSineWave_inc = amplitude_inc * sin(2 * pi * frequency * linspace(0, 1, sizeMatrix));
temporalSineWave_x = amplitude_x * sin(2 * pi * frequency * linspace(0, 1, sizeMatrix));
temporalSineWave_y = amplitude_y * sin(2 * pi * frequency * linspace(0, 1, sizeMatrix));
temporalSineWave_z = amplitude_z * sin(2 * pi * frequency * linspace(0, 1, sizeMatrix));
temporalSineWave_xy = amplitude_xy * sin(2 * pi * frequency * linspace(0, 1, sizeMatrix));

% 高斯分布的空间特性
meanValue = 0;             % 高斯分布均值
stdDeviation = stdDeviationFactor;  % 调整高斯分布标准差
stdDeviation_dec = stdDeviationFactor * averageMagneticField_dec;  % 调整高斯分布标准差
stdDeviation_inc = stdDeviationFactor * averageMagneticField_inc;  % 调整高斯分布标准差
stdDeviation_x = stdDeviationFactor * averageMagneticField_x;  % 调整高斯分布标准差
stdDeviation_y = stdDeviationFactor * averageMagneticField_y;  % 调整高斯分布标准差
stdDeviation_z = stdDeviationFactor * averageMagneticField_z;  % 调整高斯分布标准差
stdDeviation_xy = stdDeviationFactor * averageMagneticField_xy;  % 调整高斯分布标准差


spatialGaussianNoise = stdDeviation * randn(sizeMatrix) + meanValue;
spatialGaussianNoise_dec = stdDeviation_dec * randn(sizeMatrix) + meanValue;
spatialGaussianNoise_inc = stdDeviation_inc * randn(sizeMatrix) + meanValue;
spatialGaussianNoise_x = stdDeviation_x * randn(sizeMatrix) + meanValue;
spatialGaussianNoise_y = stdDeviation_y * randn(sizeMatrix) + meanValue;
spatialGaussianNoise_z = stdDeviation_z * randn(sizeMatrix) + meanValue;
spatialGaussianNoise_xy = stdDeviation * randn(sizeMatrix) + meanValue;
% 将时空特性的干扰添加到磁场数据中
magnetic_TemporalSpatial_Noise = magnetic_mod + temporalSineWave + spatialGaussianNoise;
magnetic_TemporalSpatial_Noise_dec = magnetic_dec + temporalSineWave_dec + spatialGaussianNoise_dec;
magnetic_TemporalSpatial_Noise_inc = magnetic_inc + temporalSineWave_inc + spatialGaussianNoise_inc;
magnetic_TemporalSpatial_Noise_x = magnetic_x + temporalSineWave_x + spatialGaussianNoise_x; 
magnetic_TemporalSpatial_Noise_y = magnetic_y + temporalSineWave_y + spatialGaussianNoise_y;
magnetic_TemporalSpatial_Noise_z = magnetic_z + temporalSineWave_z + spatialGaussianNoise_z;
magnetic_TemporalSpatial_Noise_xy = magnetic_mod_xy + temporalSineWave_xy + spatialGaussianNoise_xy;

%创建加入干扰后的磁地图模值纹理图
figure("Name","加入干扰源后的磁地图模值");
% 绘制磁地图模值纹理图
worldmap([bottom_lat top_lat], [left_lon right_lon]);
geoshow(lats_smooth, lons_smooth, magnetic_TemporalSpatial_Noise, 'DisplayType', 'texturemap');
view(2);
colormap("jet");
geoshow(input_lat, input_lon, "Marker","*",'Color','black', 'MarkerSize', 10);
geoshow('landareas.shp', 'FaceAlpha', 0);
%geoshow([point_x1, point_x2], [point_y, point_y], 'DisplayType', 'line', 'Color', 'black','LineWidth', 2);
% xlabel('经度');
xlabel({'经度','(a) 加入干扰源后的磁地图模值纹理图'},'FontSize',12,'Fontname', '黑体','FontWeight','bold');
ylabel('纬度');
title('加入干扰源后的磁地图模值纹理图');

N = size(magnetic_mod,1);
M = size(magnetic_mod,2);
dataset = zeros(N,M,4);

dataset(:,:,1) = magnetic_TemporalSpatial_Noise;
dataset(:,:,2) = magnetic_TemporalSpatial_Noise_dec;
dataset(:,:,3) = magnetic_TemporalSpatial_Noise_inc;
dataset(:,:,4) = magnetic_TemporalSpatial_Noise_xy;

%% 生成运动轨迹
disp("*4*生成实际运动轨迹同时计算并存储步长和方向数据*********************************")
%load('./pos_data/record2_pos2.mat');       %三角形      匹配结果B
%load('./pos_data/record3_pos2.mat');       %三角形      匹配结果A
%load('./pos_data/record4_pos2.mat');       %长方形      匹配结果C
%load('./pos_data/record5_pos2.mat');       %三角形      匹配结果B   
%load('./pos_data/record6_pos2.mat');       %三角形      匹配结果B   
%load('./pos_data/record7_pos2.mat');       %三角形      匹配结果B
load('./pos_data/record8_square.mat');      %正方形      匹配结果A
true_pos2 = [true_pos2(:,1)+bottom_lat,true_pos2(:,2)+left_lon];
true_pos1 = [true_pos2;true_pos2;true_pos2];
%true_pos1 = true_pos2;
init_pos1_lat = true_pos2(1,1);
init_pos1_lon = true_pos2(1,2);

disp("*5*沿实际轨迹在加入干扰源的地磁图上采样并构建磁指纹******************************")

%% 读取受干扰后该运动轨迹上的地磁数据——曲折闭合曲线
number = size(true_pos1,1);
magnetic_pos1_noise=ones(number,4);
for i=1:number
    lat_point = round((true_pos1(i,1)-bottom_lat)*10+1);    %读取轨迹坐标的在磁地图中的行数和列数
    lon_point = round((true_pos1(i,2)-left_lon)*10+1);
    magnetic_pos1_noise(i,1)=magnetic_mod(lat_point,lon_point);
    magnetic_pos1_noise(i,2)=magnetic_dec(lat_point,lon_point);
    magnetic_pos1_noise(i,3)=magnetic_inc(lat_point,lon_point);
    magnetic_pos1_noise(i,4)=magnetic_mod_xy(lat_point,lon_point);
end

%% 生成轨迹数据近似的步长数据矩阵和方向数据矩阵——曲折闭合曲线
M =size(magnetic_pos1_noise,1);
sl1=zeros(M-1,1);
deg1=zeros(M-1,1);
for i=1:M-1
    sl1(i)=sqrt((true_pos1(i+1,1)-true_pos1(i,1)).^2+(true_pos1(i+1,2)-true_pos1(i,2)).^2);
    y=true_pos1(i+1,1)-true_pos1(i,1);
    x=true_pos1(i+1,2)-true_pos1(i,2);
    theta = atan2d(y,x); % 两点之间的方向角，单位为度
    deg1(i) = theta; % 将方向角转换为度，存储到 deg2 中
    
end

% 设定步长和航向的初始值
sl1_IMU = sl1; % 初始步长
deg1_IMU = deg1; % 初始航向角

% 参数设置
step_noise_mean = 0.001; % 步长噪声的均值，单位可以为米
step_noise_std = 0.001; % 步长噪声的标准差
heading_noise_mean = 0.1; % 航向角噪声的均值，单位为度                    
heading_noise_std = 0.1; % 航向角噪声的标准差

step_drift_rate = 0.0002; % 步长漂移率，单位为米/步
heading_drift_rate = 0.001; % 航向角漂移率，单位为度/步

%step_drift_rate = 0; % 步长漂移率，单位为米/步
%heading_drift_rate = 0; % 航向角漂移率，单位为度/步

% 添加噪声和漂移
for i = 2:length(sl1)
    % 添加非零均值的高斯噪声
    sl1_IMU(i) = sl1(i) + step_noise_mean + step_noise_std * randn();
    deg1_IMU(i) = deg1(i) + heading_noise_mean + heading_noise_std * randn();

    % 累计漂移
    sl1_IMU(i) = sl1_IMU(i) + step_drift_rate * i;
    deg1_IMU(i) = deg1_IMU(i) + heading_drift_rate * i;
    
    % 添加随机飘移
    sl1_IMU(i) = sl1_IMU(i) + 0.01 * (rand() - 0.5);
    deg1_IMU(i) = deg1_IMU(i) + 0.1 * (rand() - 0.5);
end

sl1 = sl1_IMU;
deg1 = deg1_IMU;

% 初始化路径矩阵，预先分配空间
n_steps = length(sl1);
ins_pos = zeros(n_steps+1, 2); % 第一行为初始位置，后面为每一步的位置
ins_pos(1, :) = [init_pos1_lat, init_pos1_lon]; % 存储初始位置

% 将角度转换为弧度
radians = deg1*pi/180;
%deg1 = radians;
% 计算每一步的轨迹位置
for i = 1:n_steps
    % 当前步长和方向
    step_length = sl1(i);
    angle = radians(i);
    
    % 计算新位置
    delta_y = step_length * sin(angle); % y方向增量
    delta_x = step_length * cos(angle); % x方向增量
    
    % 更新位置
    new_y = ins_pos(i, 1) + delta_y;
    new_x = ins_pos(i, 2) + delta_x;
    
    % 存储新位置
    ins_pos(i+1, :) = [new_y, new_x];
end

%% 利用粒子滤波算法计算地磁匹配结果
number_mag = 4; 
deg1 = deg1*pi/180;
disp("*6*分别利用4种算法计算地磁匹配结果数据***************************************")
disp("*6.1*计算第一个轨迹的地磁匹配结果数据")
[pos_PF_1,pos_AOFA_1,pos_EKPF_1,pos_AIPF_1,rsm_pf1_1,rsm_pf2_1,rsm_pf3_1,rsm_pf4_1,mean_tm_1,mean_rsm_pf_1,mean_mean_pf_1] =...
 prediction_pos1(init_pos1_lon, init_pos1_lat, sl1, deg1, magnetic_pos1_noise(:, [1:number_mag]), dataset(:, :, [1:number_mag]), fx, xfx, true_pos1, bottom_lat, left_lon);

%        轨迹1                         轨迹2                      轨迹3
stdpf1_1 = std(rsm_pf1_1); 
stdpf2_1 = std(rsm_pf2_1); 
stdpf3_1 = std(rsm_pf3_1); 
stdpf4_1 = std(rsm_pf4_1); 

[min_value_1_1, index_1_1] = min(rsm_pf1_1); % 找到最小值及其索引
pos_PF_1 = pos_PF_1(:,:,index_1_1);

[min_value_2_1, index_2_1] = min(rsm_pf2_1); % 找到最小值及其索引
pos_AOFA_1 = pos_AOFA_1(:,:,index_2_1);

[min_value_3_1, index_3_1] = min(rsm_pf3_1); % 找到最小值及其索引
pos_EKPF_1 = pos_EKPF_1(:,:,index_3_1);

[min_value_4_1, index_4_1] = min(rsm_pf4_1); % 找到最小值及其索引
pos_AIPF_1 = pos_AIPF_1(:,:,index_4_1);

%平均匹配误差指的是，n次实验的均方根误差的均值；稳定性指n次实验的标准差
disp("*7*分析：****************************************************************")
disp("                        PF       AOFA    EKPF      AIPF            ")
fprintf("轨迹1：平均计算时间 ： %f       %f      %f       %f\n",mean_tm_1(1),mean_tm_1(2),mean_tm_1(3),mean_tm_1(4));

disp("                                                              ")
fprintf("轨迹1：平均匹配误差 ： %f       %f      %f       %f\n",mean_rsm_pf_1(1),mean_rsm_pf_1(2),mean_rsm_pf_1(3),mean_rsm_pf_1(4));

disp("                                                              ")

fprintf("轨迹1：平均误差均值： %f       %f      %f       %f\n",mean_mean_pf_1(1),mean_mean_pf_1(2),mean_mean_pf_1(3),mean_mean_pf_1(4));

%fprintf("轨迹2：阱的大小： %f \n",rmse_matrix(const_j*10+1,const_i*10+1));
%绘图展示四种算法num次计算得到的轨迹与真实轨迹的误差值图
plotResult(rsm_pf1_1,rsm_pf2_1,rsm_pf3_1,rsm_pf4_1);

%% 绘制轨迹1粒子滤波算法的结果
figure;
%ax = geoaxes;
worldmap([bottom_lat top_lat], [left_lon right_lon]);
% 在地理坐标轴上绘制地磁图
geoshow(X, Y, magnetic_mod, 'DisplayType', 'texturemap');
view(2);
colormap("jet");

% 在地理坐标轴上绘制真实轨迹和散点
h1 = geoshow(true_pos1(:,1), true_pos1(:,2), 'DisplayType', 'line', 'Color', 'k', 'LineWidth', 2);
h2 = geoshow(pos_PF_1(:,1),pos_PF_1(:,2),'DisplayType', 'line', 'Color', 'g','LineWidth', 2);
h3 = geoshow(pos_AOFA_1(:,1),pos_AOFA_1(:,2),'DisplayType', 'line', 'Color', 'b','LineWidth', 2);
h4 = geoshow(pos_EKPF_1(:,1),pos_EKPF_1(:,2),'DisplayType', 'line', 'Color', 'y','LineWidth', 2);
h5 = geoshow(pos_AIPF_1(:,1), pos_AIPF_1(:,2), 'DisplayType', 'line', 'Color', 'r', 'LineWidth', 2);
h6 = geoshow(ins_pos(:,1), ins_pos(:,2), 'DisplayType', 'line', 'Color', 'm', 'LineWidth', 2);
legend([h1,h2,h3,h4,h5,h6],'True', 'PF', 'AOFA', 'EKPF', 'AIPF','INS');
geoshow('landareas.shp', 'FaceAlpha', 0);
hold on; % 保持图形
xlabel({'经度','(a) 磁地图模值纹理图'},'FontSize',12,'Fontname', '黑体','FontWeight','bold');
ylabel('纬度');
hold off; % 取消保持图形

