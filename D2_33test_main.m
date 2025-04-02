% 磁地图理论参考仿真器
clear all; clc; close all;

load('./BPIT/fx.mat')
load('./BPIT/xfx.mat')
addpath('./multi_algorithm')
addpath('./cal_pos')
disp("*1*输入初始数据******************************************************")
set(0,'defaultfigurecolor','w')%显示背景设置为白色
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('./magnetic_map_data/Data.mat'); %读取参考轨迹
load('./magnetic_map_data/magnetic_mod.mat');
magnetic_mod = Z0;
min_lat = min(Data(:,1));
min_lon = min(Data(:,2));
%真实轨迹矩阵_闭合曲线
[true_pos1,init_pos1_lat,init_pos1_lon,const_i,const_j] = cal_33_pos2_10();
%[true_pos1,init_pos1_lat,init_pos1_lon] = cal_33_pos2_3();
number = size(true_pos1,1);
magnetic_pos1=ones(number,1);
for i=1:number
    lat_point = round(true_pos1(i,1)-min_lat)*10+1;    %读取轨迹坐标的在磁地图中的行数和列数
    lon_point = round(true_pos1(i,2)-min_lon)*10+1;
    magnetic_pos1(i,1)=magnetic_mod(lat_point,lon_point);
end

%% 模拟导航位置初始化误差
mmp =0.1;
noise_lat = 0;
noise_lon = 0;
init_pos1_lat = true_pos1(1,1) + mmp* noise_lat;
init_pos1_lon = true_pos1(1,2) + mmp * noise_lon;

N = size(magnetic_mod,1);
M = size(magnetic_mod,2);
dataset = zeros(N,M,1);
dataset(:,:,1) = magnetic_mod;


%% 生成轨迹数据近似的步长数据矩阵和方向数据矩阵——曲折闭合曲线
M =size(true_pos1,1);
sl1=zeros(M-1,1);
deg1=zeros(M-1,1);
for i=1:M-1
    sl1(i)=sqrt((true_pos1(i+1,1)-true_pos1(i,1)).^2+(true_pos1(i+1,2)-true_pos1(i,2)).^2);
    y=true_pos1(i+1,1)-true_pos1(i,1);
    x=true_pos1(i+1,2)-true_pos1(i,2);
    theta = atan2d(y, x); % 两点之间的方向角，单位为度
    deg1(i) = theta; % 将方向角（度）存储到 deg1 中
end

% 设定步长和航向的初始值
sl1_IMU = sl1; % 初始步长
deg1_IMU = deg1; % 初始航向角

% 参数设置
step_noise_mean = 0.02; % 步长噪声的均值，单位可以为米
step_noise_std = 0.05; % 步长噪声的标准差

heading_noise_mean = 0.3; % 航向角噪声的均值，单位为度
heading_noise_std = 0.5; % 航向角噪声的标准差

step_drift_rate = 0.001; % 步长漂移率，单位为米/步
heading_drift_rate = 0.05; % 航向角漂移率，单位为度/步

% 添加噪声和漂移
for i = 2:length(sl1)
    % 添加非零均值的高斯噪声
    sl1_IMU(i) = sl1(i) + step_noise_mean + step_noise_std * randn();
    deg1_IMU(i) = deg1(i) + heading_noise_mean + heading_noise_std * randn();
    
    % 累计漂移
    sl1_IMU(i) = sl1_IMU(i) + step_drift_rate * i;
    deg1_IMU(i) = deg1_IMU(i) + heading_drift_rate * i;
    
    % 添加随机游走
    sl1_IMU(i) = sl1_IMU(i) + 0.01 * (rand() - 0.5);
    deg1_IMU(i) = deg1_IMU(i) + 0.1 * (rand() - 0.5);
end

% 结果：sl1_IMU 和 deg1_IMU 为模拟的IMU步长和航向输出

sl1 = sl1_IMU;
deg1 = deg1_IMU;
% 结果：sl1_IMU 和 deg1_IMU 为模拟的IMU步长和航向输出

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



deg1 = deg1*pi/180;
%% 利用粒子滤波算法计算地磁匹配结果
number_mag = 1; 
disp("*6*分别利用4种算法计算地磁匹配结果数据***************************************")
disp("*6.1*计算第一个轨迹的地磁匹配结果数据")
[pos_PF_1,pos_AOFA_1,pos_EKF_1,pos_AIPF_1,rsm_pf1_1,rsm_pf2_1,rsm_pf3_1,rsm_pf4_1,mean_tm_1,mean_rsm_pf_1,mean_mean_pf_1]= ...
    prediction_pos1(init_pos1_lon,init_pos1_lat,sl1,deg1,magnetic_pos1(:,[1:number_mag]),dataset(:,:,[1:number_mag]),fx,xfx,true_pos1,min_lat,min_lon);

%        轨迹1             
stdpf1_1 = std(rsm_pf1_1);  
stdpf2_1 = std(rsm_pf2_1);  
stdpf3_1 = std(rsm_pf3_1);  
stdpf4_1 = std(rsm_pf4_1); 

[min_value_1_1, index_1_1] = min(rsm_pf1_1); % 找到最小值及其索引
pos_PF_1 = pos_PF_1(:,:,index_1_1);


[min_value_2_1, index_2_1] = min(rsm_pf2_1); % 找到最小值及其索引
pos_AOFA_1 = pos_AOFA_1(:,:,index_2_1);


[min_value_3_1, index_3_1] = min(rsm_pf3_1); % 找到最小值及其索引
pos_EKF_1 = pos_EKF_1(:,:,index_3_1);


[min_value_4_1, index_4_1] = min(rsm_pf4_1); % 找到最小值及其索引
pos_AIPF_1 = pos_AIPF_1(:,:,index_4_1);

[rsm5,MEAN5]=CountRSM(ins_pos(1:end,:)/100,true_pos1(1:end,:)/100);

%平均匹配误差指的是，n次实验的均方根误差的均值；稳定性指n次实验的标准差
disp("*7*分析：****************************************************************")
disp("                        PF       AOFA    EKPF      AIPF          INS           ")
fprintf("轨迹1：平均计算时间 ： %f       %f      %f       %f\n",mean_tm_1(1),mean_tm_1(2),mean_tm_1(3),mean_tm_1(4));

disp("                                                              ")
fprintf("轨迹1：平均匹配误差 ： %f       %f      %f       %f            %f\n",mean_rsm_pf_1(1),mean_rsm_pf_1(2),mean_rsm_pf_1(3),mean_rsm_pf_1(4),rsm5);

disp("                                                              ")
%fprintf("轨迹1：稳定性 ：      %f       %f      %f       %f\n",stdpf1_1,stdpf2_1,stdpf3_1,stdpf4_1);
%fprintf("轨迹2：稳定性 ：      %f       %f      %f       %f\n",stdpf1_2,stdpf2_2,stdpf3_2,stdpf4_2);
%fprintf("轨迹3：稳定性：       %f       %f      %f       %f\n",stdpf1_3,stdpf2_3,stdpf3_3,stdpf4_3);

fprintf("轨迹1：平均误差均值： %f       %f      %f       %f             %f\n",mean_mean_pf_1(1),mean_mean_pf_1(2),mean_mean_pf_1(3),mean_mean_pf_1(4),MEAN5);


%fprintf("轨迹2：阱的大小： %f \n",rmse_matrix(const_j*10+1,const_i*10+1));
%绘图展示四种算法num次计算得到的轨迹与真实轨迹的误差值图

plotResult(rsm_pf1_1,rsm_pf2_1,rsm_pf3_1,rsm_pf4_1);
hold off;



set(0,'defaultfigurecolor','w')%显示背景设置为白色
figure(3)
plot(Data(:,2)/100,Data(:,1)/100,'r*-','LineWidth',0.01);%横坐标为归一化经度，纵坐标为归一化纬度

figure (4)
plot3(Data(:,2)/100,Data(:,1)/100,Data(:,3),'r*-','LineWidth',0.01);%横坐标为归一化经度，纵坐标为归一化纬度，高度坐标为非归一化高度（飞机距海平面的相对高度，非海拔高度）


figure (5)
plot3(Data(:,2)/100,Data(:,1)/100,Data(:,4),'r*-','LineWidth',0.01);%横坐标为归一化经度，纵坐标为归一化纬度，高度坐标为补偿掉主磁场后的局部失真磁场（nT）


figure(6)
subplot(1,2,1)
minx0=min(Data(:,1));
maxx0=max(Data(:,1));
miny0=min(Data(:,2));
maxy0=max(Data(:,2));
% 计算网格点的数量，精度为0.1
x_points = minx0:0.1:maxx0;
y_points = miny0:0.1:maxy0;
[Y0, X0] = meshgrid(y_points, x_points);
Z0=griddata(Data(:,2),Data(:,1),Data(:,4),Y0,X0,'v4');
mesh(Y0/100,X0/100,Z0)
hold on
plot3(Data(:,2)/100,Data(:,1)/100,Data(:,4),'r.','MarkerFaceColor','r')
view(3)
xlabel('Normalized longitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Normalized latitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
zlabel('Magnetic anomaly(nT)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
colorbar

subplot(1,2,2)
mesh(Y0/100,X0/100,Z0)
hold on
plot3(Data(:,2)/100,Data(:,1)/100,Data(:,4),'r*-','MarkerFaceColor','r')
view(2)
xlabel('Normalized longitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Normalized latitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
zlabel('Magnetic anomaly(nT)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
colorbar


N = size(true_pos1,1);
Data1 = zeros(N,4);
Data1(:,1:2) = true_pos1;
Data1(:,4) = magnetic_pos1;
%% 创建磁地图的三维表面
figure(7);
mesh(Y0/100,X0/100,Z0)
hold on
% 绘制真实轨迹路径
h1=plot3(Data1(:,2)/100,Data1(:,1)/100,Data1(:,4),'k*-','MarkerFaceColor','k');
% 绘制AIPF轨迹路径
h2=plot3(pos_AIPF_1(:,2)/100, pos_AIPF_1(:,1)/100, Data1(:,4),'r-o', 'LineWidth', 0.2);
% 绘制PF轨迹路径
h3=plot3( pos_PF_1(:,2)/100, pos_PF_1(:,1)/100,Data1(:,4), 'g-*', 'LineWidth', 0.2); 
% 绘制AOFA轨迹路径
h4=plot3(pos_AOFA_1(:,2)/100, pos_AOFA_1(:,1)/100, Data1(:,4),'b-s', 'LineWidth', 0.2); 
% 绘制EKF轨迹路径
h5=plot3(pos_EKF_1(:,2)/100, pos_EKF_1(:,1)/100, Data1(:,4),'m-d', 'LineWidth', 0.2); 

h6=plot3(ins_pos(:,2)/100, ins_pos(:,1)/100, Data1(:,4),'c-p', 'LineWidth', 0.2); 
% 设置图例
%legend([h1,h2,h3,h4,h5],'True', 'AIPF', 'PF', 'AOFA', 'EKF');
legend([h1,h2,h3,h4,h5,h6],'真实轨迹', 'AIPF轨迹', 'PF轨迹', 'AOFA轨迹', 'EKPF轨迹','INS轨迹');
%legend([h1,h2,h3,h4,h6],'真实轨迹', 'AIPF轨迹', 'PF轨迹', 'AOFA轨迹','ins轨迹');
view(3)
xlabel('Normalized longitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Normalized latitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
zlabel('Magnetic anomaly(nT)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
colorbar


figure(8);

hold on

% 绘制真实轨迹路径
h1=plot(Data1(:,2)/100,Data1(:,1)/100,'k*-','MarkerFaceColor','k');
% 绘制AIPF轨迹路径
h2=plot(pos_AIPF_1(:,2)/100, pos_AIPF_1(:,1)/100,'r-o', 'LineWidth', 0.2);
% 绘制PF轨迹路径
h3=plot( pos_PF_1(:,2)/100, pos_PF_1(:,1)/100, 'g-*', 'LineWidth', 0.2); 
% 绘制AOFA轨迹路径
h4=plot(pos_AOFA_1(:,2)/100, pos_AOFA_1(:,1)/100,'b-s', 'LineWidth', 0.2); 
% 绘制EKF轨迹路径
h5=plot(pos_EKF_1(:,2)/100, pos_EKF_1(:,1)/100,'m-d', 'LineWidth', 0.2); 

h6=plot(ins_pos(:,2)/100, ins_pos(:,1)/100,'c-p', 'LineWidth', 0.2); 
% hold on
% 映射磁地图，并设置透明度
surf(Y0/100, X0/100, Z0, magnetic_mod, 'EdgeColor', 'none');
uistack(h1, 'top');
uistack(h2, 'top');
uistack(h3, 'top');
uistack(h4, 'top');
uistack(h5, 'top');
uistack(h6, 'top');
set(gca, 'SortMethod', 'childorder');  % 根据添加顺序进行叠加显示
%mesh(Y0,X0,Z0)
% 设置图例
%legend([h1,h2,h3,h4,h5],'True', 'AIPF', 'PF', 'AOFA', 'EKF');
legend([h1,h2,h3,h4,h5,h6],'True', 'AIPF', 'PF', 'AOFA', 'EKPF','INS');
%legend([h1,h2,h3,h4,h6],'真实轨迹', 'AIPF轨迹', 'PF轨迹', 'AOFA轨迹','ins轨迹');

xlabel('Normalized longitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Normalized latitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
zlabel('Magnetic anomaly(nT)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');

% 计算每个算法估计位置与真实位置之间的误差
error_PF = sqrt(sum((pos_PF_1/100 - true_pos1/100).^2, 2));
error_AOFA = sqrt(sum((pos_AOFA_1/100 - true_pos1/100).^2, 2));
error_EKF = sqrt(sum((pos_EKF_1/100 - true_pos1/100).^2, 2));
error_AIPF = sqrt(sum((pos_AIPF_1/100 - true_pos1/100).^2, 2));

% 绘制CDF图
figure(9);
hold on;

% 绘制各算法误差的CDF
cdfplot(error_PF);
cdfplot(error_AOFA);
cdfplot(error_EKF);
cdfplot(error_AIPF);

% 设置图例
legend('PF', 'AOFA', 'EKF','AIPF', 'Location', 'best');

% 设置轴标签和标题
xlabel('Error','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Cumulative Distribution Function (CDF)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
title('','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');

% 设置网格
grid on;

hold off;








