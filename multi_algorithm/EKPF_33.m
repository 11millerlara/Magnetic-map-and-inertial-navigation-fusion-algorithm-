function [pos]=EKPF_33(n,y,x,sl,deg,p,q,mag,dataset)
    sizedataset = size(mag,2);
    P1=3^2*eye(sizedataset);
    step=length(sl);
    pos=zeros(step+1,2);
    pos(1,:)=[x,y];
    %% 初始化参数
%     Pplus=diag([0,0,p(1)^2,q(1)^2]);%初始误差
    Pplus=diag([0,0]);%初始误差

%     C=[1,0,0,0; %观测矩阵
%        0,1,0,0];
      C=[1,0; %观测矩阵
         0,1,];

    R=diag([0.5^2,0.5^2]);%观测噪声协方差

    for i=1:step
        x=pos(i,1);
        y=pos(i,2);
        l=sl(i);
        theta=deg(i);
%         X=[x,y,l,theta].';%变量
%         A=[1,0,cos(theta),0;%系统矩阵
%            0,1,sin(theta),0;
%            0,0,1         ,0;
%            0,0,0         ,1];
%         P=diag([0.2^2,0.2^2,p(i)^2,q(i)^2]);%X变量误差
%         Q=A*P*A.';%系统协方差计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_e=x+l*cos(theta);
        y_e=y+l*sin(theta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X=[x_e,y_e].';
        A=[1,0,cos(theta),0;%系统矩阵
           0,1,sin(theta),0;
           0,0,1         ,0;
           0,0,0         ,1];
        P=diag([0.2^2,0.2^2,p(i)^2,q(i)^2]);%X变量误差
        Q=A*P*A.';%系统协方差计算
        Q=Q(1:2,1:2);
        A=A(1:2,1:2);
        Y=KNN(x_e,y_e,mag(i,:),dataset);
        Y=Y.';
        %% 预测步
        Xminus=A*X;
        Pminus=A*Pplus*A.'+Q;
        %% 更新步
        K=Pminus*C.'*inv(C*Pminus*C.'+R);
        Xplus=Xminus+K*(Y-C*Xminus);
%         Pplus=(eye(4)-K*C)*Pminus;
        Pplus=(eye(2)-K*C)*Pminus;
        [U,S,V]=svd(Pplus);
        phi=atan(U(2,1)/U(1,1));
        W=sqrt(S(1,1));
        H=sqrt(S(2,2));
        Rbn=[cos(phi),sin(phi);
            -sin(phi),cos(phi)];
        particle=[Xplus(1)*ones(n,1),Xplus(2)*ones(n,1)];
        w=zeros(n,1);
        for j=1:n
            param=unifrnd(-1,1,[2,1]);
            param(1)=param(1)*W;
            param(2)=param(2)*H;
            param=Rbn*param;
            particle(j,:)=particle(j,:)+param.';
            index=Find_33(dataset,particle(j,2),particle(j,1));
            if isnan(index(1))
                index(1) = -1;
                index(2) = -1;
            end
            if(index(1)==-1||index(2)==-1)
                w(j)=0;
            else
                w(j)=exp(-0.5*(mag(i,:)-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))*inv(P1)*(mag(i,:)-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))');
                if(isnan(w(j)))
                    w(j)=0;
                end
            end
        end
        %归一化
        w=w/sum(w,1);
        c=zeros(n,1);
        x_n=zeros(n,1);
        y_n=zeros(n,1);
        c(1)=w(1);
        for j=2:n
            c(j)=c(j-1)+w(j);
        end
        for j=1:n
            a=unifrnd(0,1);
            for k=1:n
                if(a<c(k))
                    x_n(j)=particle(k,1);
                    y_n(j)=particle(k,2);
                    break;
                end
            end
        end
        w=1/n*ones(n,1);
        pos(i+1,1)=x_n.'*w;
        pos(i+1,2)=y_n.'*w; 
%        pos(i+1,:)=Xplus.';
    end
    pos(:,[1 2]) = pos(:,[2 1]); % 交换第一列和第二列，第一列为纬度数据，第二列为经度数据
end
function [result]=KNN(x,y,mag,dataset)
    %Nnum=100;
    sizedataset = size(mag,2);
%     pos=zeros(Nnum,2);
    M=size(dataset,1);
    N=size(dataset,2);
    num=M*N;
    distance=zeros(num,3);
    cnt=0;
    radius=1;
    for i=1:M
        for j=1:N
            tmp(2) = (i-1)/10 +34.2;
            tmp(1) = (j-1)/10 +90;
            if(tmp(1)<x+radius&&tmp(1)>x-radius&&tmp(2)<y+radius&&tmp(2)>y-radius)
                cnt=cnt+1;
    %             
                distance(cnt,1)=sum((mag-reshape(dataset(i,j,1:sizedataset),1,sizedataset)).^2);
                distance(cnt,2)=tmp(1);
                distance(cnt,3)=tmp(2);
            end
        end
    end
    distance=distance(1:cnt,:);
    distance=sortrows(distance,1);
    pos=distance(1:end,2:3);
    result=mean(pos);
end



