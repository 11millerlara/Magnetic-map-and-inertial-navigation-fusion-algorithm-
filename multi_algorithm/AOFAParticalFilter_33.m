function pos=AOFAParticalFilter_33(n,x,y,sl,deg,p,q,mag,dataset)
    sizedataset = size(mag,2);
    P=3^2*eye(sizedataset);
w=1/n*ones(n,1);
x_p=x*ones(n,1);
y_p=y*ones(n,1);
x_n=zeros(n,1);
y_n=zeros(n,1);
pBefores=zeros(n,2);
step=length(sl);
pos=zeros(step+1,2);
pos(1,:)=[x,y];
% figure
% hold on
% ylim([-2,3])
for i=1:step
    if(q(i)<0.1)
        q(i)=0.1;
    end
    for j=1:n
        sl_t=sl(i)+normrnd(0,p(i));
        deg_t=deg(i)+normrnd(0,q(i));
        pBefores(i,:)=[x_p(j),y_p(j)];
        x_p(j)=x_p(j)+sl_t*sin(deg_t);
        y_p(j)=y_p(j)+sl_t*cos(deg_t);
        index=Find_33(dataset,x_p(j),y_p(j));
        if(index(1)==-1||index(2)==-1)
            w(j)=0;
        else
            w(j)=exp(-0.5*(mag(i,:)-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))*inv(P)*(mag(i,:)-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))');
            if(isnan(w(j)))
                w(j)=0;
            end
        end
        
    end
    
     %萤火虫算法过程
        particles=[x_p,y_p];
        
        maxIndex=1;
        maxWeight=-1;
        for j=1:n
            if(w(j,1)>maxWeight)
               maxWeight=w(j,1);
               maxIndex=j;
            end
        end
        ua=0;
        for j=1:n
            ua=ua+(w(j,1)-w(maxIndex,1))^2;
        end
        ua=ua/n;
        
        wavg=0;
        for j=1:n
            wavg=wavg+w(j,1);
        end
        wavg=wavg/n;
        %maxLength=-1;
        %maxLengtnIndex=1;
        maxX=-1;
        maxY=-1;
        for j=1:n
            %dist=sqrt((particles(i,1)-particles(maxIndex,1))^2+(particles(i,2)-particles(maxIndex,2))^2);
            distX=abs(particles(j,1)-particles(maxIndex,1));
            distY=abs(particles(j,2)-particles(maxIndex,2));
            if(distX>maxX)
                maxX=distX;
            end
            if(distY>maxY)
                maxY=distY;
            end
        end
        
        for j=1:n
           if(j~=maxIndex)
             ww=w(j,1);
             x1=particles(j,1);
             y1=particles(j,2);
             dist=power(x1-particles(maxIndex,1),2)+power(y1-particles(maxIndex,2),2);
             beta=(0.1+1/(sqrt(i)+1))*ww*exp(-1*dist);                        %1.5和5可能比0.5好一些  改进萤火虫算法优化粒子滤波的信号源定位  吸引度跟随迭代次数改变  K的平方和根号K都行 根号K少见一些

             %beta=exp(-dist);
%              x1=x1+beta*(particles(maxIndex,1)-x1)+0.5*(rand([1,1])-0.5);%0.05比0.25好一些
%              y1=y1+beta*(particles(maxIndex,2)-y1)+0.5*(rand([1,1])-0.5);

             %自适应步长
             ex=((particles(maxIndex,1)-particles(j,1))*(particles(j,1)-pBefores(j,1)))/(abs((particles(maxIndex,1)-particles(j,1))*(particles(j,1)-pBefores(j,1))));
             ey=((particles(maxIndex,2)-particles(j,2))*(particles(j,2)-pBefores(j,2)))/(abs((particles(maxIndex,2)-particles(j,2))*(particles(j,2)-pBefores(j,2))));
             %
%              dx=abs(particles(maxIndex,1)-particles(j,1))/maxX;
%              dy=abs(particles(maxIndex,2)-particles(j,2))/maxY;              %这个是论文Multi-target tracking method based on improved firefly  algorithm optimized particle filter中的方法
             
             dx=1/(1+exp((w(j,1)-wavg)/(w(maxIndex,1)-wavg)));%这个论文Firefly Algorithm With Disturbance-Factor-Based Particle Filter for Seismic Random Noise Attenuation中的方法
             dy=1/(1+exp((w(j,1)-wavg)/(w(maxIndex,1)-wavg)));
             
             if(ex>0)
                 x1=x1+beta*(particles(maxIndex,1)-x1)+(rand([1,1])-0.5)*(dx);
             else
                 x1=x1+beta*(particles(maxIndex,1)-x1)-(rand([1,1])-0.5)*(dx)*0.5;
             end
             if(ey>0)
                 y1=y1+beta*(particles(maxIndex,2)-y1)+(rand([1,1])-0.5)*(dy);
             else
                 y1=y1+beta*(particles(maxIndex,2)-y1)-(rand([1,1])-0.5)*(dy)*0.5;
             end
             
%              x1=x1+beta*(particles(maxIndex,1)-x1)+0.5*(rand([1,1])-0.5);
%              y1=y1+beta*(particles(maxIndex,2)-y1)+0.5*(rand([1,1])-0.5);
             particles(j,1)=x1;
             particles(j,2)=y1;
           end
        end
        
        ub=ua;
        
        for j=1:n
            index=Find_33(dataset,particles(j,1),particles(j,2));
            if isnan(index(1))
                index(1) = -1;
                index(2) = -1;
            end
            if(index(1)==-1||index(2)==-1)
                w(j)=0;
            else
                w(j)=exp(-0.5*(mag(i,:)-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))*inv(P)*(mag(i,:)-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))');
                if(isnan(w(j)))
                    w(j)=0;
                end
            end       
        end
    x_p=particles(:,1);   
    y_p=particles(:,2);
    %归一化
	w=w/sum(w,1);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     pos(i+1,1)=x_p.'*w;
%     pos(i+1,2)=y_p.'*w;
%     pos(i+1,:)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    c=zeros(n,1);
    c(1)=w(1);
    for j=2:n
        c(j)=c(j-1)+w(j);
    end
    for j=1:n
        a=unifrnd(0,1);
        for k=1:n
            if(a<c(k))
                x_n(j)=x_p(k); 
                y_n(j)=y_p(k);
                break;
            end
        end
    end
    w=1/n*ones(n,1);
    x_p=x_n;
    y_p=y_n;
%     plot(x_n,y_n,'o')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    pos(i+1,1)=x_p.'*w;
    pos(i+1,2)=y_p.'*w;
end

