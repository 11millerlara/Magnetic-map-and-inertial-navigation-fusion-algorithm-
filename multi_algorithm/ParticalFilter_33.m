function pos=ParticalFilter_33(n,x,y,sl,deg,p,q,mag,dataset)
    sizedataset = size(mag,2);
    P=3^2*eye(sizedataset);
    w=1/n*ones(n,1);
    x_p=x*ones(n,1);
    y_p=y*ones(n,1);
    x_n=zeros(n,1);
    y_n=zeros(n,1);
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
        %归一化
	    w=w/sum(w,1);
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
end
