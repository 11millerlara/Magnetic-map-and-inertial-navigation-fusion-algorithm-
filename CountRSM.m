function [res,MEAN]=CountRSM(a,b)
    via=sum((a-b).^2,2);
    res=sqrt(sum(via,1)/length(via));
    r=sqrt(via);
    MIN=min(r);
    MAX=max(r);
    MEAN=mean(r);
    test=sprintf("min=%f,max=%f,mean=%f,rsm=%f",MIN,MAX,MEAN,res);
    disp(test)
end
