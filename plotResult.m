function plotResult(rsm1,rsm2,rsm3,rsm4)
    figure 
    plot(rsm1,'ok:','LineWidth',1,'MarkerSize',8, ... 
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0, 255, 0]/255)
    hold on
    plot(rsm2,'ok:','LineWidth',1,'MarkerSize',8, ... 
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0, 0, 255]/255)
    hold on
    plot(rsm3,'ok:','LineWidth',1,'MarkerSize',8, ... 
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [255, 255, 0]/255)
    hold on
    plot(rsm4,'ok:','LineWidth',1,'MarkerSize',8, ... 
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [255, 0, 0]/255)
    legend('PF','AOFA','EKPF','AIPF')
    xlabel('Number of experiments','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
    ylabel('Root mean square error','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');

end