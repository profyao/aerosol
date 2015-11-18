x = [1,2,4,8,16,32];
t_local = [288,133,69,37,27,23];
t_cluster = [117,70,39,25,22,20];
t_local_n = t_local(1)./t_local;
t_cluster_n = t_cluster(1)./t_cluster;
figure
loglog(x,t_local_n,'-o','LineWidth',2),hold on
loglog(x,t_cluster_n,'-o','LineWidth',2)
xlim([1,50]),ylim([1,50])
set(gca,'xtick',[1,10,50],'xticklabel',{'10^0','10^1','50'})
set(gca,'ytick',[1,10,50],'yticklabel',{'10^0','10^1','50'})
legend({'Local Mode','Cluster Mode'})
set(gca,'FontSize',18)
grid on
xlabel('Number of Cores')
ylabel('Normalized Speed')
title('Normalized Speed v.s. Parallelism')