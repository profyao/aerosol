figure
for i=1:8
    subplot(4,2,i)
    plot(acf(squeeze(sample.theta(i,1847,:)),1000))
    title(strcat('Component Num:',num2str(const.Component_Particle(i))))
end

figure
for i=1:8
    subplot(4,2,i)
    plot(cumsum(squeeze(sample.theta(i,1847,:)))./[1:1001]')
    title(strcat('Component Num:',num2str(const.Component_Particle(i))))
end

%[sample,reg] = load_cache('2011.07.29',15,61764,59,4400,'sample','reg','CD-random');
[sample,reg] = load_cache('2011.06.02',16,60934,59,4400,'sample','reg','MCMC');
%[sample,reg] = load_cache('2011.06.04',14,60963,59,4400,'sample','reg','CD-random');

hFig = figure;
set(hFig, 'Position', [10 10 800 600])

ang_ind = 1:36;
sample_ind = 1;
total = sample.surf(ang_ind,sample_ind)+sample.atm_path(ang_ind,sample_ind)+sample.resid(ang_ind,sample_ind);

plot(sample.surf(ang_ind,sample_ind),'Color',const.cols(1,:))
hold on
plot(sample.atm_path(ang_ind,sample_ind),'Color',const.cols(3,:))
hold on
plot(total,'Color',const.cols(2,:))
xlabel('Camera Index')
ylabel('Reflectance')
title('Angular Shape of Reflectance')
set(gca,'FontSize',18)
%ylim([-0.05,0.35])

figure,plot_1d(4400,mean(sample.tau(:,200:200:end),2),x,y,jet,const,[0.02,0.2])
coord = [22,97;23,97;19,93;21,82;30,73];
[~,~,pid] = intersect(coord,[x,y],'rows');

figure
for p = pid'
    plot(cumsum(sample.tau(p,:))./[1:1001]),hold on
end
xlabel('Iteration Number')
ylabel('Cumulative Mean of \tau_p')
set(gca,'FontSize',18)

figure
for p = pid'
    plot(acf(sample.tau(p,:)',1000)),hold on
end
xlabel('Iteration Lag')
ylabel('Correlation')
set(gca,'FontSize',18)


figure('position',[300 300 800 600])
plot(tau,sample_mle_50.tau,'o','Color',const.cols(1,:)),hold on
plot(tau,mean(sample_mcmc_50.tau(:,100:20:end),2),'o','Color',const.cols(3,:))
plot(tau,sample_map_50_kappa.tau,'o','Color',const.cols(2,:)),hold on
%plot(tau,mean(sample_mcmc_00.tau(:,100:20:end),2),'o','Color',const.cols(3,:)),hold on
%plot(tau,sample_map_00.tau,'o','Color',const.cols(1,:)),hold on
%plot(tau,sample_map_00_kappa.tau,'o','Color',const.cols(2,:))
cy = [0.05,0.55];
cx = [0.05,0.55];
xlim(cx),ylim(cy)
line('XData', cx, 'YData', cy, 'LineStyle', '-','LineWidth', 1, 'Color','k')
xlabel('True AOD in Simulation')
ylabel('Retrieved AOD')
set(gca,'FontSize',18)
%legend({'MCMC','MAP+Dirichlet','MAP+uniform'},'Location','northwest')
legend({'MISR Operational','MCMC+Dirichlet','MAP+uniform'},'Location','northwest')
%legend({'MAP+GMRF+Dirichlet','MLE','MCMC+GMRF+Dirichlet','MAP+GMRF'},'Location','northwest')
%legend({'MAP+GMRF+Dirichlet','MAP+GMRF'},'Location','northwest')

figure('position',[300 300 800 600])
t = zeros(const.Component_Num,2);
t(:,1) = mean(sample_map_00.theta,2);
t(:,2) = mean(sample_map_00_kappa.theta,2);
t(:,3) = mean(squeeze(mean(sample_mcmc_00.theta(:,:,100:20:end),3)),2);
bar(t)
legend({'MAP+Dirichlet','MAP+uniform','MCMC'},'Location','northwest')
hline(0.125,'--r')
set(gca,'FontSize',18)
% s1 = scatter(tau,sample_mle_00.tau,'MarkerEdgeColor',const.cols(1,:),'MarkerFaceColor',const.cols(1,:));hold on
% s2 = scatter(tau,sample_map_00.tau,'MarkerEdgeColor',const.cols(2,:),'MarkerFaceColor',const.cols(2,:));hold on
% s3 = scatter(tau,mean(sample_mcmc_00.tau(:,100:20:end),2),'MarkerEdgeColor',const.cols(3,:),'MarkerFaceColor',const.cols(3,:));
% set(s1,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.3)
% set(s2,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.3)
% set(s3,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.3)
% legend({'MLE','MAP+GMRF+Dirichlet','MCMC+GMRF+Dirichlet'},'Location','northwest')



figure,bar([60,6280,600])
set(gca,'XTickLabel',{'MAP','MCMC','MISR'},'FontSize',18)
ylabel('Computation Time (sec)')