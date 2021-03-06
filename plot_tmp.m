cores = [1,2,4,8,16,32];
t_local = [288,133,69,37,27,23];
t_cluster = [117,70,39,25,22,20];
t_local_n = t_local(1)./t_local;
t_cluster_n = t_cluster(1)./t_cluster;
figure
loglog(cores,t_local_n,'-o','LineWidth',2),hold on
loglog(cores,t_cluster_n,'-o','LineWidth',2)
xlim([1,50]),ylim([1,50])
set(gca,'xtick',[1,10,50],'xticklabel',{'10^0','10^1','50'})
set(gca,'ytick',[1,10,50],'yticklabel',{'10^0','10^1','50'})
legend({'Local Mode','Cluster Mode'})
set(gca,'FontSize',18)
grid on
xlabel('Number of Cores')
ylabel('Normalized Speed')
title('Normalized Speed v.s. Parallelism')

figure,
plot(tau,sample_rls.tau,'o'),hold on
plot(tau,mean(sample_mcmc.tau(:,50:10:end),2),'o'),hold on,
plot(tau,sample_cd_prior.tau,'o')
xlim([0.18,0.4]),ylim([0.18,0.4])
refline(1,0)
legend({'Random Local Search','MCMC','Coordinate Descent'})
set(gca,'FontSize',18)
xlabel('True AOD in Simulation')
ylabel('Retrieved AOD')


figure
subplot(3,1,1),plot_1d(4400,tau,x,y,jet,const)
subplot(3,1,2),plot_1d(4400,sample.tau,x,y,jet,const)
subplot(3,1,3),plot_1d(4400,sample.tau-tau,x,y,jet,const)

% figure
% for p = 1:50:reg.num_reg_used
%     chisq = zeros(100,1);
%     tau_grid = linspace(0,0.5,100);
%     xp = x(p);
%     yp = y(p);
%     thetap = sample_cd.theta(:,p);
%     for i = 1:100
%         chisq(i) = get_chisq(r,tau_grid(i),thetap,sample_cd.sigmasq(:,end),xp,yp,const,ExtCroSect,CompSSA,smart,reg,reg_sim);
%     end
%     plot(tau_grid,chisq,'k'),hold on
% end
% set(gca,'FontSize',18)
% xlabel('AOD')
% ylabel('\chi^2_p','rot',0)
% title('Convexity near Optima')
% 
% figure
% for p = 1:50:reg.num_reg_used
%     chisq = zeros(100,1);
%     theta1_grid = linspace(0,1,100);
%     theta2_grid = linspace(1,0,100);
%     theta_grid = [theta1_grid;theta2_grid];
%     xp = x(p);
%     yp = y(p);
%     taup = tau(p);
%     for i = 1:100
%         chisq(i) = get_chisq(r,taup,theta_grid(:,i),sample_cd.sigmasq(:,end),xp,yp,const,ExtCroSect,CompSSA,smart,reg,reg_sim);
%     end
%     plot(theta_grid(1,:),chisq,'k'),hold on
% end
% set(gca,'FontSize',18)
% xlabel('Component 2 Percentage')
% ylabel('\chi^2_p','rot',0)
% title('Convexity near Optima')

[chisq,chisq_misr,tau_misr,theta_misr] = retri_aod_misr(theta_misr_grid,r,sample,x,y,const,ExtCroSect,CompSSA,smart,reg);

figure,plot(chisq,chisq_misr,'o'),refline(1,0)
xlabel('MAP \chi^2_p')
ylabel('MISR \chi^2_p')
set(gca,'FontSize',18)

theta1_grid = linspace(0,1,100);
theta2_grid = linspace(1,0,100);
%theta2_grid = kron((1 - theta1_grid),thetap(2:end));
theta_grid = [theta1_grid;theta2_grid];
tau_grid = linspace(0,0.5,100);
xp = x(p);
yp = y(p);
[X,Y] = meshgrid(theta1_grid,tau_grid);
Z = zeros(100,100);
for ii = 1:100
    for jj = 1:100
        chisq = 0.5 * get_chisq(r,tau_grid(ii),theta_grid(:,jj),sample.sigmasq(:,end),xp,yp,const,ExtCroSect,CompSSA,smart,reg,reg_sim);
        term_kappa = 0.5 * get_term_kappa(i1d,j1d,p,tau,tau_grid(ii));
        %term_alpha = (alpha-1)'*log(theta_grid(:,jj));
        Z(ii,jj) = chisq; + term_kappa; %- term_alpha;
    end
end
figure
surf(X,Y,Z),view([-75,30]),colorbar
hold on
contour3(X,Y,Z,50,'white')
hold off
xlabel('\theta_2')
ylabel('AOD')
zlabel('\chi^2_p')
title('Global Convexity')
set(gca,'FontSize',18)
set(get(gca,'xlabel'),'rotation',45)

figure
for j = 1:const.Component_Num
    h = subplot(4,2,j);
    p = get(h, 'pos');
    p(1) = p(1) - 0.05;
    p(3) = p(3) + 0.05;
    set(h, 'pos', p);
    
    plot_1d(r, sample_cd.theta(j,:), x, y, jet, const,[0,1])
    
    title(strcat('Component Num:',num2str(const.Component_Particle(j))))
    set(gca,'FontSize',18)
end 

equ_ref00 = sample_mle_00.atm_path + sample_mle_00.resid;
equ_ref20 = sample_mle_20.atm_path + sample_mle_20.resid;
equ_ref50 = sample_mle_50.atm_path + sample_mle_50.resid;
figure
subplot(4,1,1)
plot_1d(4400,tau,x,y,jet,const)
subplot(4,1,2)
plot_1d(4400,equ_ref00(10,:),x,y,jet,const)
subplot(4,1,3)
plot_1d(4400,equ_ref20(10,:),x,y,jet,const,[0.08,0.1])
subplot(4,1,4)
plot_1d(4400,equ_ref50(10,:),x,y,jet,const,[0.08,0.1])


retri_theta = zeros(9,1);
retri_theta(1) = nanmean(sample_mle_00.theta(1,:));
retri_theta(2) = mean(sample_mle_20.theta(1,:));
retri_theta(3) = mean(sample_mle_50.theta(1,:));
retri_theta(4) = mean(mean(squeeze(sample_mcmc_00.theta(1,:,100:20:end)),2));
retri_theta(5) = mean(mean(squeeze(sample_mcmc_20.theta(1,:,100:20:end)),2));
retri_theta(6) = mean(mean(squeeze(sample_mcmc_50.theta(1,:,100:20:end)),2));
retri_theta(7) = mean(sample_map_00.theta(1,:));
retri_theta(8) = mean(sample_map_20.theta(1,:));
retri_theta(9) = mean(sample_map_50.theta(1,:));
bar(reshape(retri_theta,3,3))
