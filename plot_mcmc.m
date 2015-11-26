% figure
% for i=1:8
%     subplot(4,2,i)
%     plot(acf(squeeze(sample.theta(i,1847,:)),1000))
%     title(strcat('Component Num:',num2str(const.Component_Particle(i))))
% end
% 
% figure
% for i=1:8
%     subplot(4,2,i)
%     plot(cumsum(squeeze(sample.theta(i,1847,:)))./[1:1001]')
%     title(strcat('Component Num:',num2str(const.Component_Particle(i))))
% end

%[sample,reg] = load_cache('2011.07.29',15,61764,59,4400,'sample','reg','CD-random');
[sample,reg] = load_cache('2011.06.02',16,60934,59,4400,'sample','reg','MCMC-G');
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
