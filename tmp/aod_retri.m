function sample = aod_retri(Date,Path,Orbit,Block,const)
    
    [~,reg,smart] = load_cache(Date,Path,Orbit,Block,const,'subreg','reg','smart');
    
    XDim_r17600 = const.XDim_r17600;
    YDim_r17600 = const.YDim_r17600;
    XDim_r = const.XDim_r;
    YDim_r = const.YDim_r;
    Band_Green = const.Band_Green;
    Model_ComponentDim = const.Model_ComponentDim;
    Component_Num = const.Component_Num;
    NChannel = const.NChannel;
    header_MIL2ASAE_filename = const.header_MIL2ASAE_filename;
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');

    tau0 = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'RegMeanSpectralOptDepth', ...
        'Index',{[Block  1  1  Band_Green],[1  1  1  1],[1  XDim_r17600  YDim_r17600  1]});
    tau0(tau0==-9999) =  mean(tau0(tau0~=-9999));
    tau0 = double(tau0);
    
    ExtCroSect = hdfread(file_aerosol, '/Component Particle Information/Data Table', 'Fields', ...
        'Spectral extinction cross section', 'FirstRecord',1 ,'NumRecords', Model_ComponentDim);
    ExtCroSect = ExtCroSect{1};
    CompSSA = hdfread(file_aerosol, '/Component Particle Information/Data Table', 'Fields', ...
        'Spectral single scattering albedo', 'FirstRecord',1 ,'NumRecords', Model_ComponentDim);
    CompSSA = CompSSA{1};
    
    [x,y] = find(reg.reg_is_used);
    
    Q = igmrfprec([XDim_r, YDim_r], 1); % precision matrix for Gaussian Markov Random Field on a given grid
    [i2d, j2d] = find(Q); % find nonzero element indices: 1-D from 1 to 4096 for i and j. size(i)=20160
    mask = (i2d ~= j2d) & reg.reg_is_used(i2d)==true & reg.reg_is_used(j2d)==true;
    i = reg.ind_used(i2d(mask));
    j = reg.ind_used(j2d(mask));
       
    iter = 50;
    %tau0_r = kron(tau0, ones(RegScale));
	%current.tau = diag(tau0_r(x,y));
    current.tau = nanmean(tau0(:))*ones(reg.num_reg_used,1);
    
    delta = 0.05;
    sample.tau = zeros(reg.num_reg_used, iter+1);
    sample.tau(:,1) = current.tau;

    theta = zeros(Component_Num,XDim_r,YDim_r);
    for ii = 1:XDim_r
        for jj = 1:YDim_r
            for kk = 1:Component_Num
                theta(kk,ii,jj) = 1/Component_Num;
            end
        end
        
    end
    
    current.theta = zeros(Component_Num,reg.num_reg_used);
    for kk = 1:Component_Num
        current.theta(kk,:) = diag(squeeze(theta(kk,x,y)));
    end
    
    sample.theta = zeros(Component_Num, reg.num_reg_used, iter+1);
    sample.theta(:,:,1) = current.theta;
    
    [current.atm_path,current.surf,current.resid] = update_resid(current, x, y, smart, reg, ExtCroSect, CompSSA);
    sample.atm_path = zeros(NChannel,reg.num_reg_used,iter+1);
    sample.surf = zeros(NChannel,reg.num_reg_used,iter+1);
    sample.resid = zeros(NChannel,reg.num_reg_used,iter+1);

    sample.atm_path(:,:,1) = current.atm_path;
    sample.surf(:,:,1) = current.surf;
    sample.resid(:,:,1) = current.resid;
    
    current.sigmasq = update_sigmasq(current.resid);
    sample.sigmasq = zeros(NChannel, iter+1);
    sample.sigmasq(:,1) = current.sigmasq;
    
    % Initialize kappa
    current.kappa = update_kappa(current.tau,i,j,reg.num_reg_used);
    %current.kappa = 8787;
    sample.kappa = zeros(iter+1,1);
    sample.kappa(1) = current.kappa;
    
    sample.loglik = zeros(iter+1,1);
    sample.loglik(1) = log_lik(current,i,j);

    figure
    for t = 1: iter
        clf
        
        show(sample,reg,2,1,t)
        
        M=getframe;

        [loglik,num] = log_lik(current,i,j);
        fprintf('Round: %d, Log-lik: %.2e, active: %d \n',t,loglik,num);

        [current.tau,current.resid] = update_tau(current, delta, i, j, x, y, smart, reg, ExtCroSect, CompSSA);

        current.kappa = update_kappa(current.tau,i,j,reg.num_reg_used);

        [current.theta,current.resid] = update_theta(current, i, j, x, y, smart, reg, ExtCroSect, CompSSA);

        current.sigmasq = update_sigmasq(current.resid);
        
        [current.atm_path,current.surf,~] = update_resid(current, x, y, smart, reg, ExtCroSect, CompSSA);

        sample.tau(:,t+1) = current.tau;
        sample.kappa(t+1) = current.kappa;
        sample.theta(:,:,t+1) = current.theta;
        sample.sigmasq(:,t+1) = current.sigmasq;
        sample.resid(:,:,t+1) = current.resid;
        sample.loglik(t+1) = log_lik(current,i,j);
        sample.atm_path(:,:,t+1) = current.atm_path;
        sample.surf(:,:,t+1) = current.surf;

    end

    [loglik,num] = log_lik(current,i,j);
    fprintf('Round: %d, Log-lik: %.2e, active: %d \n',t+1,loglik,num);
    

end