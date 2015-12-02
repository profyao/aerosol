function [sample,error_flag] = par_aod_retri(Date,Path,Orbit,Block,r,Method,core,const,delta,varargin)
    
    par = true;
    add_limit = false;
    [reg,smart] = load_cache(Date,Path,Orbit,Block,r,'reg','smart');

    if sum(reg.reg_is_used(:))==0
        error_flag = 1;
        sample = [];
        fprintf('no region is used!\n')
        return
    else
        error_flag = 0;
    end
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    
    dir_aerosol = fullfile('products/MIL2ASAE/',Date);
    file_aerosol = strcat(dir_aerosol,'/',const.header_MIL2ASAE_filename,Path,'_O',Orbit,'_F12_0022.hdf');

    tau0 = hdfread(file_aerosol, 'RegParamsAlgDiagnostics', 'Fields', 'RegMeanSpectralOptDepth', ...
        'Index',{[Block  1  1  const.Band_Green],[1  1  1  1],[1  const.XDim_r17600  const.YDim_r17600  1]});
    tau0(tau0==-9999) =  mean(tau0(tau0~=-9999));
    tau0 = double(tau0);
    
    ExtCroSect = hdfread(file_aerosol, '/Component Particle Information/Data Table', 'Fields', ...
        'Spectral extinction cross section', 'FirstRecord',1 ,'NumRecords', const.Model_ComponentDim);
    ExtCroSect = ExtCroSect{1}; % RH and band dependent
    CompSSA = hdfread(file_aerosol, '/Component Particle Information/Data Table', 'Fields', ...
        'Spectral single scattering albedo', 'FirstRecord',1 ,'NumRecords', const.Model_ComponentDim);
    CompSSA = CompSSA{1};

    [x,y] = find(reg.reg_is_used);
    
    XDim_r = const.XDim_r4400 * const.r4400/r;
    YDim_r = const.YDim_r4400 * const.r4400/r;

    Q = igmrfprec([XDim_r, YDim_r], 1); % precision matrix for Gaussian Markov Random Field on a given grid
    [i2d, j2d] = find(Q); % find nonzero element indices: 1-D from 1 to 4096 for i and j. size(i)=20160
    mask = (i2d ~= j2d) & reg.reg_is_used(i2d)==true & reg.reg_is_used(j2d)==true;
    i = reg.ind_used(i2d(mask));
    j = reg.ind_used(j2d(mask));
    
    % simulation mode
    % ratio = varargin{1};
    % data_str = varargin{2};
%     load(data_str,'tau','theta')
%     fprintf('tau and theta are loaded from %s', data_str)
%     surface = zeros(const.Cam_Dim,const.Band_Dim);
%     [reg_sim,~] = get_data(r, XDim_r, YDim_r, const, reg, smart, x, y, tau, theta, ExtCroSect,CompSSA, surface, ratio);
%     reg.mean_equ_ref = reg_sim;
%     reg.min_equ_ref = reg_sim;
    % simulation end
    
    if strcmp(Method,'CD') || strcmp(Method,'CD-noprior')
        iter = 5;
    elseif strcmp(Method,'MCMC')
        iter = 200;
    else
        iter = 20;
    end
    
    %tau0_r = kron(tau0, ones(RegScale));
	%current.tau = diag(tau0_r(x,y));
    current.tau = nanmean(tau0(:))*ones(reg.num_reg_used,1);
    %current.tau = 0.1*ones(reg.num_reg_used,1);
    %current.tau = tau;
    
    %current.theta = theta;
    
    % Initialize kappa
    if strcmp(Method,'CD-random-noprior') || strcmp(Method,'CD-noprior')
        current.kappa = 0;
        current.alpha = ones(const.Component_Num,1);
    else
        current.kappa = update_kappa(current.tau,i,j,reg.num_reg_used,Method);
        current.alpha = ones(const.Component_Num,1);
        %current.alpha = Dirichlet_mle(current.theta',2);
    end
    
    if strcmp(Method,'MCMC') %|| strcmp(Method,'CD') || strcmp(Method,'CD-random')
        theta0 = gen_theta(current.alpha,XDim_r,YDim_r);
        current.theta = theta0(:,reg.reg_is_used);
    else
        current.theta = 1/const.Component_Num * ones(const.Component_Num,reg.num_reg_used);
    end
    
    % for simulation, edit function get_resid()
    [current.atm_path,current.surf,current.resid] = par_update_resid(current.tau,current.theta, x, y, smart, reg, ExtCroSect, CompSSA, par, core, const, r, add_limit);    
    current.sigmasq = update_sigmasq(current.resid,Method);
    
    if strcmp(Method,'MCMC')

        sample.tau = zeros(reg.num_reg_used, iter+1);
        sample.tau(:,1) = current.tau;
        sample.alpha = zeros(const.Component_Num, iter+1);
        sample.alpha(:,1) = current.alpha;
        sample.theta = zeros(const.Component_Num, reg.num_reg_used, iter+1);
        sample.theta(:,:,1) = current.theta;

        sample.atm_path = zeros(const.NChannel,reg.num_reg_used,iter+1);
        sample.surf = zeros(const.NChannel,reg.num_reg_used,iter+1);
        sample.resid = zeros(const.NChannel,reg.num_reg_used,iter+1);

        sample.atm_path(:,:,1) = current.atm_path;
        sample.surf(:,:,1) = current.surf;
        sample.resid(:,:,1) = current.resid;

        sample.sigmasq = zeros(const.NChannel, iter+1);
        sample.sigmasq(:,1) = current.sigmasq;

        sample.kappa = zeros(iter+1,1);
        sample.kappa(1) = current.kappa;

    end
    
    kappa = zeros(iter+1,1);
    kappa(1) = current.kappa;
    alpha = ones(const.Component_Num,iter+1);
    alpha(:,1) = current.alpha;
    sigmasq = zeros(const.NChannel, iter+1);
    sigmasq(:,1) = current.sigmasq;
    
    [current.loglik,num] = log_lik(current,i,j,const.Channel_Used,Method);
    fprintf('Round: 0, Log-lik: %.4e, active: %d \n',current.loglik,num);
    
    tic
    
    %figure
    
    for t = 1: iter
        %clf
        %show(r,current,reg,2,1,jet(256),const)   
        %getframe;

        [current.tau,current.resid] = par_update_tau(current.tau,current.theta,current.resid,current.kappa,current.sigmasq,...
            delta,i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, r, par, core, add_limit, const);
        
        if ~(strcmp(Method,'CD-random-noprior') || strcmp(Method,'CD-noprior'))
            current.kappa = update_kappa(current.tau,i,j,reg.num_reg_used,Method);
        end

        [current.theta, current.resid] = par_update_theta(current.theta,current.tau,current.resid,current.sigmasq,current.alpha,...
            i, j, x, y, smart, reg, ExtCroSect, CompSSA, Method, r, par, core, add_limit, const);

        if strcmp(Method,'MCMC') %|| strcmp(Method,'CD') || strcmp(Method,'CD-random')
            current.alpha = update_alpha(current.alpha,current.theta',const.Component_Num,reg.num_reg_used, Method);
        end

        current.sigmasq = update_sigmasq(current.resid,Method);

        [current.atm_path,current.surf,~] = par_update_resid(current.tau,current.theta, x, y, smart, reg, ExtCroSect, CompSSA, par,core,const,r, add_limit);
        
        if strcmp(Method,'MCMC')
            sample.tau(:,t+1) = current.tau;
            sample.kappa(t+1) = current.kappa;
            sample.theta(:,:,t+1) = current.theta;
            sample.alpha(:,t+1) = current.alpha;
            sample.sigmasq(:,t+1) = current.sigmasq;
            sample.resid(:,:,t+1) = current.resid;
            sample.atm_path(:,:,t+1) = current.atm_path;
            sample.surf(:,:,t+1) = current.surf;
        else
            alpha(:,t+1) = current.alpha;
            kappa(t+1) = current.kappa;
            sigmasq(:,t+1) = current.sigmasq;
        end
        
        %[sample.loglik(t+1),num] = log_lik(current,i,j,const.Channel_Used,Method);
        %fprintf('Round: %d, Log-lik: %.4e, active: %d \n',t,sample.loglik(t+1),num);
        fprintf('.');
        
        if mod(t,50)==0
            fprintf('\n')
        end

    end
    
    toc
    
    [current.loglik,num] = log_lik(current,i,j,const.Channel_Used,Method);
    fprintf('\nRound: %d, Log-lik: %.4e, active: %d \n',t,current.loglik,num);

    if ~strcmp(Method,'MCMC')
        sample = current;
        sample.kappa = kappa;
        sample.sigmasq = sigmasq;
        sample.alpha = alpha;
    else
        sample.loglik = current.loglik;
    end

end