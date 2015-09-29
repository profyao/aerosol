function [cores,t] = par_performance(server,Method,Location,const)
    
    if strcmp(server,'SCF')
        %num_cores = str2double(getenv('NSLOTS'));
        num_cores = 32;
        cl = parcluster('local');
        cl.NumWorkers = 32;
        saveProfile(cl)
    elseif strcmp(server,'local')
        num_cores = 4;
        cl = parcluster('local');
    else
        error('server is not specified!')
    end

    iters = log2(num_cores):-1:0;
    t = NaN*ones(1,length(iters));
    cores = 2.^iters;
    
    if ~exist('timing','dir')
        mkdir('timing')
    end
    
    for i = 1:length(iters)
        core = cores(i);
        
        if core == 1
            par = 0;
        else
            parpool(core);
            par = 1;
        end
        
        %profile on
        f = @() par_aod_retri_batch(Method,Location,0,0,par,core,const,0);
        t(i) = timeit(f);
        %p = profile('info');
        %profile off
               
        if core~=1
            delete(cl.Jobs)
            delete(gcp)
        end
        
    end
    
    save(strcat('timing/runtime_',Method,'.mat'),'cores','t');

end