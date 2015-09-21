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

    iters = 0:log2(num_cores);
    t = NaN*ones(1,length(iters));
    cores = 2.^iters;
    
    for i = iters+1
        core = cores(i);
        parpool(core)
        tic
        par_aod_retri_batch(Method,Location,const)
        t(i) = toc;
        delete(cl.Jobs)
        delete(gcp)
    end

end