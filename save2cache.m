function save2cache(Date,Path,Orbit,Block,const,varargin)
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    Block = num2str(Block);
    
    dir_cache = fullfile('cache');
    dir_result = fullfile('cache/result');
    
    if ~exist(dir_cache,'dir')
        mkdir(dir_cache);
        fprintf('%s is created!\n',dir_cache)
    end
    
    if ~exist(dir_result,'dir')
        mkdir(dir_result)
        fprintf('%s is created!\n',dir_result)
    end

    file_name = strcat(dir_cache,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'.mat');
    
    flag_sample = false;
    for i = 6:nargin
        if strcmp(inputname(i),'sample')
            flag_sample = true;
        end
    end
   
    if flag_sample == true
        
        Method = varargin{nargin-8};
        kf = varargin{nargin-7};
        dy = varargin{nargin-6};
        par = varargin{nargin-5};
        
        file_name_result = strcat(dir_result,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_',Method,const.str_kf{kf+1},const.str_dy{dy+1},const.str_par{par+1},'.mat');
    else
        fprintf('no sample is saved, so no need to specify options!\n')
    end
    
    for i = 6:nargin
        
        var_name = inputname(i);
        
        if strcmp(var_name,'sample')
            sample = varargin{i-5};
            save(file_name_result,'sample')
            fprintf('sample is saved to %s!\n',file_name_result)
        elseif strcmp(var_name,'subreg')
            eval(strcat(var_name,'=varargin{i-5};'));
            save(file_name,var_name)
            fprintf('%s is saved to %s!\n',var_name,file_name)
        elseif strcmp(var_name,'reg') || strcmp(var_name,'smart')
            eval(strcat(var_name,'=varargin{i-5};'));
            save(file_name,var_name,'-append')
            fprintf('%s is saved to %s!\n',var_name,file_name)
        else
            continue
        end
    end

end