function save2cache(Date,Path,Orbit,Block,const,varargin)

    Component_Particle = const.Component_Particle;
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    Block = num2str(Block);
    
    dir_cache = fullfile('cache');
    
    if ~exist(dir_cache,'dir')
        mkdir(dir_cache);
        fprintf('%s is created!\n',dir_cache)
    else
        fprintf('directory %s exists, continue to save files!\n',dir_cache)
    end

    file_name = strcat(dir_cache,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'.mat');
    
    if nargin~=8
        Method = varargin{nargin-5};
        file_name_result = strcat(dir_cache,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_Comp',mat2str(Component_Particle),'_',Method,'.mat');
    else
        fprintf('no need for method!\n')
    end
    
    for i = 6:nargin
        
        var_name = inputname(i);
        
        if strcmp(var_name,'sample')
            sample = varargin{i-5};
            save(file_name_result,'sample')
            fprintf('sample is saved to %s!\n',file_name_result)
        elseif (strcmp(var_name,'subreg') || strcmp(var_name,'reg') || strcmp(var_name,'smart')) && i==6
            eval(strcat(var_name,'=varargin{i-5};'));
            save(file_name,var_name)
            fprintf('%s is saved to %s!\n',var_name,file_name)
        elseif (strcmp(var_name,'subreg') || strcmp(var_name,'reg') || strcmp(var_name,'smart')) && i>6
            eval(strcat(var_name,'=varargin{i-5};'));
            save(file_name,var_name,'-append')
            fprintf('%s is saved to %s!\n',var_name,file_name)
        else
            continue
        end
    end

end