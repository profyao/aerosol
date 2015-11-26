function save2cache(Date,Path,Orbit,Block,r,varargin)
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    Block = num2str(Block);
    
    dir_data = fullfile('cache/data');
    dir_result = fullfile('cache/result');
    
    if ~exist(dir_data,'dir')
        mkdir(dir_data);
        fprintf('%s is created!\n',dir_data)
    end
    
    if ~exist(dir_result,'dir')
        mkdir(dir_result)
        fprintf('%s is created!\n',dir_result)
    end

    file_name = strcat(dir_data,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_R',num2str(r),'.mat');
    
    flag_sample = false;
    flag_sim = false;
    
    var_num = 5;
    
    for i = 1:nargin - var_num
        if strcmp(inputname(i+var_num),'sample') 
            flag_sample = true;
        elseif strcmp(inputname(i+var_num),'reg_sim')
            flag_sim = true;
        end
    end
   
    if flag_sample == true
        Method = varargin{nargin-var_num};
        file_name_result = strcat(dir_result,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_R',num2str(r),'_',Method,'.mat');
    elseif flag_sim == true
        Method = varargin{nargin-var_num};
        file_name_sim = strcat(dir_data,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_R',num2str(r),'_',Method,'.mat');
    end

    for i = 1:nargin-var_num
        
        var_name = inputname(i+var_num);
        
        if strcmp(var_name,'sample') 
            sample = varargin{i};
            save(file_name_result,'sample')
            fprintf('sample is saved to %s!\n',file_name_result)
        elseif strcmp(var_name,'reg_sim')
            reg_sim = varargin{i};
            save(file_name_sim,'reg_sim')
            fprintf('reg_sim is saved to %s!\n',file_name_sim)
        elseif strcmp(var_name,'reg') 
            eval(strcat(var_name,'=varargin{i};'));
            save(file_name,var_name)
            fprintf('%s is saved to %s!\n',var_name,file_name)
        elseif strcmp(var_name,'smart')
            eval(strcat(var_name,'=varargin{i};'));
            save(file_name,var_name,'-append')
            fprintf('%s is saved to %s!\n',var_name,file_name)
        else
            continue
        end
    end

end