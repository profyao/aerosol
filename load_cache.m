function varargout = load_cache(Date,Path,Orbit,Block,r,varargin)
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    Block = num2str(Block);
    
    dir_data = fullfile('cache/data');
    file_name = strcat(dir_data,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_R',num2str(r),'.mat');
    
    flag_sample = false;
    flag_sim = false;
    
    var_num = 5;
    
    for i = 1:nargin-var_num
        if strcmp(varargin(i),'sample')
            flag_sample = true;
        elseif strcmp(varargin{i},'reg_sim')
            flag_sim = true;
        end
    end
   
    if flag_sample == true
        Method = varargin{nargin-var_num};
        file_name_result = strcat('cache/result/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_R',num2str(r),'_',Method,'.mat');
    elseif flag_sim == true
        Method = varargin{nargin-var_num};
        file_name_sim = strcat(dir_data,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_R',num2str(r),'_',Method,'.mat');
    end
    
    for i = 1:nargin-var_num
        if strcmp(varargin{i},'sample')
            load(file_name_result,'sample')
            fprintf('sample is loaded from %s!\n',file_name_result)
        elseif strcmp(varargin{i},'reg_sim')
            load(file_name_sim,'reg_sim')
            fprintf('reg_sim is loaded from %s!\n',file_name_sim)
        elseif strcmp(varargin{i},'reg') || strcmp(varargin{i},'smart')
            load(file_name,varargin{i})
            fprintf('%s is loaded from %s!\n',varargin{i},file_name)
        else
            continue
        end
        varargout{i} = eval(varargin{i});
    end

end