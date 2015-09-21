function varargout = load_cache(Date,Path,Orbit,Block,const,varargin)
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    Block = num2str(Block);
    
    dir_cache = fullfile('cache');
    file_name = strcat(dir_cache,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'.mat');
    
    flag_sample = false;
    for i = 1:nargin-5
        if strcmp(varargin(i),'sample')
            flag_sample = true;
        end
    end
   
    if flag_sample == true
        
        Method = varargin{nargin-8};
        kf = varargin{nargin-7};
        dy = varargin{nargin-6};
        par = varargin{nargin-5};
        
        file_name_result = strcat(dir_cache,'/result/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_',Method,const.str_kf{kf+1},const.str_dy{dy+1},const.str_par{par+1},'.mat');
    else
        fprintf('no sample is loaded, so no need to specify options!\n')
    end
    
    for i = 1:nargin-5
        if strcmp(varargin{i},'sample')
            load(file_name_result,'sample')
            fprintf('sample is loaded from %s!\n',file_name_result)
        elseif strcmp(varargin{i},'subreg') || strcmp(varargin{i},'reg') || strcmp(varargin{i},'smart')
            load(file_name,varargin{i})
            fprintf('%s is loaded from %s!\n',varargin{i},file_name)
        else
            continue
        end
        varargout{i} = eval(varargin{i});
    end

end