function varargout = load_cache(Date,Path,Orbit,Block,const,varargin)

    Component_Particle = const.Component_Particle;
    
    Orbit = num2str(Orbit,'%06d');
    Path = num2str(Path,'%03d');
    Block = num2str(Block);
    
    dir_cache = fullfile('cache');
    file_name = strcat(dir_cache,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'.mat');
    
    Method = varargin{nargin-5};
    file_name_result = strcat(dir_cache,'/',Date,'_P',Path,'_O',Orbit,'_B',Block,'_Comp',mat2str(Component_Particle),'_',Method,'.mat');
        
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