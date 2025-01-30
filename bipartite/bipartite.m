function bipartite(varargin)
tictime = clock(); 
script_path = fileparts(mfilename('fullpath'));
addpath(genpath([script_path, '/../../utils']));
addpath([script_path, '/../']);


%% input
ip = inputParser;
ip.KeepUnmatched = true;
ip.PartialMatching = false;

% ip.addParameter('cir', 14);
ip.addParameter('boundary', 'periodic');
ip.addParameter('shift', 0);
ip.addParameter('T', 150);
ip.addParameter('verbose', false);
% ip.addParameter('probs', [0.25, 0.25, 0.25, 0.25]);
ip.addParameter('num', 10);
ip.addParameter('plotting', 1);
ip.addParameter('update', 0);
ip.addParameter('init', 'flux_free');

ip.parse(varargin{:});
pars = ip.Results;

% Convert pars to a cell array of parameter-value pairs
fields = fieldnames(pars);      
values = struct2cell(pars);     
varargin = reshape([fields, values]', 1, []); 

%% Prepare input
cirs = [18, 24, 30, 36, 42, 48, 54, 60];  
cirs = [6, 10, 12, 14, 16, 18, 24, 30, 36, 42];  
cirs = [6, 10, 12, 14, 16, 18, 24, 30];  
% cirs = [10];
probs = [0.25, 0.25, 0.25, 0.25];


for i = 1: length(cirs)
    L = cirs(i);
    
    % if L >= 36
    %     pars.T = 300;
    % end
    
    %% save
    name = sprintf('cir_%d_T_%d_probs_%.2f_%.2f_%.2f_%.2f_%s',...
        L, pars.T, probs(1), probs(2), probs(3), ...
        probs(4), pars.boundary);
    
    fprintf('%s\n', name);                
    data_folder = sprintf([script_path, '/data_%s/'], mfilename);
    triexpf(~exist(data_folder, 'dir'), {@mkdir, data_folder}, {});
    if ~ pars.update && exist([data_folder, name, '.mat'], 'file')
        loaded = load([data_folder, name, '.mat']);
        fields_loaded = fieldnames(loaded);
        for idx = 1:length(fields_loaded)
            if strcmp(fields_loaded{idx}, 'script_path') || ...
                    strcmp(fields_loaded{idx}, 'fields_loaded') || ...
                    strcmp(fields_loaded{idx}, 'cirs')
                continue
            end
            eval([fields_loaded{idx} ' = loaded.(fields_loaded{idx});']);
        end
        disp('file exists, loaded');
    else
        disp('update or file not exist');
        
        % Compute entropy
        ls = 0: max(1, floor(L / pars.num)): L;
        ls = unique(sort([ls, floor(L / 2)]));
        entropies = zeros(size(ls));

        % Purify load
        res = purify(varargin{:}, 'T', pars.T, 'cir', L, 'probs', probs, ...
            'plotting', 0, 'update', 0);
    
        tableau = res.simul.tableau.clone();
        % qb_last = [];
        for idx = fliplr(1: length(ls))
            fprintf('%d/ %d = %.2f\n', idx, length(ls), idx/ length(ls))
            sys_size = ls(idx);
            qubits_env = res.simul.check_generator.len2qubits(L - sys_size, 'start', 1 + sys_size);
            % if qb_last
            %     assert(all(ismember(qb_last, qubits_env)))
            % end
            % qb_last = qubits_env;
            tableau.partial_trace(qubits_env);

            entropies(idx) = tableau.get_entropy();
        end 
    
        %% save
        pars_tmp = pars;
        clear pars fields_loaded
        save([data_folder, name, '.mat']);
        fprintf('data saved:\n %s\n', [data_folder, name]);
        fprintf('Time elapsed: %f s\n',  etime(clock(), tictime));
        pars = pars_tmp;
    end
end        


end
