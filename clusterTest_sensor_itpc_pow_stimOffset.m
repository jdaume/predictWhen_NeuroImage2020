%% Cluster statistics on sensor level for ITPC and power only for the enlarged disappearance analysis window
% JD 2020

clear variables

% Fieldtrip
addpath('.../fieldtrip-20170607');  % change to personal location of fieldtrip

ft_defaults;

%% files & folders
subjects = 1:23;
sessions = 1:2;

files_base = fullfile(pwd);
files_meg = fullfile(files_base, 'data_meg_itpc_stimOffset/');
files_bsl_pow = fullfile(files_base, 'data_meg_wlconv_wholeTrial/');
files_evoked_pow = fullfile(files_base, 'data_meg_evoked_pow_stimOffset/');

whichTrials = '_stratSubjCorr';
whichMethod = '_wl';

blc_pow = 1; % baseline correction for power?
twl_bsl = 0.3; %baseline length in s, starts with first time sample in baseline
bsl_offset = 4; % in number of time windows
twl = 0.1; % length of time shifts
bsl_ind = bsl_offset:twl_bsl/twl+bsl_offset;

whichDirection = 0; % 0 = both, 1 = l2r, -1 = r2l
flipSenslr = 1; % flip sensors left to right 1 = yes, 0 = no.

whichTasks = 1:3;

%% Load data and perform blc on each subject and condition

if whichDirection == 1
    subjects = subjects(mod(subjects,4)==1 | mod(subjects,4)==2);
elseif whichDirection == -1
    subjects = subjects(mod(subjects,4)==3 | mod(subjects,4)==0);
end

n_subjects = length(subjects);

for i_task = whichTasks
    
    noSubj = 0;
    
    for i_subj = subjects % loop subjects
        
        noSubj = noSubj + 1;
        
        for i_session = sessions
            
            if i_subj == 18 && i_session == 2
                continue
            end
            
            %free memory
            clear itpc
            
            % get MEG data filename for this subject
            fprintf('Loading VP_%02d_%d_Cond_%d%s%s.mat \n', i_subj,i_session,i_task,whichTrials,whichMethod);
            eval(sprintf('this_meg_data = fullfile(files_meg, [''VP_%02d_%d_Cond_%d%s%s.mat'']);',i_subj,i_session,i_task,whichTrials,whichMethod));
            if size(this_meg_data,1) > 1, error('Error! More than one MEG file!'), end
            if i_subj == subjects(1) && i_session == 1 || (i_subj == subjects(1) && length(sessions) == 1)% loads faster
                itpc = load(this_meg_data);
                [nchan,nf,nt]= size(itpc.itpc);
                eval(sprintf('cond%d = itpc;',i_task));
                eval(sprintf('cond%d.dimord    = ''subj_chan_freq_time'';',i_task));
                eval(sprintf('cond%d.itpc = [];',i_task));
                eval(sprintf('cond%d.pow = [];',i_task));
                eval(sprintf('cond%d.subj = {itpc.subj};',i_task));
                if flipSenslr
                    eval(sprintf('[~, indExcl] = intersect(cond%d.label, {''MRC53'';''MRF21'';''MLO42'';''MRP11'';''MRP54''});',i_task)); % these dont have a counterpart on the other side
                    eval(sprintf('cond%d.label(indExcl) = [];',i_task)); % exclude those without counterpart
                    eval(sprintf('nchan = size(cond%d.label,1);',i_task));
                    eval(sprintf('[~, indL] = intersect(cond%d.label, ft_channelselection(''ML*'', cond%d.label));',i_task,i_task));
                    eval(sprintf('[~, indR] = intersect(cond%d.label, ft_channelselection(''MR*'', cond%d.label));',i_task,i_task));
                    eval(sprintf('[~, indZ] = intersect(cond%d.label, ft_channelselection(''MZ*'', cond%d.label));',i_task,i_task));
                end
                eval(sprintf('cond%d.itpc = zeros(n_subjects,nchan,nf,nt);',i_task));
                eval(sprintf('cond%d.pow = zeros(n_subjects,nchan,nf,nt);',i_task));
            else
                itpc = load(this_meg_data,'itpc','subj','pow');
                eval(sprintf('cond%d.subj{noSubj} = itpc.subj;',i_task));
            end

            
            clear pow bsl pow_blc
            if blc_pow
                eval(sprintf('this_pow_data = fullfile(files_bsl_pow, [''VP_%02d_%d_Cond_%d%s_pow.mat'']);',i_subj,i_session,i_task,whichTrials)); % takes the "real" baseline from the wholeTrial analysis
                evPow = load(this_pow_data,'pow');
                pow_blc = mean(evPow.pow(:,:,bsl_ind),3);
                itpc.pow  = (itpc.pow ./ repmat(pow_blc, [1 1 nt])-1) * 100;
            end

            % Flip sensors left to right (if indicated above)
            if flipSenslr
                itpc.itpc(indExcl,:,:) = [];
                itpc.pow(indExcl,:,:) = [];
                if (mod(i_subj,4)==3 || mod(i_subj,4)==0)
                    itpc.itpc = itpc.itpc([indR;indL;indZ],:,:);
                    itpc.pow = itpc.pow([indR;indL;indZ],:,:);
                end
            end
            
            % write data into cond variable
            if i_subj == 18 % has only 1 session
                eval(sprintf('cond%d.itpc(noSubj,:,:,:) = squeeze(cond%d.itpc(noSubj,:,:,:)) + itpc.itpc;',i_task,i_task));
                eval(sprintf('cond%d.pow(noSubj,:,:,:) = squeeze(cond%d.pow(noSubj,:,:,:)) + itpc.pow;',i_task,i_task));
            else
                eval(sprintf('cond%d.itpc(noSubj,:,:,:) = squeeze(cond%d.itpc(noSubj,:,:,:)) + itpc.itpc/2;',i_task,i_task));
                eval(sprintf('cond%d.pow(noSubj,:,:,:) = squeeze(cond%d.pow(noSubj,:,:,:)) + itpc.pow/2;',i_task,i_task));
            end
            
        end
    end
    % Number of subjects correct?
    test = n_subjects - noSubj;
    if test ~= 0, error('No. of subjects not correct!'), end
end

disp('All subject data loaded!');

return


%%  T-Test between two conditions (cluster permutation)

% Compute cluster statistics
clear stats task1 task2

% cond1 = Visual temporal prediction
% cond2 = Tactile temporal prediction
% cond3 = Luminance matching

task1 = cond1;  
task2 = cond3;


cfg                  = [];
cfg.channel          = 'all';
cfg.latency          = [-1.9 1.9];
cfg.frequency        = [0.5 3]; 
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'depsamplesT';
cfg.computeprob      = 'yes';
cfg.correctm         = 'cluster'; 
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0; % two-tailed
cfg.tail             = 0; % two-tailed
cfg.alpha            = 0.025/2; % two tests: 1-3 & 2-3
cfg.numrandomization = 1000;
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'yes';
cfg.avgoverchan      = 'no';
cfg.parameter        = 'itpc'; % 'pow'


if strcmp(cfg.avgoverchan,'no')
    cfg.minnbchan    = 2;
end

% specifies with which sensors other sensors can form clusters
cfg_neighb.method           = 'template';
cfg_neighb.template         = 'CTF275_neighb.mat';
cfg_neighb.feedback         = 'no';
cfg.neighbours              = ft_prepare_neighbours(cfg_neighb, cond1);

design = zeros(2,2*n_subjects);
design(1,:) = repmat(1:n_subjects,1,2);
design(2,:) = mod(floor([0:(2*n_subjects-1)]/(n_subjects/1)),2)+1;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stats] = ft_freqstatistics(cfg, task1,task2);

    
% Clusterplot
cfg = [];
cfg.alpha  = 0.025/2;
cfg.marker  = 'off';
cfg.highlightseries = {'on', 'on', 'on', 'on', 'on'};
cfg.highlightsymbolseries = ['.','.','o','+','s'];
cfg.highlightsizeseries  = [10 10 15 15 15];
cfg.highlightcolorpos = [0 0 0];
cfg.hotkeys      = 'yes';
cfg.style        = 'straight';
cfg.gridscale    = 200;
cfg.shading      = 'interp';
cfg.interactive  = 'no';
cfg.parameter    = 'stat';
cfg.colorbar     = 'no';
cfg.zlim         = [-4 4];
cfg.layout       = 'CTF275.lay';
cfg.subplotsize  = [3 4];

ft_clusterplot(cfg, stats);

