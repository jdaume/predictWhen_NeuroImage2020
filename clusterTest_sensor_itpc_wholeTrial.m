%% Cluster statistics on sensor level for ITPC and the wholeTrial windows
% JD 2020

clear variables

% Fieldtrip
addpath('.../fieldtrip-20170607');  % change to personal location of fieldtrip

ft_defaults;

%% files & folders

subjects = 1:23;
sessions = 1:2;


files_base = fullfile(pwd);
files_meg = fullfile(files_base, 'data_meg_wlconv_wholeTrial/');

whichTrials = '_stratSubjCorr';


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
            fprintf('Loading VP_%02d_%d_Cond_%d%s_itpc.mat \n', i_subj,i_session,i_task,whichTrials);
            eval(sprintf('this_meg_data = fullfile(files_meg, [''VP_%02d_%d_Cond_%d%s_itpc.mat'']);',i_subj,i_session,i_task,whichTrials));
            if size(this_meg_data,1) > 1, error('Error! More than one MEG file!'), end
            itpc = load(this_meg_data);
            
            % round time so that we have a clear 0
            itpc.time = round(itpc.time,2);
            
            % Flip sensors left to right (if indicated above)
            if flipSenslr 
                [~, indExcl] = intersect(itpc.label, {'MRC53';'MRF21';'MLO42';'MRP11';'MRP54'}); % these dont have a counterpart on the other side
                itpc.label(indExcl) = []; %exlude them
                [~, indL] = intersect(itpc.label, ft_channelselection('ML*', itpc.label));
                [~, indR] = intersect(itpc.label, ft_channelselection('MR*', itpc.label));
                [~, indZ] = intersect(itpc.label, ft_channelselection('MZ*', itpc.label));
                itpc.itpc(indExcl,:,:) = [];
                if mod(i_subj,4)==3 || mod(i_subj,4)==0
                    itpc.itpc = itpc.itpc([indR;indL;indZ],:,:); %flip
                end
            end

            % make cond variable
            if i_subj == subjects(1) && i_session == 1 || (i_subj == subjects(1) && length(sessions) == 1)
                eval(sprintf('cond%d = itpc;',i_task));
                eval(sprintf('cond%d.dimord = ''subj_chan_freq_time'';',i_task));
                eval(sprintf('cond%d.itpc = [];',i_task));
                eval(sprintf('cond%d.subj  = {itpc.subj};',i_task));
                eval(sprintf('cond%d.itpc = zeros([n_subjects,size(itpc.itpc)]);',i_task));
            else
                eval(sprintf('cond%d.subj{noSubj}      = itpc.subj;',i_task));
            end
            
            
            % write data into cond variable
            if i_subj == 18 % has only 1 session
                eval(sprintf('cond%d.itpc(noSubj,:,:,:) = squeeze(cond%d.itpc(noSubj,:,:,:)) + itpc.itpc;',i_task,i_task));
            else
                eval(sprintf('cond%d.itpc(noSubj,:,:,:) = squeeze(cond%d.itpc(noSubj,:,:,:)) + itpc.itpc/2;',i_task,i_task));
            end
            
        end
    end
    % Number of subjects correct?
    test = n_subjects - noSubj;
    if test ~= 0, error('No. of subjects not correct!'), end
end

% Condition average
cond0 = cond1;
cond0.itpc = (cond1.itpc + cond2.itpc + cond3.itpc)/3;

disp('All subject data loaded!');

return


%% Condition (average) against baseline Clusterstat (2D)
data2plot = [];

for i = 1:3
    
    clear stats task1 task2
    
    
    task1 = cond0;  %
    
    % make task2 from baseline
    task2 = task1;
    task2.itpc = [];
    task2.itpc = repmat(mean(task1.itpc(:,:,:,bsl_ind),4),[1 1 1 length(task1.time)]);
    
    
    if i == 1
        latency = [itpc.time_stimMov(1) itpc.time_stimMov(end)];
    elseif i == 2
        latency = [itpc.time_stimOff(1) itpc.time_stimOff(end)];
    elseif i == 3
        latency = [itpc.time_stimReapp(1) itpc.time_stimReapp(end)];
    end
    
    
    cfg                  = [];
    % cfg.parameter        = 'powspctrm';
    cfg.channel          = 'all';
    cfg.latency          = latency;
    cfg.frequency        = 'all';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';
    cfg.computeprob      = 'yes';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.clustertail      = 0; 
    cfg.tail             = 0; 
    cfg.alpha            = 0.025/3; % three windows
    cfg.numrandomization = 1000;
    cfg.avgovertime      = 'no';
    cfg.avgoverfreq      = 'no';
    cfg.avgoverchan      = 'yes';
    cfg.parameter        = 'itpc';
    
    % specifies with which sensors other sensors can form clusters
    cfg_neighb.method           = 'template';
    cfg_neighb.template         = 'CTF275_neighb.mat';
    cfg_neighb.feedback         = 'no';
    cfg.neighbours          = ft_prepare_neighbours(cfg_neighb, task1);
    
    design = zeros(2,2*n_subjects);
    design(1,:) = repmat(1:n_subjects,1,2);
    design(2,:) = mod(floor([0:(2*n_subjects-1)]/(n_subjects/1)),2)+1;
    
    
    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
    
    [stats] = ft_freqstatistics(cfg, task1 , task2);
    data2plot = [data2plot,squeeze(stats.mask .* stats.stat)];
    
end

%% Plot stats mask
% Plot
figure('color',[1 1 1],'position',[300 300 1500 750]);
titleVec = {'Baseline','Movement', 'Occluder', 'Reappearance'};

freq2plot_HF = [40 100]; % in Hz
freq2plot_LF = [0.5 40]; % in Hz
resolution = 100;
time = {itpc.time_bsl, itpc.time_stimMov, itpc.time_stimOff, itpc.time_stimReapp};
time = [time,time];

for isp = 1:8 
    h(isp) = subplot(2,4,isp);
    time2plot = [time{isp}(1) time{isp}(end)];
    xlabel('Time (s)','fontsize',15);
    if mod(isp,4)==1, ylabel('Frequency (Hz)','fontsize',15,'fontname','Arial'); end
    if isp <= 4, freq2plot = freq2plot_HF; title(titleVec{isp}); else, freq2plot = freq2plot_LF; end
    set(h(isp),'xlim',time2plot,'ylim',freq2plot);
    if isp == 1 || isp == 5, continue, end
    contourf(itpc.time(9:end), itpc.freq , data2plot2,resolution,'linecolor','none'); %
    set(h(isp),'xlim',time2plot,'ylim',freq2plot);
end

