%% Complex morelet wavelet convolution script analyzing data from prolonged disappearance (stimulus offset) analysis window
% this script computes itpc, total power, and the ERF for a time window around stimulus disappearance (time 0). 

% JD 2020


clear variables

%% Fieldtrip
addpath('.../fieldtrip-20170607'); % change to personal location of fieldtrip
ft_defaults

%%

headposThresh = 5; % kick out all trials above this threshold in any coil in mm (trials are saved in plot_and_cutoff_headpos); 0 = no cutoff!
correct = 1; % 1 = correct trials, 0 = all trials
stratified = 1; % 1 = stratified number of trials, 0 = all (correct) trials
subjCorr = 1; % only works when correct = 1; 1 = subjectively correct trials, 0 = objectively correct trials
inducedDisapp = 0; % subtracts the erf from disappearance

eyetracker = 0; % 1 = uses eyetracker data, 0 = uses MEG data
twl = 0.1; % time window length for output in s
toilim = [-2 2]; % time window of interest in s; output data will be twl s shorter on each side

if correct && stratified && ~subjCorr
    whichTrials = '_stratCorr';   
elseif correct && stratified && subjCorr
    whichTrials = '_stratSubjCorr';       
elseif correct && ~stratified
    whichTrials = '_corr';
else
    whichTrials = '_all';
end

if inducedDisapp
    induced = '/inducedDisapp';
else
    induced = ''; 
end

if eyetracker
    whichData = 'eyetracker/';
else
    whichData = '';
end


%% files & folders
if eyetracker
    files_data = fullfile(pwd, 'data_eyetracking');
else
    files_meg_preproc = fullfile(pwd, 'data_meg_preprocessed');
end
eval(sprintf('files_meg_itpc = fullfile(pwd, ''data_meg_itpc_stimOffset/%s%s'');',whichData,induced));
files_headpos = fullfile(pwd, 'data_headpos');
eval(sprintf('files_logfiles = fullfile(pwd, ''logfiles/itpc_stimOffset/%s%s'');',whichData,induced));
files_strat = fullfile(pwd, 'data_stratTrials');
if ~exist(files_meg_itpc,'dir'), mkdir(files_meg_itpc); end
if ~exist(files_logfiles,'dir'), mkdir(files_logfiles); end


%% Trialinfo and stratified trials
trialinfo_predictWhen; % returns variable 'col' containing info about which column represents which info in trialinfo matrix

if correct && stratified && ~subjCorr
    eval(sprintf('load(fullfile(files_strat,''stratTrials_corr_headpos%dmm.mat''));',headposThresh));
elseif correct && stratified && subjCorr
    eval(sprintf('load(fullfile(files_strat,''stratTrials_subjCorr_headpos%dmm.mat''));',headposThresh));
end

%% loop subjects
for i_subj = 1:23
    for i_ses = 1:2
        
        if i_subj == 18 && i_ses == 2
            continue;
        end
                
        if ~exist(sprintf([files_logfiles '/VP%02d_%d%s_wl_processing.txt'],i_subj,i_ses,whichTrials),'file')
            
            system(['touch ' files_logfiles sprintf('/VP%02d_%d%s_wl_processing.txt',i_subj,i_ses,whichTrials)]);
            
            %free memory
            clear data
            
            
            % get MEG data filename for this subject
            if eyetracker
                fprintf('Loading eyetracking-data for VP_%02d_%d\n',i_subj,i_ses);
                eval(sprintf('this_meg_data = fullfile(files_data, [''VP_%02d_%d.mat'']);',i_subj,i_ses));
            else
                fprintf('Loading MEG-data for VP_%02d_%d\n',i_subj,i_ses);
                eval(sprintf('this_meg_data = fullfile(files_meg_preproc, [''VP_%02d_%d_compRejected.mat'']);',i_subj,i_ses));
            end
            if size(this_meg_data,1) > 1, error('Error! More than one MEG file!'), end
            data = load(this_meg_data);
                     
            
            %% Kick out trials above headposition threshold
            if headposThresh
                eval(sprintf('load(fullfile(files_headpos,  [''VP_%02d_%d_trialInd_headposAbove%dmm.mat'']));',i_subj,i_ses,headposThresh));
            else % no headpos cutoff
                trialInd_headpos = [];
            end
            
            cfg = [];
            cfg.trials = ~ismember(data.trialinfo(:,col.trial),trialInd_headpos);
            data       = ft_preprocessing(cfg,data);
            
            %% redefine trials to equal length (here just for timing information)
            cfg = [];
            cfg.toilim = toilim;
            
            data = ft_redefinetrial(cfg, data);
            
            %% Wavelet convolution parameters
            
            % wavelet parameters
            min_freq     = 0.5;
            max_freq     = 100;
            n_freq       = 40;
            f            = logspace(log10(min_freq),log10(max_freq),n_freq);
            wl_time      = -2.5:1/data.fsample:2.5; 
            wl_time_half = (length(wl_time)-1)/2;
            
            % FFT parameters (use next-power-of-2)
            n_samples_wl          = length(wl_time);
            n_samples_data        = size(data.time{1},2);
            n_samples_convolution = n_samples_wl+n_samples_data-1;
            n_samples_conv_pow2   = pow2(nextpow2(n_samples_convolution));
            wavelet_cycles        = logspace(log10(2),log10(10),n_freq);
            
            % Computation and FFT of wavelets
            fprintf('Computation and fft of wavelets for all %d specified frequencies...\n',n_freq)
            % compute wavelets for each frequency
            wavelets = zeros(n_freq, n_samples_wl);
            for fi = 1:n_freq
                wavelets(fi,:) = (pi*f(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*f(fi).*wl_time) .* exp(-wl_time.^2./(2*( wavelet_cycles(fi) /(2*pi*f(fi)))^2))/f(fi);
            end
            
            % Fourier transform of the wavelets
            wavelets_fft = fft(wavelets, n_samples_conv_pow2, 2);
            
            
            %% Loop conditions/tasks
            for i_cond = 1:3
                
                if correct && ~stratified
                    fprintf('Analyzing correct trials from task %d\n',i_cond)
                elseif correct && stratified && ~subjCorr
                    fprintf('Analyzing objectively correct and stratified trials from task %d\n',i_cond)
                elseif correct && stratified && subjCorr 
                    fprintf('Analyzing subjectively correct and stratified trials from task %d\n',i_cond)
                else
                    fprintf('Analyzing all trials from task %d\n',i_cond)
                end

                %% Select trials and set parameters
                
                clear data_task cond_all
                
                % Select trials and channels
                cfg = [];
                if correct && ~stratified
                    cfg.trials     = data.trialinfo(:,col.condition)==i_cond & data.trialinfo(:,col.acc)==1;
                elseif correct && stratified && ~subjCorr
                    cond_all       = find(data.trialinfo(:,col.condition)==i_cond);
                    cfg.trials     = cond_all(stratTrials{i_subj,i_ses,i_cond});
                elseif correct && stratified && subjCorr
                    cond_all       = find(data.trialinfo(:,col.condition)==i_cond);
                    cfg.trials     = cond_all(stratTrials_subjCorr{i_subj,i_ses,i_cond});
                else
                    cfg.trials     = data.trialinfo(:,col.condition)==i_cond;
                end
                data_task      = ft_preprocessing(cfg, data);
                n_trials       = size(data_task.trial,2);
                n_channels     = size(data_task.label,1);
                 
                %% Write data in 3D matric and subtract ERF if desired
                
                clear datMat
                
                % convert cell array to mat, channels x samples x trials
                disp('Writing data into 3D matrix...')
                datMat = cat(3,data_task.trial{:});  
                
                % Subtract ERF
                if inducedDisapp
                    clear erf
                    erf = nanmean(datMat,3);
                    datMat = datMat - repmat(erf, [1 1 n_trials]);
                end
                
                %% create FT output structure for itpc
                
                %free memory
                clear itpc
                
                itpc         = [];
                eval(sprintf('itpc.subj    = ''VP_%02d'';',i_subj));
                itpc.session = i_ses;
                itpc.method  = 'wlconv';
                itpc.label   = data_task.label;
                itpc.fsample = data.fsample;
                itpc.dimord  = 'chan_freq_time';
                itpc.freq    = f;
                itpc.time    = round(data_task.time{1}(1),2)+twl:twl:round(data_task.time{1}(end),2)-twl;
                itpc.grad    = data_task.grad;
                itpc.condition = i_cond;
                itpc.n_trials = n_trials;
                itpc.headpos_thresh = headposThresh;
                itpc.whichTrials = whichTrials;
                itpc.erf = nanmean(datMat,3);
                
                %% initialize output time-frequency
                itpc_data = zeros(n_channels, n_freq, n_samples_data);
                pow_data = zeros(n_channels, n_freq, n_samples_data);
                
                %% Loop trials
                for i_trial = 1:n_trials
                    % Tell which trial is processed
                    fprintf('Processing trial %d of %d\n',i_trial,n_trials)
                    
                    % cut out this trial
                    this_trial = squeeze(datMat(:,:,i_trial));
                    
                    % FFT of data (note: this doesn't change on frequency iteration)
                    fft_trial = fft(this_trial,n_samples_conv_pow2,2);
                    
                    % compute convolution for each frequency
                    for i_freq=1:n_freq
                        
                        % duplicate wavelets to match number of channel
                        wl = repmat(wavelets_fft(i_freq,:), [n_channels 1]);
                        
                        % run convolution
                        convResult = ifft(wl.*fft_trial,n_samples_conv_pow2,2);
                        convResult = convResult(:,1:n_samples_convolution); % here the extra points from the power-of-2 FFT are removed
                        convResult = convResult(:,wl_time_half+1:end-wl_time_half);
                        
                        % Put averaged data to tf-matrix
                        itpc_data(:,i_freq,:) = squeeze(itpc_data(:,i_freq,:)) + exp(1i*angle(convResult)) ./ n_trials; % calculates ITPC
                        pow_data(:,i_freq,:) = squeeze(pow_data(:,i_freq,:)) + (convResult .* conj(convResult)) ./ n_trials; % calculates power
                                
                    end %freq
                end %trial
                
                %% ITPC
                itpc_data = abs(itpc_data);
                
                %% average data acording to time window length 
                n_tw = length(itpc.time);
                twl_samples = twl * itpc.fsample;
                itpc_data_tw = zeros(n_channels, n_freq, n_tw);
                pow_data_tw = zeros(n_channels, n_freq, n_tw);
                
                for i_tw = 1:n_tw
                    itpc_data_tw(:,:,i_tw) = mean(itpc_data(:,:,(((i_tw-1)*twl_samples+1:i_tw*twl_samples)+twl_samples/2)),3);
                    pow_data_tw(:,:,i_tw)  = mean(pow_data(:,:,(((i_tw-1)*twl_samples+1:i_tw*twl_samples)+twl_samples/2)),3);
                end
                clear itpc_data pow_data
                %% write result to output structure
                itpc.itpc = itpc_data_tw;
                itpc.pow  = pow_data_tw;
                
                
                % free memory
                clear itpc_data_tw pow_data_tw fft_trial convResult wl
                
                %% save
                
                fprintf('Saving itpc output for VP_%02d_%d_wl\n',i_subj,i_ses);
                eval(sprintf('fname_out = fullfile(files_meg_itpc, [''VP_%02d_%d_Cond_%d%s_wl.mat''])',i_subj,i_ses,i_cond,whichTrials));
                save(fname_out, '-struct', 'itpc', '-v7.3');
                
                
            end %conditions
        else
            fprintf('VP_%02d_%d is/was already processed. Continue...\n',i_subj,i_ses);
        end
    end % sessions
end % subjects
disp('Done!')
clear variables
