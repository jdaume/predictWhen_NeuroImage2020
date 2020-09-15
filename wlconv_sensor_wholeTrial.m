%% Complex morelet wavelet convolution script for all events (windows) in the trial
% this script computes power and itpc and cross-spectra for four time windows:
% baseline, stimulus movement, stimulus disappearance, stimulus
% reappearance
% trials are not equally long, so that every trial needs to be cut
% individually

% JD 2020

clear variables

%% Fieldtrip
addpath('.../fieldtrip-20170607'); %% change to personal location of fieldtrip
ft_defaults

%%
headposThresh = 5; % kick out all trials above this threshold in any coil in mm (trials are saved in plot_and_cutoff_headpos); 0 = no cutoff!
correct = 1; % 1 = correct trials, 0 = all trials
stratified = 1; % 1 = stratified number of trials, 0 = all (correct) trials
subjCorr = 1; % only works when correct = 1; 1 = subjectively correct trials, 0 = objectively correct trials

save_pow = 1;
save_itpc = 1;
save_cs = 1;

twl = 0.1; % time window length for averaging samples into
twl_bsl = 0.8; % time window length for baseline in s (right before stimulus movement onset)
twl_stimMov = 1; % time window length for stimulus movement (start with onset of stimulus)
twl_stimOff = 1.3; % time window length for stimulus offset (starts with twOnset_stimOff)
twOnset_stimOff = -0.3; % offset for time window for stimulus offset (relativ to offset of stimulus, which is 0)
twl_stimReapp = 0.8; % time window length for stimulus reappearance (starts with twOnset_Reapp)
twOnset_Reapp = -0.3; % offset for time window for stimulus reapp relative to actual reappearance


if correct && stratified && ~subjCorr
    whichTrials = '_stratCorr';   
elseif correct && stratified && subjCorr
    whichTrials = '_stratSubjCorr';
elseif correct && ~stratified
    whichTrials = '_corr';
else
    whichTrials = '_all';
end

%% files & folders

files_meg_preproc = fullfile(pwd, 'data_meg_preprocessed');
files_meg_out = fullfile(pwd, 'data_meg_wlconv_wholeTrial');
files_headpos = fullfile(pwd, 'data_headpos');
files_logfiles = fullfile(pwd, 'logfiles/wlconv_wholeTrial');
files_strat = fullfile(pwd, 'data_stratTrials');
if ~exist(files_meg_out,'dir'), mkdir(files_meg_out); end
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
        
        if ~exist(sprintf([files_logfiles '/VP%02d_%d%s_processing.txt'],i_subj,i_ses,whichTrials),'file') % continue if already exists (another process is or finished working on this already)
            
            system(['touch ' files_logfiles sprintf('/VP%02d_%d%s_processing.txt',i_subj,i_ses,whichTrials)]); 
            
            %free memory
            clear data
            
            % get MEG data filename for this subject
            fprintf('Loading MEG-data for VP_%02d_%d\n',i_subj,i_ses);
            eval(sprintf('this_meg_data = fullfile(files_meg_preproc, [''VP_%02d_%d_compRejected.mat'']);',i_subj,i_ses));
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
            
            % get time for moving stimulus and write into trialinfo
            col.tMovStim = 18;
            data.trialinfo(:,col.tMovStim) = ((data.trialinfo(:,col.startpos)+20-245)/60)/5;
            
            %% find longest trial (for zero padding)
            n_trials_all = size(data.trialinfo,1);
            n_samples_trials = zeros(n_trials_all,1);
            for itrial = 1:n_trials_all
                n_samples_trials(itrial) = size(data.time{itrial},2);
            end
            
            %% Wavelet convolution parameters
            
            % wavelet parameters
            min_freq     = 0.5;
            max_freq     = 100;
            n_freq       = 40; % number of frequencies
            f            = logspace(log10(min_freq),log10(max_freq),n_freq);  % log-spacing of frequencies
            wl_time      = -2.5:1/data.fsample:2.5; % wavelet window
            wl_time_half = (length(wl_time)-1)/2;
            
            % FFT parameters (use next-power-of-2)
            n_samples_wl          = length(wl_time);
            n_samples_data        = max(n_samples_trials);
            n_samples_convolution = n_samples_wl+n_samples_data-1;
            n_samples_conv_pow2   = pow2(nextpow2(n_samples_convolution));
            wavelet_cycles        = logspace(log10(2),log10(10),n_freq); % log-spacing for number of cycles for each wavelet
            
            
            % Computation and FFT of wavelets
            fprintf('Computation and fft of wavelets for all %d specified frequencies...\n',n_freq)
            % compute wavelets for each frequency
            wavelets = zeros(n_freq, n_samples_wl);
            for fi = 1:n_freq
                wavelets(fi,:) = (pi*f(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*f(fi).*wl_time) .* exp(-wl_time.^2./(2*( wavelet_cycles(fi) /(2*pi*f(fi)))^2))/f(fi);
            end
            
            % Fourier transformation of the wavelets
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
                
                %% Write data in 3D matric
                
                clear datMat
                
                % convert cell array to mat, channels x samples x trials
                % this may take a while due to looping... necessary since
                % all trials have different length
                disp('Writing data into 3D matrix...')
                datMat = zeros(n_channels,max(n_samples_trials),n_trials);
                for itrial = 1:n_trials
                    datMat(:,1:size(data_task.trial{itrial},2),itrial) = data_task.trial{itrial};
                end
                
                %% create FT output structure for pow
                
                % power
                if save_pow
                    clear data_tf
                    
                    data_tf         = [];
                    eval(sprintf('data_tf.subj    = ''VP_%02d'';',i_subj));
                    data_tf.session = i_ses;
                    data_tf.method  = 'wlconv';
                    data_tf.label   = data_task.label;
                    data_tf.dimord  = 'chan_freq_time';
                    data_tf.freq    = f;
                    data_tf.time_bsl = (-(twl_bsl+twl_stimMov):twl:-twl_stimMov-twl)+twOnset_stimOff;
                    data_tf.time_stimMov = (-twl_stimMov:twl:0-twl)+twOnset_stimOff;
                    data_tf.time_stimOff = (0:twl:twl_stimOff-twl)+twOnset_stimOff;
                    data_tf.time_stimReapp = (twl_stimOff:twl:twl_stimOff+twl_stimReapp-twl)+twOnset_stimOff;
                    data_tf.time = (-(twl_bsl+twl_stimMov):twl:twl_stimOff+twl_stimReapp-twl)+twOnset_stimOff;
                    data_tf.grad    = data_task.grad;
                    data_tf.condition = i_cond;
                    data_tf.n_trials = n_trials;
                    data_tf.headpos_thresh = headposThresh;
                    data_tf.whichTrials = whichTrials;
                end
                
                
                % itpc
                if save_itpc
                    clear itpc
                    
                    itpc         = [];
                    eval(sprintf('itpc.subj    = ''VP_%02d'';',i_subj));
                    itpc.session = i_ses;
                    itpc.method  = 'wlconv';
                    itpc.label   = data_task.label;
                    itpc.dimord  = 'chan_freq_time';
                    itpc.freq    = f;
                    itpc.time_bsl = (-(twl_bsl+twl_stimMov):twl:-twl_stimMov-twl)+twOnset_stimOff;
                    itpc.time_stimMov = (-twl_stimMov:twl:0-twl)+twOnset_stimOff;
                    itpc.time_stimOff = (0:twl:twl_stimOff-twl)+twOnset_stimOff;
                    itpc.time_stimReapp = (twl_stimOff:twl:twl_stimOff+twl_stimReapp-twl)+twOnset_stimOff;
                    itpc.time = (-(twl_bsl+twl_stimMov):twl:twl_stimOff+twl_stimReapp-twl)+twOnset_stimOff;
                    itpc.grad    = data_task.grad;
                    itpc.condition = i_cond;
                    itpc.n_trials = n_trials;
                    itpc.headpos_thresh = headposThresh;
                    itpc.whichTrials = whichTrials;
                end
                
                
                % cross-spectrum
                if save_cs
                    clear cs
                    
                    cs         = [];
                    eval(sprintf('cs.subj= ''VP_%02d'';',i_subj));
                    cs.session = i_ses;
                    cs.method  = 'wlconv';
                    cs.label   = data_task.label;
                    cs.dimord  = 'chan_chan_freq_time';
                    cs.freq    = f;
                    cs.time_bsl = (-(twl_bsl+twl_stimMov):twl:-twl_stimMov-twl)+twOnset_stimOff;
                    cs.time_stimMov = (-twl_stimMov:twl:0-twl)+twOnset_stimOff;
                    cs.time_stimOff = (0:twl:twl_stimOff-twl)+twOnset_stimOff;
                    cs.time_stimReapp = (twl_stimOff:twl:twl_stimOff+twl_stimReapp-twl)+twOnset_stimOff;
                    cs.time = (-(twl_bsl+twl_stimMov):twl:twl_stimOff+twl_stimReapp-twl)+twOnset_stimOff;
                    cs.condition = i_cond;
                    cs.n_trials = n_trials;
                    cs.headpos_thresh = headposThresh;
                    cs.whichTrials = whichTrials;
                end
                
                %% initialize output time-frequency
                
                fs = data_task.fsample;
                if save_pow
                    pow_bsl = zeros(n_channels, n_freq, twl_bsl*fs);
                    pow_stimMov   = zeros(n_channels, n_freq, twl_stimMov*fs);
                    pow_stimOff   = zeros(n_channels, n_freq, twl_stimOff*fs);
                    pow_stimReapp = zeros(n_channels, n_freq, twl_stimReapp*fs);
                end
                if save_itpc
                    itpc_bsl = zeros(n_channels, n_freq, twl_bsl*fs);
                    itpc_stimMov   = zeros(n_channels, n_freq, twl_stimMov*fs);
                    itpc_stimOff   = zeros(n_channels, n_freq, twl_stimOff*fs);
                    itpc_stimReapp = zeros(n_channels, n_freq, twl_stimReapp*fs);
                end
                
                if save_cs
                    xspectr_bsl = zeros(n_channels,n_channels, n_freq, round(twl_bsl/twl));
                    xspectr_stimMov   = zeros(n_channels,n_channels, n_freq,  round(twl_stimMov/twl));
                    xspectr_stimOff   = zeros(n_channels,n_channels, n_freq,  round(twl_stimOff/twl));
                    xspectr_stimReapp = zeros(n_channels,n_channels, n_freq,  round(twl_stimReapp/twl));
                end
                
                %% Loop trials
                for i_trial = 1:n_trials
                    % Tell which trial is processed
                    fprintf('Processing trial %d of %d\n',i_trial,n_trials)
                    
                    % cut out this trial
                    this_trial = squeeze(datMat(:,:,i_trial));
                    
                    % FFT of data (note: this doesn't change on frequency iteration)
                    fft_trial = fft(this_trial,n_samples_conv_pow2,2);
                    
                    % find time indices for the different time windows in
                    % each trial
                    ind_tstimMov  = dsearchn(data_task.time{i_trial}',-data_task.trialinfo(i_trial,col.tMovStim))-twl/2*fs+1; % subtract half twl so that time window is center on the actual time ind
                    ind_tbsl      = ind_tstimMov-twl_bsl*fs;
                    ind_tstimOff  = dsearchn(data_task.time{i_trial}',twOnset_stimOff)-twl/2*fs+1;
                    ind_tstimReapp = dsearchn(data_task.time{i_trial}',twOnset_Reapp)+round((1.5+data_task.trialinfo(i_trial,col.timingDiff)/60)*fs)-twl/2*fs+1;
                    
                    tbsl       = ind_tbsl:ind_tbsl+twl_bsl*fs-1;
                    tstimMov   = ind_tstimMov:ind_tstimMov+twl_stimMov*fs-1;
                    tstimOff   = ind_tstimOff:ind_tstimOff+twl_stimOff*fs-1;
                    tstimReapp = ind_tstimReapp:ind_tstimReapp+twl_stimReapp*fs-1;
                    
                    % compute convolution for each frequency
                    for i_freq=1:n_freq
                        
                        % duplicate wavelets to match number of channel
                        wl = repmat(wavelets_fft(i_freq,:), [n_channels 1]);
                        
                        % run convolution
                        convResult = ifft(wl.*fft_trial,n_samples_conv_pow2,2);
                        convResult = convResult(:,1:n_samples_convolution); % here the extra points from the power-of-2 FFT are removed
                        convResult = convResult(:,wl_time_half+1:end-wl_time_half);
                        
                        if save_pow
                            pow_bsl(:,i_freq,:)       = squeeze(pow_bsl(:,i_freq,:)) + (convResult(:,tbsl) .* conj(convResult(:,tbsl))) ./ n_trials; % computes power and adds to previous trials
                            pow_stimMov(:,i_freq,:)   = squeeze(pow_stimMov(:,i_freq,:)) + (convResult(:,tstimMov) .* conj(convResult(:,tstimMov))) ./ n_trials; 
                            pow_stimOff(:,i_freq,:)   = squeeze(pow_stimOff(:,i_freq,:)) + (convResult(:,tstimOff) .* conj(convResult(:,tstimOff))) ./ n_trials; 
                            pow_stimReapp(:,i_freq,:) = squeeze(pow_stimReapp(:,i_freq,:)) + (convResult(:,tstimReapp) .* conj(convResult(:,tstimReapp))) ./ n_trials; 
                        end
                        
                        if save_itpc
                            itpc_bsl(:,i_freq,:)       = squeeze(itpc_bsl(:,i_freq,:)) + exp(1i*angle(convResult(:,tbsl))) ./ n_trials; % computes ITPC and adds to previous trials
                            itpc_stimMov(:,i_freq,:)   = squeeze(itpc_stimMov(:,i_freq,:)) + exp(1i*angle(convResult(:,tstimMov))) ./ n_trials;
                            itpc_stimOff(:,i_freq,:)   = squeeze(itpc_stimOff(:,i_freq,:)) + exp(1i*angle(convResult(:,tstimOff))) ./ n_trials;
                            itpc_stimReapp(:,i_freq,:) = squeeze(itpc_stimReapp(:,i_freq,:)) + exp(1i*angle(convResult(:,tstimReapp))) ./ n_trials;
                        end
                        
                        if save_cs
                            for i_tw = 1:length(cs.time_bsl)
                                xspectr_bsl(:,:,i_freq,i_tw) = squeeze(xspectr_bsl(:,:,i_freq, i_tw)) + ...
                                    (convResult(:,tbsl(round((i_tw-1)*twl*fs+1:i_tw*twl*fs))) * convResult(:,tbsl(round((i_tw-1)*twl*fs+1:i_tw*twl*fs)))') ./ n_trials; % computes CS and adds to previous trials
                            end
                            for i_tw = 1:length(cs.time_stimMov)
                                xspectr_stimMov(:,:,i_freq,i_tw) = squeeze(xspectr_stimMov(:,:,i_freq, i_tw)) + ...
                                    (convResult(:,tstimMov(round((i_tw-1)*twl*fs+1:i_tw*twl*fs))) * convResult(:,tstimMov(round((i_tw-1)*twl*fs+1:i_tw*twl*fs)))') ./ n_trials;
                            end
                            for i_tw = 1:length(cs.time_stimOff)
                                xspectr_stimOff(:,:,i_freq,i_tw) = squeeze(xspectr_stimOff(:,:,i_freq, i_tw)) + ...
                                    (convResult(:,tstimOff(round((i_tw-1)*twl*fs+1:i_tw*twl*fs))) * convResult(:,tstimOff(round((i_tw-1)*twl*fs+1:i_tw*twl*fs)))') ./ n_trials;
                            end
                            for i_tw = 1:length(cs.time_stimReapp)
                                xspectr_stimReapp(:,:,i_freq,i_tw) = squeeze(xspectr_stimReapp(:,:,i_freq, i_tw)) + ...
                                    (convResult(:,tstimReapp(round((i_tw-1)*twl*fs+1:i_tw*twl*fs))) * convResult(:,tstimReapp(round((i_tw-1)*twl*fs+1:i_tw*twl*fs)))') ./ n_trials;
                            end
                        end
                        
                        
                    end %freq
                end %trial
                
                %% Combining and averaging sample-wise data into longer time windows
                twl_samples = twl * data_task.fsample;
                if save_pow
                    n_tw = length(data_tf.time);
                    pow = cat(3,pow_bsl,pow_stimMov,pow_stimOff,pow_stimReapp);
                    pow_tw       = zeros(n_channels, n_freq, n_tw);
                    for i_tw = 1:n_tw
                        pow_tw(:,:,i_tw)  = mean(pow(:,:,(((i_tw-1)*twl_samples+1:i_tw*twl_samples))),3);
                    end
                end
                if save_itpc
                    n_tw = length(itpc.time);
                    itpc_data = cat(3,itpc_bsl,itpc_stimMov,itpc_stimOff,itpc_stimReapp);
                    itpc_data = abs(itpc_data);
                    itpc_data_tw = zeros(n_channels, n_freq, n_tw);
                    for i_tw = 1:n_tw
                        itpc_data_tw(:,:,i_tw)  = mean(itpc_data(:,:,(((i_tw-1)*twl_samples+1:i_tw*twl_samples))),3);
                    end
                end

                clear pow_bsl pow_stimMov pow_stimOff pow_stimReapp  itpc_data itpc_bsl itpc_stimMov itpc_stimOff itpc_stimReapp fft_trial convResult wl 
                %% write result to output structure
                if save_pow
                    data_tf.pow  = pow_tw;
                end
                if save_itpc
                    itpc.itpc = itpc_data_tw;
                end
                if save_cs
                    cs.cs_bsl       = xspectr_bsl;
                    cs.cs_stimMov   = xspectr_stimMov;
                    cs.cs_stimOff   = xspectr_stimOff;
                    cs.cs_stimReapp = xspectr_stimReapp;
                end
                % free memory
                clear itpc_data_tw pow_tw xspectr_bsl xspectr_stimMov xspectr_stimOff xspectr_stimReapp
                
                %% save
                
                if save_pow
                    fprintf('Saving data_tf output for VP_%02d_%d_pow\n',i_subj,i_ses);
                    eval(sprintf('fname_out = fullfile(files_meg_out, [''VP_%02d_%d_Cond_%d%s_pow.mat''])',i_subj,i_ses,i_cond,whichTrials));
                    save(fname_out, '-struct', 'data_tf', '-v7.3');
                end
                
                if save_itpc
                    fprintf('Saving itpc output for VP_%02d_%d_itpc\n',i_subj,i_ses);
                    eval(sprintf('fname_out = fullfile(files_meg_out, [''VP_%02d_%d_Cond_%d%s_itpc.mat''])',i_subj,i_ses,i_cond,whichTrials));
                    save(fname_out, '-struct', 'itpc', '-v7.3');
                end
                
                if save_cs
                    fprintf('Saving cs output for VP_%02d_%d_cs\n',i_subj,i_ses);
                    eval(sprintf('fname_out = fullfile(files_meg_out, [''VP_%02d_%d_Cond_%d%s_cs.mat''])',i_subj,i_ses,i_cond,whichTrials));
                    save(fname_out, '-struct', 'cs', '-v7.3');
                end
                
            end %conditions
        else
            fprintf('VP_%02d_%d is/was already processed. Continue...\n',i_subj,i_ses);
        end
    end % sessions
end % subjects
disp('Done!')
clear variables
