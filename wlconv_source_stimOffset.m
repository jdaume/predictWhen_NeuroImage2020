%% Computes ITPC (and power) for the enlarged window around stimulus disappearance in source space

% JD 2020


clear variables

% Files & Folders
addpath('.../fieldtrip-20170607');  % change to personal location of fieldtrip
addpath('.../hhtb');  % change to personal location of Hamburg MEG/EEG toolbox (G. Nolte)

startup_hhtb;
ft_defaults;

files_data = fullfile(pwd);
files_sa = fullfile(files_data, 'data_sa');

%% Switches

subjects = 1:23;
sessions = 1:2;

n_subjects = length(subjects);

TFAmethod = 'wlconv';
whichGrid = 'medium'; % medium cortex3000 cortex5124
whichFilter = 'beamformer'; % beamformer eloreta

headposThresh = 5; % kick out all trials above this threshold in any coil in mm (trials are saved in plot_and_cutoff_headpos); 0 = no cutoff!
correct = 1; % 1 = correct trials, 0 = all trials
stratified = 1; % 1 = stratified number of trials, 0 = all (correct) trials
subjCorr = 1; % only works when correct = 1; 1 = subjectively correct trials, 0 = objectively correct trials
twl = 0.1; % time window length for output in s

toilim = [-1 2]; % time window of interest around stimulus disappearance in s


if correct && stratified && ~subjCorr
    whichTrials = '_stratCorr';
elseif correct && stratified && subjCorr
    whichTrials = '_stratSubjCorr';
elseif correct && ~stratified
    whichTrials = '_corr';
else
    whichTrials = '_all';
end

% files output
files_headpos = fullfile(pwd, 'data_headpos');
files_strat = fullfile(pwd, 'data_stratTrials');
files_meg_preproc = fullfile(pwd, 'data_meg_preprocessed');
files_cs = fullfile(pwd, 'data_meg_wlconv_wholeTrial');
eval(sprintf('files_itpc = fullfile(pwd, ''data_meg_source/%s/%s/%s/itpc_stimOffset'');',TFAmethod,whichFilter,whichGrid));
eval(sprintf('files_pow = fullfile(pwd, ''data_meg_source/%s/%s/%s/pow_stimOffset'');',TFAmethod,whichFilter,whichGrid));
eval(sprintf('files_logfiles = fullfile(pwd, ''logfiles/wlconv_sourceitpc/%s/%s_%s'');',whichWlFreqs,whichGrid,whichFilter));
eval(sprintf('files_filter = fullfile(pwd, ''data_meg_source/%s/%s/%s/filter_A'');',TFAmethod,whichFilter,whichGrid));
files_psych = fullfile(pwd, 'data_psychThresh'); % output from psychometric function

if ~exist(files_itpc,'dir'), mkdir(files_itpc); end
if ~exist(files_pow,'dir'), mkdir(files_pow); end
if ~exist(files_logfiles,'dir'), mkdir(files_logfiles); end
if ~exist(files_filter,'dir'), mkdir(files_filter); end

%% Trialinfo and stratified trials
trialinfo_predictWhen; % returns variable 'col' containing info about which column represents which info in trialinfo matrix

if correct && stratified && ~subjCorr
    eval(sprintf('load(fullfile(files_strat,''stratTrials_corr_headpos%dmm.mat''));',headposThresh));
elseif correct && stratified && subjCorr
    eval(sprintf('load(fullfile(files_strat,''stratTrials_subjCorr_headpos%dmm.mat''));',headposThresh));
end

load(fullfile(files_psych,'psychThresh_bothSessions.mat'));

%% Loop subjects


for i_subj = subjects
    for i_ses = sessions
        
        if i_subj == 18 && i_ses == 2
            continue;
        end
        
        if ~exist(sprintf([files_logfiles '/%s_VP%02d_%d_source_itpc_processing.txt'],whichTrials,i_subj,i_ses),'file')
            
            system(['touch ' files_logfiles sprintf('/%s_VP%02d_%d_source_itpc_processing.txt',whichTrials,i_subj,i_ses)]); % creates a logfile 
            
            %% load or compute spatial filter
            clear sa L A
            
            if eval(sprintf('exist(fullfile(files_filter, ''VP_%02d_%d_A%s.mat''))',i_subj,i_ses,whichTrials))
                disp('Loading common beamformer filter...')
                eval(sprintf('load(fullfile(files_filter, ''VP_%02d_%d_A%s.mat''))',i_subj,i_ses,whichTrials));
            else
                fprintf('Loading VP_%02d_%d_sa.mat...\n',i_subj,i_ses);
                eval(sprintf('sa_path = fullfile(files_sa, [''VP_%02d_%d_sa.mat'']);',i_subj,i_ses));
                sa = load(sa_path);
                fprintf('Computing leadfield for grid %s\n',whichGrid)
                eval(sprintf('L = grid2L(sa.grid_%s_indi,sa.fp_indi);',whichGrid));
                
                % Loading all conditions for common filter
                disp('Loading all tasks for computation of common beamformer filter');
                clear cs_data cs_avg
                for itask = 1:3
                    clear cs
                    fprintf('Loading VP_%02d_%d_Cond_%d%s_cs.mat...\n',i_subj,i_ses,itask,whichTrials)
                    eval(sprintf('cs_path = fullfile(files_cs, [''VP_%02d_%d_Cond_%d%s_cs.mat'']);',i_subj,i_ses,itask,whichTrials));
                    cs = load(cs_path);
                    
                    xspectrm = cat(4,cs.cs_bsl, cs.cs_stimMov, cs.cs_stimOff, cs.cs_stimReapp);
                    cs = rmfield(cs,{'cs_bsl','cs_stimMov','cs_stimOff','cs_stimReapp'});
                    
                    %% Initialize and write to cs_data
                    if itask == 1
                        chanInds = sa.inds(ismember(sa.sens_indi.label(strncmp(sa.sens_indi.label,'M',1)),cs.label));
                        L = L(chanInds,:,:); % take only the channels that are relevant here
                        n_freq = size(xspectrm,3);
                        n_tw   = size(xspectrm,4);
                        n_chan = size(xspectrm,1);
                        n_grid = size(L,2);
                        cs_data = zeros([size(xspectrm),3]);
                    end
                    cs_data(:,:,:,:,itask) = xspectrm;
                end
                cs_avg = squeeze(mean(squeeze(mean(cs_data,5)),4));
                disp('Computing Beamformer filter for all frequencies...')
                A = zeros(n_chan,n_grid,n_freq);
                for ifreq = 1:n_freq
                    filt = beamformer_jh(cs_avg(:,:,ifreq),L);
                    A(:,:,ifreq) = filt';
                end
                clear cs_avg cs_data
                disp('Saving common beamformer filter...')
                eval(sprintf('save(fullfile(files_filter, ''VP_%02d_%d_A%s.mat''),''A'',''-v7.3'')',i_subj,i_ses,whichTrials));
            end
              
            %% get MEG data filename for this subject
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
            
            data = ft_preprocessing(cfg,data);
            
            %% redefine trials to equal length
            cfg = [];
            cfg.toilim = [-2 2.2]; % first only roughly and longer that time of interest in order to avoid edge artifacts later for the itpc analysis
            data = ft_redefinetrial(cfg, data);
            
            %% Wavelet convolution parameters
            
            % wavelet parameters
            min_freq     = 0.5;
            max_freq     = 100;
            n_freq       = 40;
            f            = logspace(log10(min_freq),log10(max_freq),n_freq);
            wl_time      = -2.5:1/data.fsample:2.5;
            wl_time_half = (length(wl_time)-1)/2;
            
            f_cutoffInd  = dsearchn(f',5); % only compute wavelets up to 5 Hz to save space and time
            f            = f(1:f_cutoffInd);
            
            % FFT parameters (use next-power-of-2)
            n_samples_wl          = length(wl_time);
            n_samples_data        = size(data.time{1},2);
            n_samples_convolution = n_samples_wl+n_samples_data-1;
            n_samples_conv_pow2   = pow2(nextpow2(n_samples_convolution));
            wavelet_cycles        = logspace(log10(2),log10(10),n_freq);
            wavelet_cycles        = wavelet_cycles(1:f_cutoffInd);
            
            n_freq = f_cutoffInd;
            
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
                n_vox          = size(A,2);
                  
                %% Write into 3D matrix
                
                clear datMat
                
                % convert cell array to mat, channels x samples x trials
                disp('Writing data into 3D matrix...')
                datMat = cat(3,data_task.trial{:});
                
                %% create FT output structure for itpc
                
                %free memory
                clear srcitpc
                
                srcitpc         = [];
                eval(sprintf('srcitpc.subj    = ''VP_%02d'';',i_subj));
                srcitpc.session = i_ses;
                srcitpc.method  = 'wlconv';
                srcitpc.label   = data_task.label;
                srcitpc.dimord  = 'vox_freq_time';
                srcitpc.freq    = f;
                srcitpc.time    = round(toilim(1),2)+twl:twl:round(toilim(end),2)-twl;
                srcitpc.condition = i_cond;
                srcitpc.n_trials = n_trials;
                srcitpc.headpos_thresh = headposThresh;
                srcitpc.whichTrials = whichTrials;
                srcitpc.grid = whichGrid;
                srcitpc.filter = whichFilter;
                srcitpc.n_voxel = n_vox;
                
                %free memory
                clear srcpow
                
                srcpow         = srcitpc;
                
                %% initialize output time-frequency
                itpc_data     = zeros(n_vox, n_freq, n_samples_data);
                pow_data     = zeros(n_vox, n_freq, n_samples_data);
                
                %% Loop trials
                for i_trial = 1:n_trials
                    % Tell which trial is processed
                    fprintf('Processing trial %d of %d\n',i_trial,n_trials)
                    
                    % cut out this trial
                    this_trial = squeeze(datMat(:,:,i_trial));
                    
                    % FFT of data (note: this doesn't change on frequency iteration)
                    fft_trial = fft(this_trial,n_samples_conv_pow2,2);
                    clear this_trial
                    
                    % compute convolution for each frequency
                    for i_freq=1:n_freq
                        
                        clear convResult convResult_source
                        
                        % duplicate wavelets to match number of channel
                        wl = repmat(wavelets_fft(i_freq,:), [n_channels 1]);
                        
                        % run convolution
                        convResult = ifft(wl.*fft_trial,n_samples_conv_pow2,2);
                        convResult = convResult(:,1:n_samples_convolution); % here the extra points from the power-of-2 FFT are removed
                        convResult = convResult(:,wl_time_half+1:end-wl_time_half);
                        
                        % project result to source space
                        convResult_source = complex(zeros(n_vox,n_samples_data));
                        for i_vox = 1:n_vox
                            convResult_source(i_vox,:) = A(:,i_vox,i_freq)'*convResult;
                        end
                        
                        % Put averaged data to tf-matrix
                        itpc_data(:,i_freq,:) = squeeze(itpc_data(:,i_freq,:)) + exp(1i*angle(convResult_source)) ./ n_trials;
                        pow_data(:,i_freq,:) = squeeze(pow_data(:,i_freq,:)) + (convResult_source .* conj(convResult_source)) ./ n_trials;
                        
                    end %freq
                end %trial
                
                %% cut data to toilim (the shorter time window of interest defined in the beginning, to save space)
                toi = dsearchn(data_task.time{1}',toilim');
                itpc_data = itpc_data(:,:,toi(1):toi(2));
                pow_data = pow_data(:,:,toi(1):toi(2));
                
                %% Calculation ITPC
                itpc_data = abs(itpc_data);
                
                %% average data acording to time window length
                n_tw = length(srcitpc.time);
                twl_samples = twl * data_task.fsample;
                itpc_data_tw = zeros(n_vox, n_freq, n_tw);
                pow_data_tw = zeros(n_vox, n_freq, n_tw);
                
                for i_tw = 1:n_tw
                    itpc_data_tw(:,:,i_tw) = mean(itpc_data(:,:,(((i_tw-1)*twl_samples+1:i_tw*twl_samples)+twl_samples/2)),3);
                    pow_data_tw(:,:,i_tw) = mean(pow_data(:,:,(((i_tw-1)*twl_samples+1:i_tw*twl_samples)+twl_samples/2)),3);
                end
                clear itpc_data pow_data
                %% write result to output structure
                srcitpc.itpc = itpc_data_tw;
                srcpow.pow = pow_data_tw;
                
                % free memory
                clear itpc_data_tw fft_trial convResult wl pow_data_tw
                
                %% save
                
                fprintf('Saving srcitpc output for VP_%02d_%d\n',i_subj,i_ses);
                eval(sprintf('fname_out = fullfile(files_itpc, [''VP_%02d_%d_Cond_%d%s.mat''])',i_subj,i_ses,i_cond,whichTrials));
                save(fname_out, '-struct', 'srcitpc', '-v7.3');
                
                fprintf('Saving srcpow output for VP_%02d_%d\n',i_subj,i_ses);
                eval(sprintf('fname_out = fullfile(files_pow, [''VP_%02d_%d_Cond_%d%s.mat''])',i_subj,i_ses,i_cond,whichTrials));
                save(fname_out, '-struct', 'srcpow', '-v7.3');
                
                
            end %conditions
        else
            fprintf('VP_%02d_%d is/was already processed. Continue...\n',i_subj,i_ses);
        end
    end % sessions
end % subjects
disp('Done!')


