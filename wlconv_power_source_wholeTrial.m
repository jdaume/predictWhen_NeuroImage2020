%% Compute power in source space, wholeTrial analysis (all windows)

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
whichTrials = '_stratSubjCorr'; % all = '_all', correct = '_corr', stratified, subjectively correct = '_stratSubjCorr'
whichGrid = 'medium';  
whichFilter = 'beamformer'; 

% files output
eval(sprintf('files_in = fullfile(files_data, ''data_meg_%s_wholeTrial/'');',TFAmethod));
eval(sprintf('files_filter = fullfile(files_data, ''data_meg_source/%s/%s/%s/filter_A'');',TFAmethod,whichFilter,whichGrid));
eval(sprintf('files_power = fullfile(files_data, ''data_meg_source/%s/%s/%s/power'');',TFAmethod,whichFilter,whichGrid));
eval(sprintf('files_logfiles = fullfile(files_data, ''logfiles/wlconv_wholeTrial_sourcepow/%s_%s'');',whichGrid,whichFilter));
if ~exist(files_power,'dir'), mkdir(files_power); end
if ~exist(files_filter,'dir'), mkdir(files_filter); end
if ~exist(files_logfiles,'dir'), mkdir(files_logfiles); end


%% Loop subjects


for i_subj = subjects
    for i_ses = sessions
        
        if i_subj == 18 && i_ses == 2
            continue;
        end
        
        if ~exist(sprintf([files_logfiles '/%s%s_VP%02d_%d_source_pow_processing.txt'],TFAmethod,whichTrials,i_subj,i_ses),'file')
            
            system(['touch ' files_logfiles sprintf('/%s%s_VP%02d_%d_source_pow_processing.txt',TFAmethod,whichTrials,i_subj,i_ses)]);
            
            clear sa L A
            % Load/compute spatial filter
            
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
                    eval(sprintf('cs_path = fullfile(files_in, [''VP_%02d_%d_Cond_%d%s_cs.mat'']);',i_subj,i_ses,itask,whichTrials));
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
                clear cs_avg
                disp('Saving common beamformer filter...')
                eval(sprintf('save(fullfile(files_filter, ''VP_%02d_%d_A%s.mat''),''A'',''-v7.3'')',i_subj,i_ses,whichTrials));
            end
            
            
            
            %% trial loop
            for i_task = 1:3
                
                %
                clear cs_cond
                if ~exist('cs_data','var')
                    clear cs
                    fprintf('Loading VP_%02d_%d_Task_%d%s_cs.mat...\n',i_subj,i_ses,i_task,whichTrials)
                    eval(sprintf('wl_path = fullfile(files_in, [''VP_%02d_%d_Cond_%d%s_cs.mat'']);',i_subj,i_ses,i_task,whichTrials));
                    cs = load(wl_path);
                    xspectrm = cat(4,cs.cs_bsl, cs.cs_stimMov, cs.cs_stimOff, cs.cs_stimReapp);
                    cs = rmfield(cs,{'cs_bsl','cs_stimMov','cs_stimOff','cs_stimReapp'});
                                   
                    cs_cond = xspectrm;
                    
                    n_freq = length(cs.freq);
                    n_tw   = length(cs.time);
                    n_grid = size(A,2);
                else
                    cs_cond = cs_data(:,:,:,:,i_task);
                end
                
                % output variable
                clear source_pow pow
                source_pow = cs;
                source_pow = rmfield(source_pow,{'label','dimord'});
                source_pow.dimord = 'vox_freq_time';
                source_pow.grid = whichGrid;
                source_pow.filter = whichFilter;
                source_pow.n_voxel = n_grid;
                
                %% Compute power in source space
                pow = zeros(n_grid,n_freq,n_tw);
                
                fprintf('Computing source power for all frequencies and all time windows in task %d...\n',i_task);
                
                for i_freq = 1:n_freq % frequency
                    for i_tw = 1:n_tw % time window
                        pow(:,i_freq,i_tw) = abs(diag(A(:,:,i_freq)'*cs_cond(:,:,i_freq,i_tw)*A(:,:,i_freq)));
                    end
                end

                source_pow.pow = pow;
                
                %% save output
                eval(sprintf('fname_out = fullfile(files_power, [''VP_%02d_%d_Cond_%d%s.mat'']);',i_subj,i_ses,i_task,whichTrials));
                fprintf('Saving power file to %s\n',fname_out);
                save(fname_out, '-struct', 'source_pow', '-v7.3');
                
            end
        else
            fprintf('VP_%02d_%d is/was already processed. Continue...\n',i_subj,i_ses);
        end
    end
end

disp('Done!')

