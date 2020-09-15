%% This script simulates a phase reset and an evoked potential (CNV) at time 0 and computes low-frequency ITPC as well as power
% pinknoise.m requires Matlab 2019b or newer
% Otherwise use: Hristo Zhivomirov (2020). Pink, Red, Blue and Violet Noise Generation with Matlab (https://www.mathworks.com/matlabcentral/fileexchange/42919-pink-red-blue-and-violet-noise-generation-with-matlab)


clear all
close all

numrand = 1000;

% initialize output time-frequency
itpc_data = zeros(2, 1601,numrand);
pow_data = zeros(2, 1601,numrand);
pow_erf = zeros(2, 1601,numrand);


for i = 1:numrand
    
    disp(i)
    
    %% Parameters
    freq = 0.67; % frequency for simulating phase reset and cnv
    fs = 400;
    time = -2:1/fs:2;
    ntrials = 100;
    blc = 1; % baseline correction?
    twl_bsl = [-1.5 -0.5];
    fn_denominator = 30; % determines the degree of frequency noise
    noise =  1.5*pinknoise(length(time),ntrials);
    scale_factor = 0.5;
    phaseresetjitter = 0; % jitter in s, e.g. 0.2 for a jitter of ±200 ms around t=0
    phaseresetnoise = -phaseresetjitter:1/fs:phaseresetjitter;
    CNV_resolution = 50;
    CNV_oscillator = 0; % oscillator = 1, CNV-like ramp = 0
    
    %% Simulate phase reset and CNV
    datMat_pr = zeros(length(time),ntrials);
    datMat_cnv = zeros(length(time),ntrials);
    
    for itrial = 1:ntrials
        d = phaseresetnoise(randperm(length(phaseresetnoise),1));
        time1 = -2:1/fs:d-1/fs;
        time2 = 0:1/fs:2-d;
    
        %  ITPC
        datMat_pr(1:length(time1),itrial) = scale_factor*sin(2*pi*(freq+(rand*2-1)/fn_denominator).*time1+rand*2*pi);
        datMat_pr(length(time1)+1:end,itrial) = scale_factor*sin(2*pi*(freq+(rand*2-1)/fn_denominator).*time2);

        % evoked
        if CNV_oscillator
            % oscillator
            datMat_cnv(length(time1)+1:end,itrial) = datMat_pr(length(time1)+1:end,itrial);
        else
            % sawtooth / ramp
            datMat_cnv(length(time1)+1:end,itrial) = scale_factor*sawtooth(2*pi*(freq+(rand*2-1)/fn_denominator).*time2)+scale_factor;
            [m,ind_max] = max(datMat_cnv(length(time1):end,itrial));
            datMat_cnv(length(time1)+ind_max:length(time1)+ind_max+CNV_resolution-1,itrial) = linspace(m,0,CNV_resolution);
            datMat_cnv(length(time1)+ind_max+CNV_resolution:end,itrial) = 0;
        end
        


    end
    if ~CNV_oscillator
        datMat_cnv = -datMat_cnv;
    end
    
    %% Plot first run
    if i == 1
        h = figure('position',[100 300 1250 1000]);
        subplot(521),
        line([0 0],[-1.2 1.2],'linewidth',2,'color',[0 0 0],'linestyle','--');
        line([1.5 1.5],[-1.2 1.2],'linewidth',2,'color',[0 0 0],'linestyle','--');
        hold on
        plot(time,datMat_pr), title('phase reset'), set(gca,'ylim',[-1.2 1.2]); ylabel('Amplitude (a.u.)')
        hold off
        subplot(522),
        line([0 0],[-1.2 1.2],'linewidth',2,'color',[0 0 0],'linestyle','--');
        line([1.5 1.5],[-1.2 1.2],'linewidth',2,'color',[0 0 0],'linestyle','--');
        hold on
        plot(time,datMat_cnv), if CNV_oscillator, title('evoked'), else, title('evoked (CNV-like ramp)'),end, set(gca,'ylim',[-1.2 1.2]);
        hold off
    end
    
    %% Add noise and plot first run
    datMat_pr = datMat_pr + noise;
    datMat_cnv = datMat_cnv + noise;
    
    if i == 1
        subplot(523), plot(time,datMat_pr), set(gca,'ylim',[-7 7]); ylabel('Amplitude (a.u.)')
        subplot(524), plot(time,datMat_cnv), set(gca,'ylim',[-7 7]);
    end
    
    %% wavelet analysis
    % wavelet parameters
    wlfreq       = freq;
    wl_time      = -2.5:1/fs:2.5; %-2:1/data.fsample:2;
    wl_time_half = (length(wl_time)-1)/2;
    
    % FFT parameters (use next-power-of-2)
    n_samples_wl          = length(wl_time);
    n_samples_data        = size(time,2);
    n_samples_convolution = n_samples_wl+n_samples_data-1;
    n_samples_conv_pow2   = pow2(nextpow2(n_samples_convolution));
    f_dummy               = logspace(log10(0.5),log10(100),40); % this is to ensure equal parameters as in the wavelet analysis used in the manuscript
    fInd_dummy            = dsearchn(f_dummy',freq);
    wavelet_cycles_dummy  = logspace(log10(2),log10(10),40);
    wavelet_cycles        = wavelet_cycles_dummy(fInd_dummy);
    
    % Computation and FFT of wavelets
    wavelet = (pi*wlfreq*sqrt(pi))^-.5 * exp(2*1i*pi*wlfreq.*wl_time) .* exp(-wl_time.^2./(2*( wavelet_cycles /(2*pi*wlfreq))^2))/wlfreq;
    
    % Fourier transform of the wavelets
    wl = fft(wavelet, n_samples_conv_pow2, 2);
    
    
    % Loop trials
    for i_trial = 1:ntrials
        % Tell which trial is processed
        
        % cut out this trial
        this_trial = squeeze(datMat_pr(:,i_trial))';
        this_trial(2,:) = squeeze(datMat_cnv(:,i_trial))';
        
        % FFT of data (note: this doesn't change on frequency iteration)
        fft_trial = fft(this_trial,n_samples_conv_pow2, 2);
        
        % run convolution
        convResult = ifft(wl.*fft_trial,n_samples_conv_pow2, 2);
        convResult = convResult(:,1:n_samples_convolution); % here the extra points from the power-of-2 FFT are removed
        convResult = convResult(:,wl_time_half+1:end-wl_time_half);
        
        % Put averaged data to tf-matrix
        itpc_data(:,:,i) = itpc_data(:,:,i) + exp(1i*angle(convResult)) ./ ntrials;
        pow_data(:,:,i) = pow_data(:,:,i) + (convResult .* conj(convResult)) ./ ntrials; % calculates power
        
        
    end %trial
    
    
    
    %% Compute erf power
    
    erf_pr = mean(datMat_pr,2);
    erf_cnv = mean(datMat_cnv,2);
    % datMat_pr = datMat_pr - erf_pr;
    
    
    % cut out this trial
    
    % FFT of data (note: this doesn't change on frequency iteration)
    fft_trial = fft(erf_pr',n_samples_conv_pow2, 2);
    fft_trial(2,:) = fft(erf_cnv',n_samples_conv_pow2, 2);
    
    % run convolution
    convResult = ifft(wl.*fft_trial,n_samples_conv_pow2, 2);
    convResult = convResult(:,1:n_samples_convolution); % here the extra points from the power-of-2 FFT are removed
    convResult = convResult(:,wl_time_half+1:end-wl_time_half);
    
    % Put averaged data to tf-matrix
    pow_erf(:,:,i) = (convResult .* conj(convResult)); % calculates power
    
    
end
%% Plot
% Calculation ITPC
itpc_data = abs(itpc_data);

% Average
itpc_data = mean(itpc_data,3);
pow_data = mean(pow_data,3);
pow_erf = mean(pow_erf,3);

% baseline correction
if blc
    bsl_ind = dsearchn(time',twl_bsl');
    pow_blc = mean(pow_data(:,bsl_ind(1):bsl_ind(end)),2);
    pow_data  = (pow_data ./ repmat(pow_blc, [1 length(time)])-1) * 100;
    pow_blc_erf = mean(pow_erf(:,bsl_ind(1):bsl_ind(end)),2);
    pow_erf  = (pow_erf ./ repmat(pow_blc_erf, [1 length(time)])-1) * 100;
end

subplot(525); plot(time,pow_data(1,:),'linewidth',3), ylabel('Total power (%)'), set(gca,'ylim',[-50 50]);
subplot(527); plot(time,itpc_data(1,:),'linewidth',3), ylabel('ITPC');  set(gca,'ylim',[-0.05 1.05]);
subplot(526); plot(time,pow_data(2,:),'linewidth',3), set(gca,'ylim',[-50 50]);
subplot(528); plot(time,itpc_data(2,:),'linewidth',3), set(gca,'ylim',[-0.05 1.05]);


subplot(529); plot(time,pow_erf(1,:),'linewidth',3), ylabel('ERF power (%)'); set(gca,'ylim',[-50 4400]); xlabel('Time (s)')
subplot(5,2,10); plot(time,pow_erf(2,:),'linewidth',3), set(gca,'ylim',[-50 4400]);xlabel('Time (s)')

set(h.Children,'FontSize',20,'fontname','arial','linewidth',2)


set(h.Children(2).YLabel,'fontsize',16)
set(h.Children(5).YLabel,'fontsize',16)
set(h.Children(6).YLabel,'fontsize',16)
set(h.Children(8).YLabel,'fontsize',16)
set(h.Children(10).YLabel,'fontsize',16)
set(h.Children(1).XLabel,'fontsize',18)
set(h.Children(2).XLabel,'fontsize',18)
set(h.Children(9),'Box','on')
set(h.Children(10),'Box','on')
for i = [1 3 4 7 9]
    position = get(h.Children(i),'position');
    position(1) = position(1)-0.05;
    set(h.Children(i),'position',position)
end





%saveas(gcf,'ITPCvsCNV_simulation.tif','tiff');





