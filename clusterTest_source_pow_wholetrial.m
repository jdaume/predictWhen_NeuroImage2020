%% Cluster permutation tests of power in source space (wholeTrial windows analysis)
% JD 2020

clear variables

% Files & Folders
addpath('.../fieldtrip-20170607');  % change to personal location of fieldtrip
addpath('.../hhtb');  % change to personal location of Hamburg MEG/EEG toolbox (G. Nolte)

startup_hhtb;
ft_defaults;

load sa_template_JD;
sa = sa_template; clear sa_template

%% Switches

subjects = 1:23;
sessions = 1:2;
n_sessions = length(sessions);

TFAmethod = 'wlconv';
whichTrials = '_stratSubjCorr'; 
whichGrid = 'medium'; 
whichFilter = 'beamformer'; 

blc = 1; % baseline correction?
twl_bsl = 0.3; %baseline length in s, starts with first time sample in baseline
bsl_offset = 4; % in number of time windows
twl = 0.1; % length of time shifts

bsl_ind = bsl_offset:twl_bsl/twl+bsl_offset;


eval(sprintf('grid = sa.grid_%s;',whichGrid));
n_vox = size(grid,1);
[xrange,yrange,zrange] = grid2cube(grid); 
xgrid = length(xrange);
ygrid = length(yrange);
zgrid = length(zrange);
n_voxel_cube = xgrid*ygrid*zgrid;

whichDirection = 0; % 0 = both, 1 = l2r, -1 = r2l
flipVoxlr = 1; % flip voxels for plotting? 1 = yes

grid_mirror = [-grid(:,1),grid(:,2:3)];

% files in
eval(sprintf('files_in = fullfile(pwd, ''data_meg_source/%s/%s/%s/power'')',TFAmethod,whichFilter,whichGrid));

%% Loop subjects, trials

if whichDirection == 1
    subjects = subjects(mod(subjects,4)==1 | mod(subjects,4)==2);
elseif whichDirection == -1
    subjects = subjects(mod(subjects,4)==3 | mod(subjects,4)==0);
end
n_subjects = length(subjects);

subjNo = 0;

for i_subj = subjects
    subjNo = subjNo + 1;
    
    clear pow_ses1 pow_ses2
    
    for i_ses = sessions
        
        if i_subj == 18 && i_ses == 2
            ises = 1;
        else
            ises = i_ses;
        end
         
        for i_task = 1:3
            
            
            clear source_pow pow
            fprintf('Loading VP_%02d_%d_Cond_%d%s.mat...\n',i_subj,ises,i_task,whichTrials)
            eval(sprintf('source_path = fullfile(files_in, [''VP_%02d_%d_Cond_%d%s.mat'']);',i_subj,ises,i_task,whichTrials));
            source_pow = load(source_path);
            
            [n_grid, n_freq, n_tw] = size(source_pow.pow);
            
            pow = source_pow.pow;
            
            %% Baseline correction
            if blc
                clear pow_bl
                pow_bl = repmat(mean(pow(:,:,bsl_ind),3),[1 1 n_tw]);
                pow = 100*(pow-pow_bl)./pow_bl;
            end
            
            eval(sprintf('pow_ses%d = pow;',i_ses))
              
            %% write into matrix
            if i_subj == subjects(1) && i_ses == 1
                eval(sprintf('pow%d = zeros([n_subjects size(pow)]);',i_task))
            end
            
            if i_ses == 2
                eval(sprintf('pow%d(subjNo,:,:,:) = (pow_ses1 + pow_ses2)/2;',i_task))
            end
        end
    end
end

pow0 = (pow1+pow2+pow3)./3;

disp('Done!')

return;


%% Clusterstat for regular grids
clear stats cond1 cond2 vox

% pow1 = Vis. temp. prediction
% pow2 = Tact. temp. prediction
% pow3 = Luminance matching

test1 = pow1;
test2 = pow3;

toi = [source_pow.time_stimOff(1) source_pow.time_stimOff(end)]; % [-.3 0.9] time of interest in s / effects on sensor level in beta band
foi = [15 30]; % frequencies of interest in Hz

avgovertime = 1;
avgoverfreq = 1;

% Build fieldtrip-like cube
mni_cube = zeros(n_voxel_cube,3);
counter = 0;
for ix = xgrid:-1:1
    for iy = ygrid:-1:1
        for iz = zgrid:-1:1
            counter = counter + 1;
            mni_cube(counter,:) = [xrange(ix) yrange(iy) zrange(iz)];
        end
    end
end


cond1.dimord = 'pos_subj_freq_time';
cond1.dim = [xgrid,ygrid,zgrid];
cond1.pos = mni_cube*10;


for igrid = 1:n_grid
    dummy = round(grid(igrid,:)*100)/100;
    [~,idx] = ismember(dummy(1),xrange);
    [~,idy] = ismember(dummy(2),yrange);
    [~,idz] = ismember(dummy(3),zrange);
    vox(igrid) = sub2ind([xgrid,ygrid,zgrid],idx,idy,idz);
end
cond1.inside = false(n_voxel_cube,1);
cond1.inside(vox) = true;

freqidx = dsearchn(source_pow.freq',foi');
timeidx = dsearchn(source_pow.time',toi');

if avgovertime
    cond1.time = 1;
else
    cond1.time    = timeidx(1):timeidx(2);
end

if avgoverfreq
    cond1.freq = 1;
else
    cond1.freq    = freqidx(1):freqidx(2);
end

cond1.dimord = 'pos_subj_freq_time';
cond1.pow = zeros(n_voxel_cube,n_subjects,length(cond1.freq),length(cond1.time));
cond2 = cond1;

test1 = test1(:,:,freqidx(1):freqidx(2),timeidx(1):timeidx(2));
test2 = test2(:,:,freqidx(1):freqidx(2),timeidx(1):timeidx(2));

if avgovertime && avgoverfreq
    test1 = mean(mean(test1,4),3);
    test2 = mean(mean(test2,4),3);
elseif avgovertime
    test1 = mean(test1,4);
    test2 = mean(test2,4);
elseif avgoverfreq
    test1 = mean(test1,3);
    test2 = mean(test2,3);
end

if flipVoxlr
    subj_ind = 1:n_subjects;
    subj_rl = subj_ind(mod(subjects,4)==3 | mod(subjects,4)==0);
    for t = 1:length(cond1.time)
        for f = 1:length(cond1.freq)
            for j = subj_rl
                test1(j,:,f,t) = spatfiltergauss(test1(j,:,f,t)',grid_mirror,.1,grid);
                test2(j,:,f,t) = spatfiltergauss(test2(j,:,f,t)',grid_mirror,.1,grid);
            end
        end
    end
end

cond1.pow(vox,:,:,:) = permute(test1,[2,1,3,4]);
cond2.pow(vox,:,:,:) = permute(test2,[2,1,3,4]);

cfg                  = [];
cfg.dim              = [xgrid,ygrid,zgrid];
cfg.latency          = 'all';
cfg.parameter        = 'pow';
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'depsamplesT';
cfg.computeprob      = 'yes';
cfg.correctm         = 'cluster'; 
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0; 
cfg.tail             = 0; 
cfg.alpha            = 0.025/2;
cfg.numrandomization = 1000;

design = zeros(2,2*n_subjects);
design(1,:) = repmat(1:n_subjects,1,2);
design(2,:) = mod(floor([0:(2*n_subjects-1)]/(n_subjects/1)),2)+1;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stats] = ft_sourcestatistics(cfg, cond1, cond2);



%% Plot with hhtb
t = 1;
f = 1;

clear tstat
tstat_cube = stats.stat.*stats.mask;
tstat = tstat_cube(vox,f,t);

if sum(tstat)==0, error('Nothing to plot!'); end
tstat_smoothed = spatfiltergauss(tstat,grid,0.5,sa.cortex10K.vc);
para.myviewdir=[-1 -1 0];
figure; showsurface(sa.cortex10K,para,tstat_smoothed);
colormap(parula_transp)
set(gca,'clim',[-4 4])
para.colormaps={'parula'};
para.orientation='axial';
para.colorlimits = [-7 7];


%% Plot surface (left,right)
t = 1;
f = 1;

clear tstat
tstat_cube = stats.stat.*stats.mask;
tstat = tstat_cube(vox,f,t);
figure('position',[300 300 1000 800])

% plot 1
if sum(tstat)==0, error('Nothing to plot!'); end
tstat_smoothed = spatfiltergauss(tstat,grid,0.5,sa.cortex.vc);
para.myviewdir=[-1 0 0];
loc=[0 0 0]; dir=[-1 0 0];
ccut=cutsurface(sa.cortex,loc,dir);
subplot(2,2,1); showsurface(ccut,para,tstat_smoothed);
colormap(parula_transp)
set(gca,'clim',[-4 4])
colorbar('off')

% plot 2
para.myviewdir= [1 0 0];
dir=[1 0 0];
ccut=cutsurface(sa.cortex,loc,dir);
subplot(2,2,2); showsurface(ccut,para,tstat_smoothed);
colormap(parula_transp)
set(gca,'clim',[-4 4])
colorbar('off')
position = get(gca,'position');
position(1) = position(1) - 0.13;
position = set(gca,'position',position);

% plot 3
para.myviewdir= [1 0 0];
dir=[-1 0 0];
ccut=cutsurface(sa.cortex,loc,dir);
subplot(2,2,3); showsurface(ccut,para,tstat_smoothed);
colormap(parula_transp)
set(gca,'clim',[-4 4])
colorbar('off')
position = get(gca,'position');
position(2) = position(2) + 0.18;
position = set(gca,'position',position);

% plot 4
para.myviewdir= [-1 0 0];
dir=[1 0 0];
ccut=cutsurface(sa.cortex,loc,dir);
subplot(2,2,4); showsurface(ccut,para,tstat_smoothed);
colormap(parula_transp)
set(gca,'clim',[-4 4])
colorbar('off')
position = get(gca,'position');
position(1) = position(1) - 0.12;
position(2) = position(2) + 0.18;
position = set(gca,'position',position);

