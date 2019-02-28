%{
    file_name : preproLangPOSTTG.m
    author : Francesco Mantegna
    affiliation : Max Planck Institut for Psycholinguistics
    project : Music&Poetry
    date : 20/01/2018
%}

%% add fieldtrip path

full_path_to_fieldtrip = '/Users/francesco/Documents/MATLAB/fieldtrip-20170414';
addpath(full_path_to_fieldtrip)
ft_defaults

%% defining variables

twin      =    [-0.80 1.3];
rawDir    =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/RawData';
outDir    =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/Preprocessing/Data0.1-30POSTTGOSC';
subjnum   =    31; % subject 6 and 24 were excluded

%% Loop through the whole dataset    
    
for s = 1:subjnum
    
    %% moving to the raw data folder
    
    cd(rawDir)
    
    %% loading raw files

    cfg = [];
    if s <= 9 
        cfg.headerfile   = ['0' int2str(s) '_language.vhdr']; 
        cfg.datafile     = ['0' int2str(s) '_language.eeg'];
    else
        cfg.headerfile   = [int2str(s) '_language.vhdr']; 
        cfg.datafile     = [int2str(s) '_language.eeg'];
    end

    %% Detrend

    cfg.detrend = 'yes';    %   Remove linear trends
    % If you do not detrend, the power of the center frequency of the linear trend
    % (this frequency is very low), can bleed into other bins.
    % http://www.fieldtriptoolbox.org/faq/why_does_my_tfr_look_strange_part_ii
    
    %% Don't demean
    
    cfg.demean = 'no';
    % if you're just using the default in ft_freqanalysis, your data will be automatically demeaned 
    % (you don't have to do in separately in ft_preprocessing)
    % if you do not demean, the 0Hz bin can bleed into all other frequency
    % bins.
    % http://www.fieldtriptoolbox.org/faq/why_does_my_tfr_look_strange
     
    %% Fitering options TFR analysis
     
    cfg.continuous      = 'yes';
     
    % High-pass filters should be applied only to continuous data and not
    % to epoched data. This is because the edge artifact of a 0.5-Hz filter 
    % may last up to 6 s, which is probably longer than your epochs.
    cfg.lpfilter        = 'yes';
    cfg.lpfreq          = 30;
    cfg.hpfilter        = 'yes';    
    cfg.hpfreq          = 0.1; % is 1Hz too restrictive?
    cfg.hpfiltord       = 4; % defines how steep or shallow the transition between passband and stopband is. The higher the steeper.
    
    %% Re-referencing options - see explanation below
   
    cfg.reref         = 'yes';
    cfg.implicitref   = 'REF';
    cfg.refchannel    = {'RM' 'REF'};

    datacnt = ft_preprocessing(cfg);
   
    %% defining trials
    
    cfg.trialfun     = 'francesco_trialfun_TGlang';
    
    cfg.trialdef.pre  = abs(twin(1));
    cfg.trialdef.post = twin(2);
    cfg = ft_definetrial(cfg);
   if s >= 9 & s <= 10
       cfg.trl(45,:) = [];
   end
    
    data = ft_redefinetrial(cfg, datacnt);
    
    clear datacnt
    
    %% Eye electrodes
    
    % we remove the original eye electrodes, and we append
    % later on the new ones, collapsed 
    
    % vertical
    cfgeog              = [];
    cfgeog.channel      = {'LBEOG' 'C64'};
    cfgeog.reref        = 'yes';
    cfgeog.refchannel   = {'LBEOG'};
    EOGV             = ft_preprocessing(cfgeog, data);

    % only keep one channel, and rename to EOGV
    cfgeog              = [];
    cfgeog.channel      = 'C64';
    EOGV             = ft_selectdata(cfgeog, EOGV); 
    EOGV.label       = {'EOGV'};

    % horizontal
    cfgeog              = [];
    cfgeog.channel      = {'REOG' 'LEOG'};                                    
    cfgeog.reref        = 'yes';
    cfgeog.refchannel   = {'REOG'};                                             
    EOGH             = ft_preprocessing(cfgeog, data);

    % only keep one channel, and rename to EOGH
    cfgeog              = [];
    cfgeog.channel      = 'LEOG';                                               
    EOGH             = ft_selectdata(cfgeog, EOGH); 
    EOGH.label       = {'EOGH'};

    % only keep all non-EOG channels

    cfgeog         = [];
    cfgeog.channel = setdiff(1:60, 32);              % you can use either strings or numbers as selection
    data        = ft_selectdata(cfgeog, data); 

    % append the EOGH and EOGV channel to the 60 selected EEG channels 

    cfgeog = [];
    data = ft_appenddata(cfgeog, data, EOGH, EOGV);

    %% Neighborhood selection

    cfgneigh = [];
    cfgneigh.neighbourdist = 67;    
    cfgneigh.layout = 'mpi_customized_acticap64_fran';
    cfgneigh.method = 'distance';
    neighbours = ft_prepare_neighbours(cfgneigh);
    
    for i=1:length(neighbours)
        neighborhood(i,:) = size(neighbours(i).neighblabel,1);
    end

    mean(neighborhood)

    %% Finding GND,REF and RM and deleting them from the neighbourhood

    for i=1:length(neighbours)
        [pos,~,ground] = find(strcmp(neighbours(i).neighblabel, 'GND'));
        if ground  ~= 0
            neighbours(i).neighblabel(pos)= [];
        end
        [pos,~,leftref] = find(strcmp(neighbours(i).neighblabel, 'REF'));
        if leftref ~= 0  
            neighbours(i).neighblabel(pos)= [];
        end
        [pos,~,rightref] = find(strcmp(neighbours(i).neighblabel, 'RM'));
        if rightref ~= 0     
            neighbours(i).neighblabel(pos)= [];
        end
    end

    clear ground leftref rightref

    %% Channel Visual inspection
    
    % look for noisy or flat channels and decide whether to interpolate them or not
    cfgch = [];       
    cfgch.method = 'channel';
    cfgch.channel  = setdiff(1:59, [22 26 53]);
    cfgch.keepchannel = 'yes';
    data = ft_rejectvisual(cfgch, data);

    %% Topographic Interpolation

    load ('fran_59CH_elec')
    interpolation = true;
    prompt1       = {'Enter the number of the channels that you want to interpolate'};
    prompt2       = {'Do you want to go on with this channel? (1=yes, 0=no)'};
    name          = 'Topographic Interpolation';
    while interpolation
        for i = 1:5
            candidate = inputdlg(prompt1,name);
            if  isnan(str2double(candidate)) | i == 6
                interpolation = false;
                break
            end
            disp('###Have a look at the layout neighborhood first###')
            neighborhood(str2double(candidate))
            answer = inputdlg(prompt2,name);
            if  answer{1} == '0'
                break
            elseif answer{1} == '1'
                if ismember(str2double(candidate),[22 26 54])
                    disp('###Do not interpolate occipital channels!###')
                    disp('###Press a key to go on with the next channel###')
                    pause;
                else
                    ch = str2double(candidate);
                    cfginter = [];
                    cfginter.badchannel = ['C' int2str(ch)]; % channel name
                    cfginter.neighbours = neighbours;
                    cfginter.method = 'weighted';
                    cfginter.elec = elec;
                    data = ft_channelrepair(cfginter, data);
                end
            end
        end
    end
    
    clear cfgeog cfgch cfginter cfgneigh interpolation name ch prompt answer

    %% Identifying data segments exceeding voltage step threshold 
    
    cfgz                                 = [];
    cfgz.trl                             = cfg.trl;
    cfgz.headerfile                      = cfg.headerfile;
    cfgz.datafile                        = cfg.datafile;
    cfgz.continuous                      = 'yes';
    cfgz.artfctdef.threshold.channel     = setdiff(1:59, [22 26 53]);
    cfgz.artfctdef.threshold.bpfilter    = 'no';
    cfgz.artfctdef.threshold.derivative  = 'yes';
    cfgz.artfctdef.threshold.max         = 50;
   
    [~, artifact_THR1] = ft_artifact_threshold(cfgz,data);  
    
    %% Rejecting bad trials completely
   
    cfgz=[]; 
    cfgz.artfctdef.reject = 'complete'; % this rejects complete trials
    cfgz.artfctdef.crittoilim = [-0.5 1];
    cfgz.artfctdef.threshold.artifact = artifact_THR1;
    data_no_voltagestep = ft_rejectartifact(cfgz,data);
    
    %% Identifying data segments exceeding absolute difference threshold 
    
    cfgz                                 = [];
    cfgz.trl                             = data_no_voltagestep.cfg.trl;
    cfgz.headerfile                      = cfg.headerfile;
    cfgz.datafile                        = cfg.datafile;
    cfgz.continuous                      = 'yes'; 
    cfgz.artfctdef.threshold.channel     = setdiff(1:59, [22 26 53]);
    cfgz.artfctdef.threshold.bpfilter    = 'no';
    cfgz.artfctdef.threshold.absdiff     = 'yes';
    cfgz.artfctdef.threshold.max         = 100;
   
    [~, artifact_THR2] = ft_artifact_threshold(cfgz, data_no_voltagestep);  
    
    %% Rejecting bad trials completely
   
    cfgz=[]; 
    cfgz.artfctdef.reject = 'complete'; % this rejects complete trials
    cfgz.artfctdef.crittoilim = [-0.5 1];
    cfgz.artfctdef.threshold.artifact = artifact_THR2;
    data_no_absdiff = ft_rejectartifact(cfgz,data_no_voltagestep);
    
    %% Identifying data segments exceeding absolute difference threshold 
    
    cfgz                                 = [];
    cfgz.trl                             = data_no_absdiff.cfg.trl;
    cfgz.headerfile                      = cfg.headerfile;
    cfgz.datafile                        = cfg.datafile;
    cfgz.continuous                      = 'yes'; 
    cfgz.artfctdef.threshold.channel     = setdiff(1:59, [22 26 53]);
    cfgz.artfctdef.threshold.bpfilter    = 'no';
    cfgz.artfctdef.threshold.range        = 200; % +-100/muV
   
    [~, artifact_THR3] = ft_artifact_threshold(cfgz, data_no_absdiff);
    
    %% Rejecting bad trials completely
   
    cfgz=[]; 
    cfgz.artfctdef.reject = 'complete'; % this rejects complete trials
    cfgz.artfctdef.threshold.artifact = artifact_THR3;
    cfgz.artfctdef.crittoilim = [-0.5 1];
    data_no_ampl = ft_rejectartifact(cfgz,data_no_absdiff);
    
    %% Identifying data segments exceeding 4 z-value cutoff
    
    % EOG
    cfgz.trl        = data_no_ampl.cfg.trl;
    cfgz.headerfile = cfg.headerfile;
    cfgz.datafile   = cfg.datafile;
    cfgz.continuous = 'yes'; 
 
    % channel selection, cutoff and padding
    cfgz.artfctdef.zvalue.channel     = 61;
    cfgz.artfctdef.zvalue.cutoff      = 4.5;
    cfgz.artfctdef.zvalue.trlpadding  = 0;
    cfgz.artfctdef.zvalue.artpadding  = 0.1;
    cfgz.artfctdef.zvalue.fltpadding  = 0;
  
    % algorithmic parameters
    cfgz.artfctdef.zvalue.bpfilter   =  'yes';
    cfgz.artfctdef.zvalue.bpfilttype =  'but';
    cfgz.artfctdef.zvalue.bpfreq     =  [1 15];
    cfgz.artfctdef.zvalue.bpfiltord  =  4;
    cfgz.artfctdef.zvalue.hilbert    =  'yes';
 
    % feedback
    cfgz.artfctdef.zvalue.interactive = 'yes';
 
    [~, artifact_EOG] = ft_artifact_zvalue(cfgz,data_no_ampl); 
   
    %% Rejecting bad trials completely
   
    cfgz=[]; 
    cfgz.artfctdef.reject = 'complete'; % this rejects complete trials
    cfgz.artfctdef.eog.artifact = artifact_EOG;
    cfgz.artfctdef.crittoilim = [-0.5 1];
    data_no_artifacts = ft_rejectartifact(cfgz,data_no_ampl);
    
    clear cfgz data_no_voltagestep artifact_THR1 artifact_THR2 artifact_THR3 artifact_EOG data_no_ampl data_no_absdiff
    
    %% Reject trials based on variance 
    
    cfg=[];
    cfg.metric = 'zvalue';  % use by default zvalue method
    cfg.method = 'summary';
    data_no_artifacts = ft_rejectvisual(cfg, data_no_artifacts); % choose a couple of methods to make sure tehre are no outliers left

    %% Visual/manual artifact rejection per trial

    cfg = [];
    cfg.method = 'trial';
    data_no_artifacts = ft_rejectvisual(cfg, data_no_artifacts);

    %% Check number of trials per condition

    conditions = unique(data_no_artifacts.trialinfo);
    trlXcond = [conditions,histc(data_no_artifacts.trialinfo(:),conditions)];
    fprintf('the overall number of trials after pre-processing is the following:\n');
    fprintf('%g %g %g\n',trlXcond);
    
   %% exclude irrelevant and noisy electrodes from the analyses
    
    % exclude peripherical occipital electrodes C22 C26 C54 and
    % ocular electrodes EOGV EOGH
    cfg                      = []; 
    cfg.channel              = setdiff(1:60, [22 26 53]);       
    data_no_artifacts        = ft_selectdata(cfg, data_no_artifacts);

    %% Time-frequency analysis post-stimulus for each condition Low Frequencies (2-30 Hz)
    
    for icond = 1:3
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = 'C*';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.trials       = find(data_no_artifacts.trialinfo==icond); 
        cfg.foi          = 2:1:30;                         % analysis 2 to 30 Hz in steps of 1 Hz
        cfg.pad = ceil(max(cellfun(@numel, data.time)/data.fsample));
        % This changes the frequencies you can estimate from your data 
        % by artificially creating longer T's (and thus changing the Rayleigh frequency). 
        % Importantly though, this does not (!) change the intrinsic frequency resolution of your data.
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.4;   % length of time window = 0.5 sec
        %cfg.t_ftimwin    = 3./cfg.foi; % it is recommended to not use less than 3 cycles
        cfg.toi          = -0.5:0.01:1;   % time window "slides" from -0.5 to 1.5 sec in steps of 0.01 sec (10 ms)
        if icond == 1
            POSTlowTFRhann1 = ft_freqanalysis(cfg,data_no_artifacts);
        elseif icond == 2
            POSTlowTFRhann2 = ft_freqanalysis(cfg,data_no_artifacts);
        elseif icond == 3
            POSTlowTFRhann3 = ft_freqanalysis(cfg,data_no_artifacts);
        end
    end

%     %% Time-frequency analysis post-stimulus for each condition High Frequencies (25-80 Hz)
%     
%     for icond = 1:3
%         cfg              = [];
%         cfg.output       = 'pow';
%         cfg.channel      = 'C*';
%         cfg.method       = 'mtmconvol';
%         cfg.trials       = find(data_no_artifacts.trialinfo==icond);    
%         cfg.foi          = 25:2.5:80;                         % analysis 30 to 80 Hz in steps of 1 Hz 
%         cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.4;   % length of time smoothing window = 0.5 sec
%         % when people use multitapers to analyze gamma they usually 
%         % use fixed rather than scaling time windows and smoothing.
%         cfg.pad = ceil(max(cellfun(@numel, data.time)/data.fsample));
%         % This changes the frequencies you can estimate from your data 
%         % by artificially creating longer T's (and thus changing the Rayleigh frequency). 
%         % Importantly though, this does not (!) change the intrinsic frequency resolution of your data.
%         %cfg.tapsmofrq    = 5;   % fixed 5 Hz frequency smoothing window 
%         cfg.tapsmofrq    = 0.2*cfg.foi; % frequency smoothing window increase with frequency
%         cfg.toi          = -0.5:0.01:1;   % time window "slides" from -0.5 to 1 sec in steps of 0.01 sec (10 ms)
%         if icond == 1
%             POSThighTFRmult1 = ft_freqanalysis(cfg,data_no_artifacts);
%         elseif icond == 2
%             POSThighTFRmult2 = ft_freqanalysis(cfg,data_no_artifacts);
%         elseif icond == 3
%             POSThighTFRmult3 = ft_freqanalysis(cfg,data_no_artifacts);
%         end
%     end
%     
%     %% Multiplot TFR Low Frequencies (2-30 Hz) per condition
%     
    condname1 = 'CONGRUENT';
    condname2 = 'INTERMEDIATE';
    condname3 = 'INCONGRUENT';
%     
%     for icond = 1:3
%         figure;
%         cfg = [];
%         cfg.baseline     = [-0.5 -0.15]; 
%         cfg.baselinetype = 'relchange';
%         cfg.colormap     = jet(148);
%         cfg.zlim         = 'maxmin'; %[-10e-1 10e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
%         cfg.fontsize     = 12;
%         cfg.showlabels   = 'yes';
%         cfg.interactive  = 'yes';
%         cfg.comment      = ' ';
%         cfg.layout       = 'mpi_customized_acticap64_fran.mat';
%         ft_multiplotTFR(cfg, eval(['POSTlowTFRhann', num2str(icond)]));
%         title(['Subj', num2str(s),' TFR LF multiplot ', eval(['condname',num2str(icond)])]);
%     end
%     
    %% Singleplot on C30 TFR Low Frequencies (2-30 Hz) per condition
    
    figure;
    for icond = 1:3   
        cfg = [];
        cfg.baseline     = [-0.50 -0.15];
        cfg.baselinetype = 'relchange';  
        cfg.maskstyle    = 'opacity';
        cfg.colormap     = jet(148);
        cfg.fontsize     = 12;
        cfg.zlim         = 'maxmin'; %[-10e-1 10e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
        cfg.channel      = 'C35';
        cfg.colorbar       = 'no';
        if icond == 1
            h1 = subplot(1,3,1); ft_singleplotTFR(cfg, POSTlowTFRhann1);
            originalSize1 = get(gca, 'Position');
            title({['Subj', num2str(s), ' TFR LF C35'],['\color{black}', condname1]});
        elseif icond == 2
            h2 = subplot(1,3,2); ft_singleplotTFR(cfg, POSTlowTFRhann2);
            originalSize2 = get(gca, 'Position');
            title({['Subj', num2str(s), ' TFR LF C35'],['\color{blue}', condname2]});
        elseif icond == 3
            cfg.colorbar       = 'no';
            h3 = subplot(1,3,3); ft_singleplotTFR(cfg, POSTlowTFRhann3);
            originalSize3 = get(gca, 'Position');
            c = colorbar;
            c.Label.String = '\muV';
            c.Label.FontSize = 12;
            set(c,'YTick', -1:1:1);
            title({['Subj', num2str(s), ' TFR LF C35'],['\color{red}', condname3]});
        end
        hold on;
        yL = get(gca, 'YLim');
        plot([0 0], yL, 'k-');
        set(eval(['h',num2str(icond)]),'YTick', 5:5:30);
        set(eval(['h',num2str(icond)]), 'Position', eval(['originalSize',num2str(icond)]));
    end
   
%     %% Multiplot TFR High Frequencies (2-30 Hz) per condition
%     
%    for icond = 1:3
%        figure;
%        cfg = [];
%        cfg.baseline     = [-0.5 -0.15]; 
%        cfg.baselinetype = 'relchange';
%        cfg.colormap     = jet(148);
%        cfg.zlim         = 'maxmin'; %[-6e-1 6e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
%        cfg.fontsize     = 12;
%        cfg.showlabels   = 'yes';
%        cfg.interactive  = 'yes';
%        cfg.comment      = ' ';
%        cfg.layout       = 'mpi_customized_acticap64_fran.mat';
%        ft_multiplotTFR(cfg, eval(['POSThighTFRmult', num2str(icond)]));
%        title(['Subj', num2str(s),' TFR multiplot ', eval(['condname',num2str(icond)])]);
%    end
%     
%     %% Singleplot on C30 TFR High Frequencies (25-80 Hz) per condition
%      
%     figure;
%     for icond = 1:3   
%        cfg = [];
%        cfg.baseline     = [-0.50 -0.15];
%        cfg.baselinetype = 'relchange';  
%        cfg.maskstyle    = 'opacity';
%        cfg.colormap     = jet(148);
%        cfg.fontsize     = 12;
%        cfg.ylim         = [30 80];
%        cfg.zlim         = 'maxmin'; %[-9e-1 9e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
%        cfg.channel      = 'C35';
%        cfg.colorbar       = 'no';
%        if icond == 1
%            h4 = subplot(1,3,1); ft_singleplotTFR(cfg, POSThighTFRmult1);
%            originalSize1 = get(gca, 'Position');
%            title({['Subj', num2str(s), ' TFR LF C35'],['\color{black}', condname1]});
%        elseif icond == 2
%            h5 = subplot(1,3,2); ft_singleplotTFR(cfg, POSThighTFRmult2);
%            originalSize2 = get(gca, 'Position');
%            title({['Subj', num2str(s), ' TFR LF C35'],['\color{blue}', condname2]});
%        elseif icond == 3
%            h6 = subplot(1,3,3); ft_singleplotTFR(cfg, POSThighTFRmult3);
%            originalSize3 = get(gca, 'Position');
%            c = colorbar;
%            c.Label.String = '\muV';
%            c.Label.FontSize = 12;
%            set(c,'YTick', -1:1:1);
%            title({['Subj', num2str(s), ' TFR LF C35'],['\color{red}', condname3]});
%        end
%        hold on;
%        yL = get(gca, 'YLim');
%        plot([0 0], yL, 'k-');
%        set(eval(['h',num2str(icond+3)]),'YTick', 30:10:80);
%        set(eval(['h',num2str(icond+3)]), 'Position', eval(['originalSize',num2str(icond)]));
%    end
% 
%    clear condname1 condname2 condname3
%    
    %% Save data
    
   % moving to the output folder
   
   cd(outDir)
%    if s<=6
       outputfile = ['preprolangPOSTTG_SUBJ' int2str(s) '.mat'];
%    elseif s>=7 & s<=23
%        outputfile = ['preprolangPOSTTG_SUBJ' int2str(s-1) '.mat'];
%    elseif s>=23
%        outputfile = ['preprolangPOSTTG_SUBJ' int2str(s-1) '.mat'];
%    end
   save(outputfile,'data_no_artifacts','trlXcond','POSTlowTFRhann1', ...
       'POSTlowTFRhann2', 'POSTlowTFRhann3')%'POSThighTFRmult1', ...
       %'POSThighTFRmult2', 'POSThighTFRmult3')
   
   disp('###Press a key to go on with the next subject###')  % Press a key here
   pause;
   clf
   clearvars -except rawDir outDir subjnum twin s
   close all
end

%% Removing fieldtrip path

rmpath(genpath('/Users/francesco/Documents/MATLAB/fieldtrip-20170414'))
