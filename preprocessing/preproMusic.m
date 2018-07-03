%{
    file_name : preproMusic.m
    author : Francesco Mantegna
    institution : Max Planck Institut for Psycholinguistics
    project : Music&Poetry
    date : 10/01/2018
%}

%% add fieldtrip path

full_path_to_fieldtrip = '/Users/francesco/Documents/MATLAB/fieldtrip-20170414';
addpath(full_path_to_fieldtrip)
ft_defaults

%% defining variables

twin      =    [-0.5 1];
rawDir    =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/RawData';
outDir    =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/Preprocessing/Data0.1-30MUSIC';
subjnum   =    31;

for s = 1:subjnum
    
    %% moving to the raw data folder
    
    cd(rawDir)
    
    %% creating a configuration

    cfg = [];
    cfg.trialfun     = 'francesco_trialfun_ERPmusic';
    if s <= 9 
        cfg.headerfile   = ['0' int2str(s) '_music.vhdr']; 
        cfg.datafile     = ['0' int2str(s) '_music.eeg'];
    else
        cfg.headerfile   = [int2str(s) '_music.vhdr']; 
        cfg.datafile     = [int2str(s) '_music.eeg'];
    end

    %% defining trials

    cfg.trialdef.pre  = abs(twin(1));
    cfg.trialdef.post = twin(2);
    cfg = ft_definetrial(cfg);
    if s == 1
        cfg.trl(47:48,:) = [];
    elseif s == 2
        cfg.trl(96,:) = [];
    elseif s == 9
        cfg.trl(96,:) = [];
    elseif s == 16
        cfg.trl(48,:) = [];
    elseif s == 24
        cfg.trl(48,:) = [];
    elseif s == 26
        cfg.trl(48,:) = [];
    elseif s == 31
        cfg.trl(96,:) = [];
    end
        
    %% Baseline-correction options

%     cfg.demean          = 'yes';
%     cfg.baselinewindow  = [-0.2 0];

    %% Detrend

%     cfg.detrend = 'yes';    %   Remove linear trends form data

    %% Fitering options

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
    
    data = ft_preprocessing(cfg);
    
    % for subject 2 the electrode sets were inverted during the recording
    % of the music part
    if s ==2
        data.label{1} = 'C33';
        data.label{2} = 'C34';
        data.label{3} = 'C35';
        data.label{4} = 'C36';
        data.label{5} = 'C37';
        data.label{6} = 'C38';
        data.label{7} = 'C39';
        data.label{8} = 'C40';
        data.label{9} = 'C41';
        data.label{10} = 'C42';
        data.label{11} = 'C43';
        data.label{12} = 'C44';
        data.label{13} = 'C45';
        data.label{14} = 'C46';
        data.label{15} = 'C47';
        data.label{16} = 'C48';
        data.label{17} = 'C49';
        data.label{18} = 'C50';
        data.label{19} = 'C51';
        data.label{20} = 'C52';
        data.label{21} = 'C53';
        data.label{22} = 'C54';
        data.label{23} = 'C55';
        data.label{24} = 'C56';
        data.label{25} = 'C57';
        data.label{26} = 'C58';
        data.label{27} = 'C59';
        data.label{28} = 'C60';
        data.label{29} = 'LEOG';
        data.label{30} = 'LBEOG';
        data.label{31} = 'REOG';
        data.label{32} = 'C64';
        data.label{33} = 'C1';
        data.label{34} = 'C2';
        data.label{35} = 'C3';
        data.label{36} = 'C4';
        data.label{37} = 'C5';
        data.label{38} = 'C6';
        data.label{39} = 'C7';
        data.label{40} = 'C8';
        data.label{41} = 'C9';
        data.label{42} = 'C10';
        data.label{43} = 'C11';
        data.label{44} = 'C12';
        data.label{45} = 'C13';
        data.label{46} = 'C14';
        data.label{47} = 'C15';
        data.label{48} = 'C16';
        data.label{49} = 'C17';
        data.label{50} = 'C18';
        data.label{51} = 'C19';
        data.label{52} = 'C20';
        data.label{53} = 'C21';
        data.label{54} = 'C22';
        data.label{55} = 'C23';
        data.label{56} = 'C24';
        data.label{57} = 'C25';
        data.label{58} = 'C26';
        data.label{59} = 'C27';
        data.label{60} = 'C28';
        data.label{61} = 'C29';
        data.label{62} = 'C30';
        data.label{63} = 'C31';
        data.label{64} = 'RM';
    end    

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
    if s == 2
        cfgeog.channel = setdiff(1:64, [29 30 31 32 64]);
    else
        cfgeog.channel = setdiff(1:60, 32);              % you can use either strings or numbers as selection
    end
    data        = ft_selectdata(cfgeog, data); 

    % append the EOGH and EOGV channel to the 60 selected EEG channels 

    cfgeog = [];
    data = ft_appenddata(cfgeog, data, EOGH, EOGV);


    %% Channel layout
%     cfg = [];
%     cfg.layout   = 'mpi_customized_acticap64_fran.mat';
%     ft_layoutplot(cfg);

    %% Visual inspection

    cfgch = [];       % see all the electrodes for all trials, check if it's messy or flat and decide whether to interpolate or not
    cfgch.method = 'channel';
    data = ft_rejectvisual(cfgch, data);

    %% Neighborhood selection
    if s == 2
        
    else
        cfgneigh = [];
        cfgneigh.neighbourdist = 67;    
        cfgneigh.layout = 'mpi_customized_acticap64_fran';
        %cfgneigh.method = 'distance';
        cfgneigh.method = 'triangulation';
        neighbours = ft_prepare_neighbours(cfgneigh);

        %% find GND,REF and RM to delete them from the neighbourhood

        for i=1:length(neighbours)
            [pos,clm,ground] = find(strcmp(neighbours(i).neighblabel, 'GND'));
            if ground  ~= 0
                neighbours(i).neighblabel(pos)= [];
            end
            [pos,clm,leftref] = find(strcmp(neighbours(i).neighblabel, 'REF'));
            if leftref ~= 0  
                neighbours(i).neighblabel(pos)= [];
            end
            [pos,clm,rightref] = find(strcmp(neighbours(i).neighblabel, 'RM'));
            if rightref ~= 0     
                neighbours(i).neighblabel(pos)= [];
            end
        end

        clear ground
        clear leftref
        clear rightref

        %% Interpolation

        load ('fran_59CH_elec')
        interpolation = true;
        prompt={'Enter the number of the channels that you want to interpolate'};
        name = 'Topographic Interpolation';
    while interpolation
        for i = 1:5
            answer = inputdlg(prompt,name);
            if  isnan(str2double(answer)) || i == 6
                interpolation = false;
                break
            end
            if ismember(str2double(answer),[22 26 54])
                disp('###Do not interpolate occipital channels!###')
                disp('###Press a key to go on with the next channel###')
                pause;
            else
                ch = str2double(answer);
                cfginter = [];
                cfginter.badchannel = ['C' int2str(ch)]; % channel name
                cfginter.neighbours = neighbours;
                cfginter.method = 'weighted';
                cfginter.elec = elec;
                data = ft_channelrepair(cfginter, data);
            end
        end
    end
    
    clear cfgeog cfgch cfginter cfgneigh interpolation name ch prompt answer
    end

    %% Identifying data segments exceeding voltage step threshold 
    
    cfgz                                 = [];
    cfgz.trl                             = cfg.trl;
    cfgz.headerfile                      = cfg.headerfile;
    cfgz.datafile                        = cfg.datafile;
    cfgz.continuous                      = 'yes';
    if s == 2
        cfgz.artfctdef.threshold.channel     = setdiff(1:59, [50 54 22]);
    else
        cfgz.artfctdef.threshold.channel     = setdiff(1:59, [22 26 53]);
    end
    cfgz.artfctdef.threshold.bpfilter    = 'no';
    cfgz.artfctdef.threshold.derivative  = 'yes';
    cfgz.artfctdef.threshold.max         = 50;
   
    [~, artifact_THR1] = ft_artifact_threshold(cfgz,data);  
    
    %% Rejecting bad trials completely
   
    cfgz=[]; 
    cfgz.artfctdef.reject = 'complete'; % this rejects complete trials
    cfgz.artfctdef.threshold.artifact = artifact_THR1;
    cfgz.artfctdef.crittoilim = [-0.2 1];
    data_no_voltagestep = ft_rejectartifact(cfgz,data);
    
    %% Identifying data segments exceeding absolute difference threshold 
    
    cfgz                                 = [];
    cfgz.trl                             = data_no_voltagestep.cfg.trl;
    cfgz.headerfile                      = cfg.headerfile;
    cfgz.datafile                        = cfg.datafile;
    cfgz.continuous                      = 'yes'; 
    if s == 2
        cfgz.artfctdef.threshold.channel     = setdiff(1:59, [50 54 22]);
    else
        cfgz.artfctdef.threshold.channel     = setdiff(1:59, [22 26 53]);
    end
    cfgz.artfctdef.threshold.bpfilter    = 'no';
    cfgz.artfctdef.threshold.absdiff     = 'yes';
    cfgz.artfctdef.threshold.max         = 100;
   
    [~, artifact_THR2] = ft_artifact_threshold(cfgz,data_no_voltagestep);  
    
    %% Rejecting bad trials completely
   
    cfgz=[]; 
    cfgz.artfctdef.reject = 'complete'; % this rejects complete trials
    cfgz.artfctdef.threshold.artifact = artifact_THR2;
    cfgz.artfctdef.crittoilim = [-0.2 1];
    data_no_absdiff = ft_rejectartifact(cfgz,data_no_voltagestep);
    
    %% Identifying data segments exceeding absolute difference threshold 
    
    cfgz                                 = [];
    cfgz.trl                             = data_no_absdiff.cfg.trl;
    cfgz.headerfile                      = cfg.headerfile;
    cfgz.datafile                        = cfg.datafile;
    cfgz.continuous                      = 'yes'; 
    if s == 2
        cfgz.artfctdef.threshold.channel     = setdiff(1:59, [50 54 22]);
    else
        cfgz.artfctdef.threshold.channel     = setdiff(1:59, [22 26 53]);
    end
    cfgz.artfctdef.threshold.bpfilter    = 'no';
    cfgz.artfctdef.threshold.range        = 200; % +-150/muV
   
    [~, artifact_THR3] = ft_artifact_threshold(cfgz,data_no_absdiff);
    
    %% Rejecting bad trials completely
   
    cfgz=[]; 
    cfgz.artfctdef.reject = 'complete'; % this rejects complete trials
    cfgz.artfctdef.threshold.artifact = artifact_THR3;
    cfgz.artfctdef.crittoilim = [-0.2 1];
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
    cfgz.artfctdef.crittoilim = [-0.2 1];
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

    %% Timelockanalysis to compute the ERPs

%     cfg = [];
%     cfg.trials = find(data_no_artifacts.trialinfo==1);
%     avgcond1 = ft_timelockanalysis(cfg, data_no_artifacts);
% 
%     cfg = [];
%     cfg.trials = find(data_no_artifacts.trialinfo==2);
%     avgcond2 = ft_timelockanalysis(cfg, data_no_artifacts);
% 
%     cfg = [];
%     cfg.trials = find(data_no_artifacts.trialinfo==3);
%     avgcond3 = ft_timelockanalysis(cfg, data_no_artifacts);
    
%     %% customize figure in a publishable format
%     figure
%     
%     % first axis
%     ax1 = gca; % current axe
%     ax1.XAxisLocation = 'origin';
%     ax1_pos = ax1.Position; % position of first axes
% 
%     % ROI rectangles and zero line
% %     r1 = rectangle('Position',[0.10 -14 0.09 27.98],'FaceColor',[0.8 0.8 0.8]);
% %     set(r1,'EdgeColor','none');
% %     hold on;
% %     r2 = rectangle('Position',[0.20 -14 0.09 27.98],'FaceColor',[0.8 0.8 0.8]);
% %     set(r2,'EdgeColor','none');
% %     hold on;
% %     r3 = rectangle('Position',[0.30 -14 0.20 27.98],'FaceColor',[0.8 0.8 0.8]);
% %     set(r3,'EdgeColor','none');
% %     hold on;
%     yL = get(gca, 'YLim');
%     plot([0 0], yL, '--');
%     hold on;
%     set(gca, 'Layer', 'top');
%     
%     % figure settings
%     title(['musician no1 average across trials for one channel']);
%     %title(['Subj', num2str(s), ' average across trials for the 30th channel']);
%     box off
%     xlim([-0.2 1]);
%     ylim([-14 14]);
%     set(gca,'FontSize',16);
%     set(gca,'TickLength',[0.015 0.3]);
%     set(gca,'Tickdir', 'both');
%     set(gca,'Ydir','reverse');
%     set(gca,'Ytick',(-12:4:12));
%     set(gca,'Xticklabel',[]);
%     ylabel('Evoked Response [\muV]'); % y-axis label
%     
%     % plotting ERPs per condition
%     fig1 = plot(avgcond1.time, avgcond1.avg(10,:),'k', 'LineWidth', 2);
%     hold on;
%     fig2 = plot(avgcond2.time, avgcond2.avg(10,:),'b', 'LineWidth', 2);
%     hold on;
%     fig3 = plot(avgcond3.time, avgcond3.avg(10,:),'r', 'LineWidth', 2);
%     % costumize legend
%     legend([fig1, fig2, fig3], {'congruent','intermediate','incongruent'}, 'Location', 'northwest');
% 
%     % second axis
%     ax2 = axes('Position',[0.1380 0.1100 0.7750 0.8150],'XAxisLocation','bottom','YAxisLocation','left','Color','none');
%     axis([-0.2 1 -6 4]);
%     ax2.XRuler.Axle.LineStyle = 'none';
%     ax2.YRuler.Axle.LineStyle = 'none';
%     set(gca,'FontSize',16);
%     set(gca,'Xticklabel',(-0.2:0.2:1))
%     set(gca,'TickLength',[0 0.5]);
%     set(gca,'YTick',[]);
%     set(gca, 'Layer', 'bottom');
%     xlabel('Time [s]');
%    
%     %% Multiplot 
%     figure;
%     cfg = [];
%     cfg.layout = 'mpi_customized_acticap64_fran.mat';
%     cfg.graphcolor  = 'kbr';
%     cfg.interactive = 'yes';
%     cfg.showoutline = 'yes';
%     ft_multiplotER(cfg, avgcond1, avgcond2, avgcond3)
%     title(['Subj', num2str(s), ' multiplot averaged across trials']);

    %% Save data
    
    % moving to the output folder
    
    cd(outDir)

    outputfile = ['prepromusic_SUBJ' int2str(s) '.mat']; 
    save(outputfile,'data_no_artifacts','trlXcond'); %,'avgcond1','avgcond2','avgcond3')
    disp('###Press a key to go on with the next subject###')  % Press a key here
    pause;
    clf
    clearvars -except rawDir outDir subjnum twin s
    close all
end

rmpath(genpath('/Users/francesco/Documents/MATLAB/fieldtrip-20170414'))