%{
    file_name : TFALangPOSTTG.m
    author : Francesco Mantegna
    affiliation : Max Planck Institut for Psycholinguistics
    project : Music&Poetry
    date : 27/01/2018
%}

%% add fieldtrip path

full_path_to_fieldtrip2017 = '/Users/francesco/Documents/MATLAB/fieldtrip-20170414';
addpath(full_path_to_fieldtrip2017)
ft_defaults

% defining variables

twin      =    [-0.80 1.30];
inputDir  =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/Preprocessing/Data0.1-30POSTTGOSC';
outDir    =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/OSCanalysis/ResultsTFA0.1Hzcutoff';
subjnum   =    29; % subject 6 and 24 were excluded

% loading pre-processed files

cd(inputDir)

for s = 1:subjnum
    
    load(['preprolangPOSTTG_SUBJ' int2str(s) '.mat'])
    
%     clearvars -except twin inputDir outDir subjnum s EEGdataset data_no_artifacts trlXcond
    
%     exclude peripherical occipital electrodes C22 C26 C54
%     cfg                  = []; 
%     cfg.channel          = setdiff(1:59, [22 26 53]);
%     data_artifacts_free  = ft_selectdata(cfg, data_no_artifacts);
      data_artifacts_free  = data_no_artifacts;
    
    %% Time-frequency analysis post-stimulus for each condition Low Frequencies (2-30 Hz)
    
    for icond = 4:4
        cfg              = [];
        cfg.output       = 'pow';
%         cfg.output       = 'fourier';
        cfg.channel      = 'C*';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        if icond == 4
            cfg.trials   = 'all';
        elseif icond == 5
            cfg.trials   = find(data_artifacts_free.trialinfo==1);
        else
            cfg.trials   = find(data_artifacts_free.trialinfo==icond); 
        end
        cfg.foi          = 2:1:30;                         % analysis 2 to 30 Hz in steps of 1 Hz
        cfg.pad          = ceil(max(cellfun(@numel, data_artifacts_free.time)/data_artifacts_free.fsample));
        % This changes the frequencies you can estimate from your data 
        % by artificially creating longer T's (and thus changing the Rayleigh frequency). 
        % Importantly though, this does not (!) change the intrinsic frequency resolution of your data.
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.4;   % length of time window = 0.5 sec
        %cfg.t_ftimwin    = 3./cfg.foi; % it is recommended to not use less than 3 cycles
        cfg.toi          = -0.5:0.01:1;   % time window "slides" from -0.5 to 1.5 sec in steps of 0.01 sec (10 ms)
        if icond == 1
            POSTlowTFRhann1 = ft_freqanalysis(cfg,data_artifacts_free);
        elseif icond == 2
            POSTlowTFRhann2 = ft_freqanalysis(cfg,data_artifacts_free);
        elseif icond == 3
            POSTlowTFRhann3 = ft_freqanalysis(cfg,data_artifacts_free);
        elseif icond == 4
            POSTlowTFRhannALL = ft_freqanalysis(cfg,data_artifacts_free);
        elseif icond == 5
            POSTlowFOURIER1 = ft_freqanalysis(cfg,data_artifacts_free);
        end
    end
% 
%     %% Time-frequency analysis post-stimulus for each condition High Frequencies (25-80 Hz)
%     
%     for icond = 4:4
%         cfg              = [];
%         cfg.output       = 'pow';
%         cfg.channel      = 'C*';
%         cfg.method       = 'mtmconvol';
%         if icond == 4
%             cfg.trials   = 'all';
%         else
%             cfg.trials   = find(data_artifacts_free.trialinfo==icond); 
%         end    
%         cfg.foi          = 25:2.5:80;                         % analysis 30 to 80 Hz in steps of 2.5 Hz 
%         cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.4;   % length of time smoothing window = 0.5 sec
%         % when people use multitapers to analyze gamma they usually 
%         % use fixed rather than scaling time windows and smoothing.
%         cfg.pad          = ceil(max(cellfun(@numel, data_artifacts_free.time)/data_artifacts_free.fsample));
%         % This changes the frequencies you can estimate from your data 
%         % by artificially creating longer T's (and thus changing the Rayleigh frequency). 
%         % Importantly though, this does not (!) change the intrinsic frequency resolution of your data.
%         %cfg.tapsmofrq    = 5;   % fixed 5 Hz frequency smoothing window 
%         cfg.tapsmofrq    = 0.2*cfg.foi; % frequency smoothing window increase with frequency
%         cfg.toi          = -0.5:0.01:1;   % time window "slides" from -0.5 to 1 sec in steps of 0.01 sec (10 ms)
%         if icond == 1
%             POSThighTFRmult1 = ft_freqanalysis(cfg,data_artifacts_free);
%         elseif icond == 2
%             POSThighTFRmult2 = ft_freqanalysis(cfg,data_artifacts_free);
%         elseif icond == 3
%             POSThighTFRmult3 = ft_freqanalysis(cfg,data_artifacts_free);
%         elseif icond == 4
%             POSThighTFRmultALL = ft_freqanalysis(cfg,data_artifacts_free);
%         end
%     end
    
    EEGdataset.data_artifacts_free(s) = {data_artifacts_free};
%     EEGdataset.data_artifacts_free(s) = {data_no_artifacts};
    EEGdataset.TFALF.congruent(s) = {POSTlowTFRhann1};
    EEGdataset.TFALF.intermediate(s) = {POSTlowTFRhann2};
    EEGdataset.TFALF.incongruent(s) = {POSTlowTFRhann3};
    EEGdataset.TFALF.all(s) = {POSTlowTFRhannALL};
%     EEGdataset.TFAHF.congruent(s) = {POSThighTFRmult1};
%     EEGdataset.TFAHF.intermediate(s) = {POSThighTFRmult2};
%     EEGdataset.TFAHF.incongruent(s) = {POSThighTFRmult3};
%     EEGdataset.TFAHF.all(s) = {POSThighTFRmultALL};
    EEGdataset.trlxcond(s) = {trlXcond};
%     EEGdataset.TFALF.congruent(s) = {POSTlowFOURIER1};
    
    disp(['###',num2str(s), '###']);
    clearvars -except twin inputDir outDir subjnum s EEGdataset
end

%% Estimating data loss

for s = 1:subjnum
    loss(s,1) = EEGdataset.trlxcond{s}(1,2);
    loss(s,2) = EEGdataset.trlxcond{s}(2,2);
    loss(s,3) = EEGdataset.trlxcond{s}(3,2);
    loss(s,4) = 135-(loss(s,1) + loss(s,2) + loss(s,3));
end

lossperc = sum(loss(:,4))*100/(135*29);
freecond1 = mean(loss(:,1));
freecond2 = mean(loss(:,2));
freecond3 = mean(loss(:,3));
avgloss1 = 45 - mean(loss(:,1));
avgloss2 = 45 - mean(loss(:,2));
avgloss3 = 45 - mean(loss(:,3));
stdloss1 = std(loss(:,1));
stdloss2 = std(loss(:,2));
stdloss3 = std(loss(:,3));

condperc = loss(:,1:3);
for i = 1:29; ranova_input(i,:) = [condperc(i,1), 1, i]; ranova_input(i+29,:) = [condperc(i, 2), 2, i]; ranova_input(i+58,:) = [condperc(i, 3), 3, i]; end
[F1, P1, RMAOV1out] = RMAOV1(ranova_input, 0.05);
t_overview = struct('losspercentage', lossperc,'m_congruent', freecond1, 'm_intermediate', freecond2, 'm_incongruent', freecond3, 's_congruent', stdloss1, 's_intermediate', stdloss2, 's_incongruent', stdloss3, 'fvalue', F1, 'pvalue', P1, 'ranova1', {RMAOV1out});
fprintf('In the time window following the ##TARGET##: \nThe AVERAGE LOSS is : %f \nThe averaged trial number in the ##CONGRUENT## CONDITION is: %f ; STD : %f \nThe averaged trial number in the ##INTERMEDIATE## CONDITION is: %f ; STD: %f \nThe averaged trial number in the ##INCONGRUENT## CONDITION is: %f ; STD : %f \nThe One-Way repeated measures ANOVA for the differences in trial number \nbetween conditions resulted in the following ##F-VALUE##: %f and ##P-VALUE## : %f \n' , t_overview.losspercentage,  t_overview.m_congruent, t_overview.s_congruent, t_overview.m_intermediate, t_overview.s_intermediate, t_overview.m_incongruent, t_overview.s_incongruent, t_overview.fvalue, t_overview.pvalue);

clearvars -except twin inputDir outDir subjnum s EEGdataset t_overview

%% Grand averaged TFAs Low Frequency per condition

cfg = [];
cfg.channel        = 'all';
cfg.parameter      = 'powspctrm'; 
cfg.foilim         = [2 30];
cfg.toilim         = [-0.5 1];
GTFAL_C            = ft_freqgrandaverage(cfg, EEGdataset.TFALF.congruent{:}); 
GTFAL_INT          = ft_freqgrandaverage(cfg, EEGdataset.TFALF.intermediate{:}); 
GTFAL_IC           = ft_freqgrandaverage(cfg, EEGdataset.TFALF.incongruent{:});
GTFAL_ALL          = ft_freqgrandaverage(cfg, EEGdataset.TFALF.all{:});

%% moving to output directory

cd(inputDir)

%% Multiplot TFR Low Frequencies (2-30 Hz) per condition
    
condname1 = 'CONGRUENT';
condname2 = 'INTERMEDIATE';
condname3 = 'INCONGRUENT';

for icond = 2
    figure;
    cfg = [];
    cfg.baseline     = [-0.5 -0.15]; 
    cfg.baselinetype = 'relchange';
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
%     cfg.colormap     = jet(148);
    cfg.zlim         = [-3e-1 4e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
    cfg.fontsize     = 12;
    cfg.showlabels   = 'yes';
    cfg.interactive  = 'yes';
    cfg.comment      = ' ';
    cfg.layout       = 'mpi_customized_acticap64_fran.mat';
    if icond == 1
        condlabel = 'C';
    elseif icond == 2
        condlabel = 'INT';
    else
        condlabel = 'IC';
    end
    ft_multiplotTFR(cfg, eval(['GTFAL_', condlabel]));
    title(['TFR LF multiplot ', eval(['condname',num2str(icond)])]);
end
    
    %% Singleplot on C36 TFR Low Frequencies (2-30 Hz) per condition
    
figure;
for icond = 1:3   
    cfg = [];
    cfg.baseline     = [-0.50 -0.15];
    cfg.baselinetype = 'relchange';  
    cfg.maskstyle    = 'opacity';
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
%     cfg.colormap     = jet(148);
    cfg.fontsize     = 12;
%     cfg.zlim         = [-3e-1 4e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
    cfg.zlim         = [-3e-1 3e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
    cfg.channel      = {'C28', 'C41', 'C35', 'C48', 'C36', 'C42'};
    cfg.colorbar       = 'no';
    if icond == 1
        figure; ft_singleplotTFR(cfg, GTFAL_C);
%         h1 = subplot(1,3,1); ft_singleplotTFR(cfg, GTFAL_C);
%         originalSize1 = get(gca, 'Position');
%         title({' ','TFR(2-30 Hz) C36',['\color{black}', condname1]});
        title('TFR \color{black}CONGRUENT');
    elseif icond == 2
        figure; ft_singleplotTFR(cfg, GTFAL_INT);
%         h2 = subplot(1,3,2); ft_singleplotTFR(cfg, GTFAL_INT);
%         originalSize2 = get(gca, 'Position');
%         title({' ','TFR(2-30 Hz) C36',['\color{blue}', condname2]});
        title('TFR \color{blue}INTERMEDIATE');
    elseif icond == 3
        figure; ft_singleplotTFR(cfg, GTFAL_IC);
        title('TFR \color{red}INCONGRUENT');
    end
%         cfg.colorbar       = 'no';
%         h3 = subplot(1,3,3); ft_singleplotTFR(cfg, GTFAL_IC);
%         title({' ','TFR(2-30 Hz) C36',['\color{red}', condname3]});
%         originalSize3 = get(gca, 'Position');
        c = colorbar;
        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
        c.Label.String = 'rel. power change';
        set(c.Label,'Rotation',270);
        c.Label.FontSize = 17;
        set(c,'YTick', -0.3:0.6:0.3);
        a = get(gca,'YTickLabel');
        set(gca,'YTickLabel',a,'fontsize',15)
%     end
    hold on;
    yL = get(gca, 'YLim');
    plot([0 0], yL, 'k-');
    set(gca,'YTick', 5:5:30);
%     set(eval(['h',num2str(icond)]),'YTick', 5:5:30);
%     set(eval(['h',num2str(icond)]), 'Position', eval(['originalSize',num2str(icond)]));
end

%% Grand averaged TFAs High Frequency per condition

cfg = [];
cfg.channel      = 'all';
cfg.parameter    = 'powspctrm'; 
cfg.foilim       = [25 80];
cfg.toilim       = [-0.5 1];
GTFAH_C          = ft_freqgrandaverage(cfg, EEGdataset.TFAHF.congruent{:}); 
GTFAH_INT        = ft_freqgrandaverage(cfg, EEGdataset.TFAHF.intermediate{:}); 
GTFAH_IC         = ft_freqgrandaverage(cfg, EEGdataset.TFAHF.incongruent{:});
GTFAH_ALL        = ft_freqgrandaverage(cfg, EEGdataset.TFAHF.all{:});

%% Multiplot TFR High Frequencies (2-30 Hz) per condition
    
for icond = 1:3
   figure;
   cfg = [];
   cfg.baseline     = [-0.5 -0.15]; 
   cfg.baselinetype = 'relchange';
   cfg.colormap     = jet(148);
   cfg.ylim         = [30 80];
%    cfg.zlim         = [-2e-1 3e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
   cfg.fontsize     = 12;
   cfg.showlabels   = 'yes';
   cfg.interactive  = 'yes';
   cfg.comment      = ' ';
   cfg.layout       = 'mpi_customized_acticap64_fran.mat';
   if icond == 1
        condlabel = 'C';
    elseif icond == 2
        condlabel = 'INT';
    else
        condlabel = 'IC';
    end
   ft_multiplotTFR(cfg, eval(['GTFAH_', condlabel]));
   title(['TFR multiplot ', eval(['condname',num2str(icond)])]);
end
    
%% Singleplot on C30 TFR High Frequencies (25-80 Hz) per condition
     
figure;
for icond = 3:3   
   cfg = [];
   cfg.baseline     = [-0.50 -0.15];
   cfg.baselinetype = 'relchange';  
   cfg.maskstyle    = 'opacity';
   cfg.colormap     = jet(148);
   cfg.fontsize     = 12;
   cfg.ylim         = [30 80];
   cfg.zlim         = [-3e-1 3e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
   cfg.channel      = {'C28', 'C41', 'C40', 'C36', 'C35', 'C42'};
   cfg.colorbar       = 'no';
   if icond == 1
       figure; ft_singleplotTFR(cfg, GTFAH_C);
%        h4 = subplot(1,3,1); ft_singleplotTFR(cfg, GTFAH_C);
%        originalSize1 = get(gca, 'Position');
       title({' ', 'TFR(30-80 Hz) C36',['\color{black}', condname1]});
   elseif icond == 2
       figure; ft_singleplotTFR(cfg, GTFAH_INT);
%        h5 = subplot(1,3,2); ft_singleplotTFR(cfg, GTFAH_INT);
%        originalSize2 = get(gca, 'Position');
       title({' ', 'TFR(30-80 Hz) C36',['\color{blue}', condname2]});
   elseif icond == 3
       figure; ft_singleplotTFR(cfg, GTFAH_IC);
       title({' ', 'TFR(30-80 Hz) C36',['\color{red}', condname3]});
%        h6 = subplot(1,3,3); ft_singleplotTFR(cfg, GTFAH_IC);
%        originalSize3 = get(gca, 'Position');
   end
       c = colorbar;
       c.Label.String = 'rel. power change';
       set(c.Label,'Rotation',270);
       c.Label.FontSize = 15;
       set(c,'YTick', -0.2:0.5:0.3);

%    end
   hold on;
   yL = get(gca, 'YLim');
   plot([0 0], yL, 'k-');
   set(gca,'YTick', 30:10:80);
%    set(eval(['h',num2str(icond+3)]),'YTick', 30:10:80);
%    set(eval(['h',num2str(icond+3)]), 'Position', eval(['originalSize',num2str(icond)]));
end

%% Contrasts between condition for low frequency range (2-30 Hz)

cfg = [];
cfg.parameter    = 'powspctrm';
cfg.operation    = '(x1-x2)/x3'; %'(x1/x3)-(x2/x3)'; %'(x1-x2)'; 
% TFRL_HCvsLC       = ft_math(cfg, GTFAL_C, GTFAL_INT, GTFAL_ALL);
% TFRL_LCvsIC       = ft_math(cfg, GTFAL_INT, GTFAL_IC, GTFAL_ALL);
% TFRL_HCvsIC       = ft_math(cfg, GTFAL_C, GTFAL_IC, GTFAL_ALL); 
% cfg = [];
% cfg.parameter    = 'powspctrm';
% cfg.operation    = '(x1-x2)'; %'(x1-x2)'; 
TFRL_LCvsHC       = ft_math(cfg, GTFAL_INT, GTFAL_C, GTFAL_ALL);
TFRL_ICvsLC       = ft_math(cfg, GTFAL_IC, GTFAL_INT, GTFAL_ALL);
TFRL_ICvsHC       = ft_math(cfg, GTFAL_IC, GTFAL_C, GTFAL_ALL);


figure;
for icond = 2   
   cfg = [];
   cfg.marker  = 'on';  
   cfg.maskstyle    = 'opacity';
   ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
   colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
%    cfg.colormap     = jet(148);
   cfg.fontsize     = 12;
   cfg.ylim         = [2 30];
   cfg.zlim         = [-1.5e-1 1.5e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
   %cfg.zlim         = 'maxmin';
   cfg.channel      = {'C28', 'C41', 'C35', 'C48', 'C36', 'C42'};
   cfg.fontsize = 20;
   cfg.colorbar       = 'no';
   if icond == 1
       figure; ft_singleplotTFR(cfg, TFRL_LCvsHC);
%        h4 = subplot(1,3,1); ft_singleplotTFR(cfg, TFRL_HCvsLC);
%        originalSize1 = get(gca, 'Position');
       title({'TFR LCvsHC'});
   elseif icond == 2
       figure;  ft_singleplotTFR(cfg, TFRL_ICvsLC);
%        h5 = subplot(1,3,2); ft_singleplotTFR(cfg, TFRL_LCvsIC);
%        originalSize2 = get(gca, 'Position');
       title({'TFR \color{red}IC \color{black}vs \color{blue}LC'});
   elseif icond == 3
       figure;  ft_singleplotTFR(cfg, TFRL_ICvsHC);
%        h6 = subplot(1,3,3); ft_singleplotTFR(cfg, TFRL_HCvsIC);
%        originalSize3 = get(gca, 'Position');
       title({'TFR ICvsHC'});
   end
   c = colorbar;
   ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
   colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
   c.Label.String = 'rel. power change';
   set(c.Label,'Rotation',270);
   c.Label.FontSize = 20;
   set(c,'YTick', -0.15:0.3:0.15);
   hold on;
   yL = get(gca, 'YLim');
   plot([0 0], yL, 'k-');
   set(gca,'XTick', -0.5:0.5:1);
   a = get(gca,'XTickLabel');
   set(gca,'XTickLabel',a,'fontsize',20)
   set(gca,'YTick', 5:5:30);
   b = get(gca,'YTickLabel');
   set(gca,'YTickLabel',b,'fontsize',20)   
%    set(eval(['h',num2str(icond+3)]),'YTick', 5:5:30);
%    set(eval(['h',num2str(icond+3)]), 'Position', eval(['originalSize',num2str(icond)]));
end

clear h1 h2 h3 h4 h5 h6 icond c ans cfg originalSize1 originalSize2 originalSize3 yL

%% Multiplot Contrasts TFR Low Frequencies (2-30 Hz)
    
for icond = 1:3
   figure;
   cfg = [];
   cfg.marker       = 'on';
   cfg.colormap     = jet(148);
   cfg.ylim         = [2 30];
%    cfg.zlim         = [-2e-1 3e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
   cfg.zlim         = 'maxmin';   
   cfg.fontsize     = 12;
   cfg.showlabels   = 'yes';
   cfg.interactive  = 'yes';
   cfg.comment      = ' ';
   cfg.layout       = 'mpi_customized_acticap64_fran.mat';
   if icond == 1
        condlabel = 'LCvsHC';
    elseif icond == 2
        condlabel = 'ICvsHC';
    else
        condlabel = 'ICvsLC';
    end
   ft_multiplotTFR(cfg, eval(['TFRL_', condlabel]));
   title(['TFR multiplot ', condlabel]);
end


%% Contrasts between condition for high frequency range (30-80 Hz)


cfg = [];
cfg.parameter    = 'powspctrm';
cfg.operation    = '(x1-x2)/x3'; %'(x1/x3)-(x2/x3)';
TFRH_HCvsLC       = ft_math(cfg, GTFAH_C, GTFAH_INT, GTFAH_ALL);
TFRH_LCvsIC       = ft_math(cfg, GTFAH_INT, GTFAH_IC, GTFAH_ALL);
TFRH_HCvsIC       = ft_math(cfg, GTFAH_C, GTFAH_IC, GTFAH_ALL); 
% cfg = [];
% cfg.parameter    = 'powspctrm';
% cfg.operation    = '(x1-x2)'; 
% TFRH_LCvsHC       = ft_math(cfg, GTFAH_INT, GTFAH_C);
% TFRH_ICvsLC       = ft_math(cfg, GTFAH_IC, GTFAH_INT);
% TFRH_ICvsHC       = ft_math(cfg, GTFAH_IC, GTFAH_C); 

figure;
for icond = 1:3   
   cfg = [];
   cfg.marker  = 'on';
   cfg.maskstyle    = 'opacity';
%    cfg.colormap     = jet(148);
   ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
   colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
   cfg.fontsize     = 12;
   cfg.ylim         = [30 80];
   cfg.zlim         = [-1.5e-1 1.5e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
   %cfg.zlim         = 'maxmin';
   cfg.channel      = {'C28', 'C41', 'C40', 'C36', 'C35', 'C42'};
   cfg.colorbar       = 'no';
   if icond == 1
       figure; ft_singleplotTFR(cfg, TFRH_HCvsLC);
%        h4 = subplot(1,3,1); ft_singleplotTFR(cfg, TFRH_HCvsLC);
%        originalSize1 = get(gca, 'Position');
       title({'TFR 30-80 Hz HCvsLC C36'});
   elseif icond == 2
       figure; ft_singleplotTFR(cfg, TFRH_LCvsIC);
%        h5 = subplot(1,3,2); ft_singleplotTFR(cfg, TFRH_LCvsIC);
%        originalSize2 = get(gca, 'Position');
       title({'TFR 30-80 Hz LCvsIC C36'});
   elseif icond == 3
       figure; ft_singleplotTFR(cfg, TFRH_HCvsIC);
%        h6 = subplot(1,3,3); ft_singleplotTFR(cfg, TFRH_HCvsIC);
%        originalSize3 = get(gca, 'Position');
       title({'TFR 30-80 Hz HCvsIC C36'});
   end
   c = colorbar;
   c.Label.String = 'rel. power change';
   set(c.Label,'Rotation',270);
   c.Label.FontSize = 15;
   set(c,'YTick', -0.15:0.3:0.15);
   hold on;
   yL = get(gca, 'YLim');
   plot([0 0], yL, 'k-');
   set(gca,'YTick', 30:10:80);
%    set(eval(['h',num2str(icond+3)]),'YTick', 30:10:80);
%    set(eval(['h',num2str(icond+3)]), 'Position', eval(['originalSize',num2str(icond)]));
end

clear h1 h2 h3 h4 h5 h6 icond c ans cfg originalSize1 originalSize2 originalSize3 yL

%% Multiplot Contrasts TFR High Frequencies (30-80 Hz)
    
for icond = 1:3
   figure;
   cfg = [];
   cfg.marker       = 'on';
   cfg.colormap     = jet(148);
   cfg.ylim         = [2 30];
%    cfg.zlim         = [-2e-1 3e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
   cfg.zlim         = 'maxmin';
   cfg.fontsize     = 12;
   cfg.showlabels   = 'yes';
   cfg.interactive  = 'yes';
   cfg.comment      = ' ';
   cfg.layout       = 'mpi_customized_acticap64_fran.mat';
   if icond == 1
        condlabel = 'LCvsHC';
    elseif icond == 2
        condlabel = 'ICvsHC';
    else
        condlabel = 'ICvsLC';
    end
   ft_multiplotTFR(cfg, eval(['TFRL_', condlabel]));
   title(['TFR multiplot ', condlabel]);
end

%% cluster based permtation test

for ifrq = 1:1
    
%     % specifies with which other sensors can form clusters

    % creating clusters based on neighbours DISTANCE in millimeters
    
    cfg_neighb               = [];
    cfg_neighb.neighbourdist = 59; %61/1000 57/1000
    cfg_neighb.channel       = 'all';
    cfg_neighb.enableedit    = 'yes';
    cfg_neighb.layout        = 'mpi_customized_acticap64_fran.mat';
    cfg_neighb.method        = 'distance';
    cfg.feedback             = 'yes';

    % creating clusters based on pre-defined TEMPLATES
    
%     cfg_neighb               = [];
% %     cfg_neighb.enableedit    = 'yes';
%     cfg_neighb.method        = 'template';
%     cfg_neighb.channel       =  'all';
%     cfg_neighb.template      = 'francesco9_59_neighb.mat'; % 5 8888
% %     cfg_neighb.feedback      = 'yes';
%     cfg_neighb.layout        = 'mpi_customized_acticap64_fran.mat';

    if ifrq == 1
        neighbours = ft_prepare_neighbours(cfg_neighb, GTFAL_C);
%         ft_neighbourplot(cfg_neighb, GTFAL_C)
    else
        neighbours = ft_prepare_neighbours(cfg_neighb, GTFAH_C);
%          ft_neighbourplot(cfg_neighb, GTFAH_C)
    end
    
    %estimating neighborhood size
    
    for i=1:length(neighbours)
        neighborhood(i,:) = size(neighbours(i).neighblabel,1);
    end

    mean(neighborhood)

    %%% all time and frequency points from critical word onset to 1000 ms after word onset 
    %%% were submitted to the cluster-based permutation tests, 
    %%% without preselection of time windows, frequency bands, or channels.
    cfg                  =  [];
    cfg.channel          =  'all';
    cfg.latency          =  [0 1];
    if ifrq ==1
        cfg.frequency        =  [2 30];
    else
        cfg.frequency        =  [30 80];
    end
    cfg.parameter        =  'powspctrm';
    cfg.method           =  'montecarlo';
    cfg.statistic        =  'ft_statfun_depsamplesT'; 

    cfg.correctm         =  'cluster';
    cfg.neighbours       =  neighbours;
    cfg.clusteralpha     =  0.025; %  thresholding per tail
    %%% Maris et al. 2007 ("this threshold does not affect the FA rate of the statistical test.
    %%% However, this threshold does affect the sensitivity of the test. For
    %%% example, weak but long-lasting effects are not detected when the threshold
    %%% is large")
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization =  1000; 
    %%% Warning: The p-value confidence interval of positive cluster #1 includes 
    %%% 0.025 - consider increasing the number of permutations! 
    %%% double this number if it turns out that the p-value differs
    %%% from the critical alpha-level (0.05 or 0.01) by less than 0.02.
    cfg.minnbchan        =  2;
    %%% minimum number of neighborhood channels that is
    %%% required for a selected sample
    cfg.clustertail      =  0; % two-sided test
    %cfg.clustercritval   =  2.048; % for parametric thresholding
    cfg.clusterthreshold =  'parametric'; %'nonparametric_individual'; 'nonparametric_common';

    %%% parametric : uses the known parametric distribution of the statistic,
    %%% this has to be supported by the statfun.
    %%% nonparametric_common : estimates a threshold from the randomization
    %%% distribution. The threshold is common to all channe-time-frequency
    %%% points.
    %%% nonparametric_individual : estimates an individual threshold from the
    %%% randomization distribution for each channe-time-frequency point.

    cfg.computecritval   =  'yes';
    cfg.alpha            =  0.05;
    cfg.tail             =  0;           % two-sided test
    cfg.correcttail      =  'no';
    %%% Actually, cfg.alpha is superfluous. In the output,
    %%% every cluster has a p-value assigned to it, and if there is one or more
    %%% cluster with a p-value less than your critical alpha-level (mostly 0.05 for
    %%% a one-sided and 0.025 for a two-sided test), then you have found a
    %%% significant difference.

    %%% It does affect the output.mask field, which can be used for plotting. 
    %%% This will have a 0 for a time-frequency-channel triplet if its p-value 
    %%% is above cfg.alpha, and a 1 when it is below it.

    % experimental design: WHITHIN SUBJECTS experiment
    cfg.design(1,1:2*subjnum)  = [ones(1,subjnum) 2*ones(1,subjnum)];
    cfg.design(2,1:2*subjnum)  = [1:subjnum 1:subjnum];
    cfg.uvar     = 2;
    cfg.ivar     = 1;

    for cid = 1:3
        if cid == 1
            if     ifrq == 1
                 [LstatLCvsHC] = ft_freqstatistics(cfg, EEGdataset.TFALF.intermediate{:}, EEGdataset.TFALF.congruent{:}); %
%                    [LstatHCvsLC] = ft_freqstatistics(cfg, EEGdataset.TFALF.congruent{:}, EEGdataset.TFALF.intermediate{:}); %
            elseif ifrq == 2
%                  [HstatLCvsHC] = ft_freqstatistics(cfg, EEGdataset.TFAHF.intermediate{:}, EEGdataset.TFAHF.congruent{:}); %
                   [HstatHCvsLC] = ft_freqstatistics(cfg, EEGdataset.TFAHF.congruent{:}, EEGdataset.TFAHF.intermediate{:}); %
            end
        elseif cid == 2
            if     ifrq == 1
                 [LstatICvsLC] = ft_freqstatistics(cfg, EEGdataset.TFALF.incongruent{:}, EEGdataset.TFALF.intermediate{:}); % 
%                    [LstatLCvsIC] = ft_freqstatistics(cfg, EEGdataset.TFALF.intermediate{:}, EEGdataset.TFALF.incongruent{:}); % 
            elseif ifrq == 2
%                  [HstatICvsLC] = ft_freqstatistics(cfg, EEGdataset.TFAHF.incongruent{:}, EEGdataset.TFAHF.intermediate{:}); % 
                   [HstatLCvsIC] = ft_freqstatistics(cfg, EEGdataset.TFAHF.intermediate{:}, EEGdataset.TFAHF.incongruent{:}); % 
            end
        elseif cid == 3
            if     ifrq == 1
                 [LstatICvsHC] = ft_freqstatistics(cfg, EEGdataset.TFALF.incongruent{:}, EEGdataset.TFALF.congruent{:}); % 
%                    [LstatHCvsIC] = ft_freqstatistics(cfg, EEGdataset.TFALF.congruent{:}, EEGdataset.TFALF.incongruent{:}); % 
            elseif ifrq == 2 
%                  [HstatICvsHC] = ft_freqstatistics(cfg, EEGdataset.TFAHF.incongruent{:},  EEGdataset.TFAHF.congruent{:}); %
                   [HstatHCvsIC] = ft_freqstatistics(cfg, EEGdataset.TFAHF.congruent{:}, EEGdataset.TFAHF.incongruent{:}); %
            end
        end
    end
end

clear neighborhood neighbours cfg_neighb ans cfg

for cid = 1:3
    if cid ==1
        %condlabel = 'LstatHCvsLC';
        condlabel = 'LstatLCvsHC';
    elseif cid == 2
%         condlabel = 'LstatLCvsIC';
        condlabel = 'LstatICvsLC';
    elseif cid == 3
%         condlabel = 'LstatHCvsIC';
        condlabel = 'LstatICvsHC';
    elseif cid == 4
        condlabel = 'HstatHCvsLC';
    elseif cid == 5
        condlabel = 'HstatLCvsIC';
    elseif cid == 6
        condlabel = 'HstatHCvsIC';
    end
    p_min = min(min(min(eval([condlabel, '.prob']))));
    fprintf('In the following contrast: %s the lowest pvalue observed was: %f \n', condlabel(6:end), p_min);
    if p_min < 0.025
        fprintf('There was ##SOME## significant effect in %s \n', condlabel)
    else
        fprintf('There was ##NO## significant effect in %s \n', condlabel)
    end
end 

clear p_min i ifrq condlabel cid

%{
In the following contrast: LCvsHC the lowest pvalue observed was: 0.023976 
There was ##SOME## significant effect in LstatLCvsHC 
In the following contrast: ICvsLC the lowest pvalue observed was: 0.009990 
There was ##SOME## significant effect in LstatICvsLC 
In the following contrast: ICvsHC the lowest pvalue observed was: 0.057942 
There was ##NO## significant effect in LstatICvsHC
%}

%% plot significant clusters

for ifrq = 1:1
    for icond = 1:3
       figure;
       cfg = [];
       cfg.marker        = 'on';
       cfg.parameter     = 'stat';  % plot the t-value 
       cfg.mask          = 'yes';
       cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
       cfg.maskstyle     = 'saturation';
       cfg.colormap      = jet(148);
       if ifrq ==1
           cfg.ylim          = [2 30];
       else 
           cfg.ylim          = [30 80];
       end
       %cfg.zlim         = [-2e-1 3e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
       cfg.zlim          = 'maxmin';   
       cfg.fontsize      = 12;
       cfg.channel       = 'C*';
       cfg.showlabels    = 'yes';
       cfg.interactive   = 'yes';
       cfg.comment       = ' ';
       cfg.layout        = 'mpi_customized_acticap64_fran.mat';
       if icond == 1
            condlabel = 'HCvsLC';
       elseif icond == 2
            condlabel = 'HCvsIC';
       else
            condlabel = 'LCvsIC';
       end
       if ifrq == 1 
           ft_multiplotTFR(cfg, eval(['Lstat', condlabel]));
       else
           ft_multiplotTFR(cfg, eval(['Hstat', condlabel]));
       end
       title(['TFR multiplot ', condlabel]);
    end
end

figure;
for icond = 1   
   cfg = [];
%    if icond == 2
%        TFRL_HCvsIC.mask = logical(zeros(size(TFRL_HCvsIC.powspctrm)));
%        TFRL_HCvsIC.mask(:,:,51:151) = LstatHCvsIC.mask(:,:,:);
%    elseif icond == 1
%        TFRH_LCvsIC.mask = logical(zeros(size(TFRH_LCvsIC.powspctrm)));
%        TFRH_LCvsIC.mask(:,3:23,51:151) = HstatLCvsIC.mask(:,:,:);
%    elseif icond == 3
%        TFRL_HCvsLC.mask = logical(zeros(size(TFRL_HCvsLC.powspctrm)));
%        TFRL_HCvsLC.mask(:,:,51:151) = LstatHCvsLC.mask(:,:,:);
%    end
   cfg.marker        = 'on';
   cfg.parameter     = 'stat';  % plot the t-value 
   cfg.mask          = 'yes';
   cfg.maskparameter = 'posclusterslabelmat';  % use the thresholded probability to mask the data
%    cfg.maskparameter = 'mask';
   cfg.maskstyle     = 'saturation';
   cfg.colormap      = jet(148);
%    cfg.colormap      = brewermap(64,'RdBu');
   cfg.fontsize      = 25;
   cfg.xlim          = [0 1];
   if icond == 1
       cfg.ylim          = [2 30];
   else
       cfg.ylim          = [30 80];
   end
   cfg.zlim          = [-3 3]; % asymmetric color scales can be used to highlight increases or decreases in activity.
%    cfg.zlim         = 'maxmin';
   cfg.channel       = 'C48';
   cfg.colorbar      = 'no';
%    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
   if icond == 1
       figure; ft_singleplotTFR(cfg, LstatICvsLC);
%        h4 = subplot(1,3,1); ft_singleplotTFR(cfg, TFRH_LCvsIC);
%        originalSize1 = get(gca, 'Position');
%        title({'TFR 30-80 Hz ICvsLC C36'});
         title('C48');
   elseif icond == 3
       figure; ft_singleplotTFR(cfg, HstatHCvsLC);
%        h6 = subplot(1,3,2); ft_singleplotTFR(cfg, TFRL_HCvsLC);
%        originalSize3 = get(gca, 'Position');
%        title({'TFR 2-30 Hz LCvsHC C36'});
         title({' '});
   elseif icond == 2
       figure; ft_singleplotTFR(cfg, LstatHCvsIC);
%        h5 = subplot(1,3,3); ft_singleplotTFR(cfg, TFRL_HCvsIC);
%        originalSize2 = get(gca, 'Position');
%        title({'TFR 2-30 Hz ICvsHC C36'});
         title({' '});
   end
%    c = colorbar;
%    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
%    c.Label.String = 't-value';
%    set(c.Label,'Rotation',270);
%    c.Label.FontSize = 20;
%    set(c,'YTick', -3:6:3);
   hold on;
%    yL = get(gca, 'YLim');
%    plot([0 0], yL, 'k-');
   if icond == 1
       set(gca,'YTick', 2:28:30);
       b = get(gca,'YTickLabel');
       set(gca,'YTickLabel',b,'fontsize',20)
       set(gca,'XTick', 0:1:1);
       a = get(gca,'XTickLabel');
       set(gca,'XTickLabel',a,'fontsize',20)
%        set(eval(['h',num2str(icond+3)]),'YTick', 5:5:30);
   else
       set(gca,'YTick', 30:10:80);
%        set(eval(['h',num2str(icond+3)]),'YTick', 30:10:80);
   end
%    set(eval(['h',num2str(icond+3)]), 'Position', eval(['originalSize',num2str(icond)]));
end

%% topoplot contrasts

% surface Laplacian (FT_SCALPCURRENTDENSITY)
% will increase topographical localization and highlight local spatial features.

toi1 = [0.8 1];
toi2 = [0.1 0.4];
foi1 = [2 5];
foi2 = [70 80];

for fig=1:1
    
    cfg = [];
    %cfg.marker                = 'on';  
    cfg.layout                 = 'mpi_customized_acticap64_fran.mat';
    cfg.highlight              = 'on';
    cfg.highlightsymbol        = '*';
    cfg.highlightcolor         = 'k';
    cfg.highlightsize          = 15;
    cfg.baselinetype           = 'relchange';
    cfg.comment                = ' ';
    cfg.renderer               = 'painters';
    cfg.style                  = 'straight';
    cfg.gridscale              = 250; 
    cfg.parameter              = 'powspctrm';
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
%     cfg.colormap               = jet(148);
    set(gca,'FontSize',25);
    
    if fig == 1
       cfg.xlim             = toi1;
       cfg.ylim             = foi1;
       cfg.zlim             = [-1.5e-1 1.5e-1]; 
       contrast             = 'LCvsIC';
       filename             = 'thetaeffect2-30Hz.png';
       for i = 1:56
       testch = find(abs(LstatICvsLC.stat(i,1:4,71:101))>3);
       if ~isempty(testch); ch(i,1) = i; else ch(i,1) = 0; end
       end
       cfg.highlightchannel = ch(find(ch~=0));
%        cfg.highlightchannel = cfg.highlightchannel(1:17);
    elseif fig == 2
       cfg.xlim             = toi2;
       cfg.ylim             = foi2;
       cfg.zlim             = [-1.5e-1 1.5e-1]; 
       contrast             = 'LCvsIC';
       filename             = 'gammaeffect30-80Hz.png';
       for i = 1:56
       testch = find(abs(HstatLCvsIC.stat(i,17:21,11:41))>2.8);
       if ~isempty(testch); ch(i,1) = i; else ch(i,1) = 0; end
       end
       cfg.highlightchannel = ch(find(ch~=0));
    end
    timewin = ([num2str(cfg.xlim(1)*1000), '-', num2str(cfg.xlim(2)*1000),' ms']);
    freqwin = ([num2str(cfg.ylim(1)), '-', num2str(cfg.ylim(2)), 'Hz']);
    if fig == 1
        figure; ft_topoplotTFR(cfg, TFRL_ICvsLC);
        set(gca,'FontSize',17);
        title({'800-1000 ms'});%{'DELTA-THETA EFFECT', 'INCONGRUENT - INTERMEDIATE'})
%         text(5, 1.5, {freqwin, timewin},'FontSize',15)
        c = colorbar;
        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
        set(c,'YTick', -0.15:0.3:0.15);
    elseif fig == 2  
        figure; ft_topoplotTFR(cfg, TFRH_LCvsIC);
        set(gca,'FontSize',17);
        title({'GAMMA EFFECT', 'INTERMEDIATE - INCONGRUENT'})
        text(5, 1.5, {freqwin, timewin},'FontSize',15)
        c = colorbar;
        set(c,'YTick', -0.15:0.3:0.15);
    end
    c.Label.String = 'rel. power change';
    set(c.Label,'Rotation',270);
    c.Label.FontSize = 19;
    c.Location = 'EastOutside';
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',17)
    disp('###Press a key to go on with the next plot###')  % Press a key here
    pause;
    print(filename, '-dpng')
    clf
    close
end

%{
deltatheta t value >3.2
gamma t value>3.1
%}

%% Power Time Course plot

cfg = [];
cfg.baseline     = [-0.5 -0.15];
cfg.baselinetype = 'relchange';
TFRL_basC        = ft_freqbaseline(cfg, GTFAL_C);
TFRL_basINT      = ft_freqbaseline(cfg, GTFAL_INT);
TFRL_basIC       = ft_freqbaseline(cfg, GTFAL_IC);
% TFRH_basC        = ft_freqbaseline(cfg, GTFAH_C);
% TFRH_basINT      = ft_freqbaseline(cfg, GTFAH_INT);
% TFRH_basIC       = ft_freqbaseline(cfg, GTFAH_IC);

cfg = [];
cfg.baseline     = [-0.5 -0.15];
cfg.baselinetype = 'db';
for sid = 1:subjnum
    newfreqcongL   = ft_freqbaseline(cfg, EEGdataset.TFALF.congruent{1, sid});
    newfreqinterL  = ft_freqbaseline(cfg, EEGdataset.TFALF.intermediate{1, sid});
    newfreqincongL = ft_freqbaseline(cfg, EEGdataset.TFALF.incongruent{1, sid});
    EEGdataset.baselineLF.congruent(sid) = {newfreqcongL};
    EEGdataset.baselineLF.intermediate(sid) = {newfreqinterL};
    EEGdataset.baselineLF.incongruent(sid) = {newfreqincongL};
%     newfreqcongH   = ft_freqbaseline(cfg, EEGdataset.TFAHF.congruent{1, sid});
%     newfreqinterH  = ft_freqbaseline(cfg, EEGdataset.TFAHF.intermediate{1, sid});
%     newfreqincongH = ft_freqbaseline(cfg, EEGdataset.TFAHF.incongruent{1, sid});
%     EEGdataset.baselineHF.congruent(sid) = {newfreqcongH};
%     EEGdataset.baselineHF.intermediate(sid) = {newfreqinterH};
%     EEGdataset.baselineHF.incongruent(sid) = {newfreqincongH};
    clear newfreqcongL newfreqinterL newfreqincongL newfreqcongH newfreqinterH newfreqincongH
end

for sid = 1:subjnum
    slicedeltathetaC(sid,:) = mean(EEGdataset.baselineLF.congruent{1, sid}.powspctrm(27, 1:4, 51:150),2);
    slicedeltathetaINT(sid,:) = mean(EEGdataset.baselineLF.intermediate{1, sid}.powspctrm(27, 1:4, 51:150),2);
    slicedeltathetaIC(sid,:) = mean(EEGdataset.baselineLF.incongruent{1, sid}.powspctrm(27, 1:4, 51:150),2);
%     slicealphaC(sid,:) = mean(EEGdataset.baselineLF.congruent{1, sid}.powspctrm(35, 7:11, 51:150),2);
%     slicealphaINT(sid,:) = mean(EEGdataset.baselineLF.intermediate{1, sid}.powspctrm(35, 7:11, 51:150),2);
%     slicealphaIC(sid,:) = mean(EEGdataset.baselineLF.incongruent{1, sid}.powspctrm(35, 7:11, 51:150),2);
%     slicegammaC(sid,:) = mean(EEGdataset.baselineHF.congruent{1, sid}.powspctrm(35, 19:23, 51:150),2);
%     slicegammaINT(sid,:) = mean(EEGdataset.baselineHF.intermediate{1, sid}.powspctrm(35, 19:23, 51:150),2);
%     slicegammaIC(sid,:) = mean(EEGdataset.baselineHF.incongruent{1, sid}.powspctrm(35, 19:23, 51:150),2);
end

% Calculating SEM and 95% CIs
N = 29; % number of subjects

M1 = mean(slicedeltathetaC,1);
E1 = (std(slicedeltathetaC,1)./sqrt(N)); %*1.96; % standard error of the mean for congruent
M2 = mean(slicedeltathetaINT,1);
E2 = (std(slicedeltathetaINT,1)./sqrt(N)); %*1.96; % standard error of the mean for C30 intermediate
M3 = mean(slicedeltathetaIC,1);
E3 = (std(slicedeltathetaIC,1)./sqrt(N)); %*1.96; % standard error of the mean for C30 incongruent
% M1 = mean(slicealphaC,1);
% E1 = (std(slicealphaC,1)./sqrt(N)); %*1.96; % standard error of the mean for congruent
% M2 = mean(slicealphaINT,1);
% E2 = (std(slicealphaINT,1)./sqrt(N)); %*1.96; % standard error of the mean for C30 intermediate
% M3 = mean(slicealphaIC,1);
% E3 = (std(slicealphaIC,1)./sqrt(N)); %*1.96; % standard error of the mean for C30 incongruent
% M1 = mean(slicegammaC,1);
% E1 = (std(slicegammaC,1)./sqrt(N)); %*1.96; % standard error of the mean for congruent
% M2 = mean(slicegammaINT,1);
% E2 = (std(slicegammaINT,1)./sqrt(N)); %*1.96; % standard error of the mean for C30 intermediate
% M3 = mean(slicegammaIC,1);
% E3 = (std(slicegammaIC,1)./sqrt(N)); %*1.96; % standard error of the mean for C30 incongruent


% figure settings
figure;
title('Power time course Delta-Theta (2-5 Hz) C29');
box off
xlim([0 1])
% ylim([0.25 0.28])
set(gca,'FontSize',16);
set(gca,'TickDir','out');
% set(gca,'YTick',(0.2:0.01:0.28));
set(gca,'XTick',(0:0.2:1));
ylabel('Power [dB]');
xlabel('Time [s]');

% Plotting CIs filled polygons
% M1 = smooth(M1,10); % smoothing the mean
% M1 = M1'; % re-transposing the mean in an horizontal vector
p1 = patch([GTFAL_C.time(51:150),GTFAL_C.time(150):-0.01:GTFAL_C.time(50)],[M1+E1,M1(length(GTFAL_C.time(51:150)):-1:1)-E1(length(GTFAL_C.time(51:150)):-1:1)], [0.3 0.3 0.3]); % creating a Confidence interval for MLT12 (meaningful channel)
set(p1,'EdgeColor','none');
alpha(p1,0.3);
hold on
theta1 = plot(GTFAL_C.time(51:150), M1, 'k','LineWidth', 3);
% M2 = smooth(M2,10);
% M2 = M2';
p2 = patch([GTFAL_INT.time(51:150),GTFAL_INT.time(150):-0.01:GTFAL_INT.time(50)],[M2+E2,M2(length(GTFAL_INT.time(51:150)):-1:1)-E2(length(GTFAL_INT.time(51:150)):-1:1)], [0.2 0.4 1]); % creating a Confidence interval for MRF41 (flat channel)
set(p2,'EdgeColor','none');
alpha(p2,0.3);
hold on
theta2 = plot(GTFAL_INT.time(51:150), M2, 'b','LineWidth', 3);
% M3 = smooth(M3,10);
% M3 = M3';
p3 = patch([GTFAL_IC.time(51:150),GTFAL_IC.time(150):-0.01:GTFAL_IC.time(50)],[M3+E3,M3(length(GTFAL_IC.time(51:150)):-1:1)-E3(length(GTFAL_IC.time(51:150)):-1:1)], [1 0.3 0.3]); % creating a Confidence interval for the GApc (absolute baseline)
set(p3,'EdgeColor','none');
alpha(p3,0.3);
hold on
theta3 = plot(GTFAL_IC.time(51:150), M3, 'r','LineWidth', 3);

legend([theta1,theta2,theta3],{'congruent \pm s.e.m','intermediate \pm s.e.m.','incongruent \pm s.e.m.'}, 'Location', 'northwest');

%% Topoplot per condition Delta-Theta & Gamma

toi1 = [0.1 0.9];
toi2 = [0.1 0.9];
toi3 = [0.1 0.9];
foi1 = [2 5];
foi2 = [70 80];
foi3 = [8 12];

for fig=1:2
    
    cfg = [];
    %cfg.marker                = 'on';  
    cfg.layout                 = 'mpi_customized_acticap64_fran.mat';
%     cfg.highlight              = 'on';
%     cfg.highlightsymbol        = '.';
%     cfg.highlightcolor         = 'k';
%     cfg.highlightsize          = 30;
%     cfg.highlightchannel       = 35;
    cfg.baselinetype           = 'relchange';
    cfg.comment                = ' ';
    cfg.renderer               = 'painters';
    cfg.style                  = 'straight';
    cfg.gridscale              = 250; 
    cfg.parameter              = 'powspctrm';
    cfg.colormap               = jet(148);
    if fig == 1
       cfg.xlim             = toi1;
       cfg.ylim             = foi1;
       cfg.zlim             = [-1.5e-1 1.5e-1]; 
       filename             = 'thetaconditions.png';       
    elseif fig == 2
       cfg.xlim             = toi2;
       cfg.ylim             = foi2;
       cfg.zlim             = [-1.5e-1 1.5e-1]; 
       filename             = 'gammaconditions.png';
    elseif fig == 3
       cfg.xlim             = toi3;
       cfg.ylim             = foi3;
       cfg.zlim             = [-2e-1 2e-1]; 
       filename             = 'alphaconditions.png';
        
    end
    timewin = ([num2str(cfg.xlim(1)*1000), '-', num2str(cfg.xlim(2)*1000),' ms']);
    freqwin = ([num2str(cfg.ylim(1)), '-', num2str(cfg.ylim(2)), 'Hz']);
    if fig == 1 || fig == 3
        figure; 
        subplot(1,3,1); ft_topoplotTFR(cfg, TFRL_basC);
        set(gca,'FontSize',15);
        title('\color{black}congruent')
        subplot(1,3,2); ft_topoplotTFR(cfg, TFRL_basINT);
        set(gca,'FontSize',15);
        title('\color{blue}intermediate')
        h3 = subplot(1,3,3); ft_topoplotTFR(cfg, TFRL_basIC);
        set(gca,'FontSize',15);
        title('\color{red}incongruent')
        c = colorbar;
        set(c,'YTick', -0.15:0.3:0.15);
    elseif fig == 2
        figure;
        subplot(1,3,1); ft_topoplotTFR(cfg, TFRH_basC);
        set(gca,'FontSize',15);
        title('\color{black}congruent')
        subplot(1,3,2); ft_topoplotTFR(cfg, TFRH_basINT);
        set(gca,'FontSize',15);
        title('\color{blue}intermediate')
        h3 = subplot(1,3,3); ft_topoplotTFR(cfg, TFRH_basIC);
        set(gca,'FontSize',15);
        title('\color{red}incongruent')
        c = colorbar;
        set(c,'YTick', -0.15:0.3:0.15);
    end
    originalSize3 = get(gca, 'Position');
    set(h3, 'Position', originalSize3);
    c.Label.String = 'rel. power change';
    set(c.Label,'Rotation',270);
    c.Label.FontSize = 15;
    c.Location = 'EastOutside';
    set(c,'position',[.93 .35 .02 .27]);
    if fig ==1
        t = suptitle('DELTA-THETA band (2-5 Hz) 100-900 ms');
        set(t,'FontSize',17);
    elseif fig ==2
        t = suptitle('GAMMA band (70-80 Hz) 100-900 ms');
        set(t,'FontSize',17);
    elseif fig == 3
        t = suptitle('ALPHA band (8-12 Hz) 100-900 ms');
        set(t,'FontSize',17);
    end
    disp('###Press a key to go on with the next plot###')  % Press a key here
    pause;
    print(filename, '-dpng')
    clf
    close
end


cfg=[];
cfg.freqlow    = [2 5];
cfg.freqhigh   = [70 80];
cfg.channel    = 'C35';
cfg.method     = 'plv';
crossfreq = ft_crossfrequencyanalysis(cfg, GTFAL_C, GTFAH_C);
