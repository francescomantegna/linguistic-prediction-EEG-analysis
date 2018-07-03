%{
    file_name : TFALangPRE.m
    author : Francesco Mantegna
    affiliation : Max Planck Institut for Psycholinguistics
    project : Music&Poetry
    date : 14/02/2018
%}

%% add fieldtrip path

full_path_to_fieldtrip2017 = '/Users/francesco/Documents/MATLAB/fieldtrip-20170414';
addpath(full_path_to_fieldtrip2017)
ft_defaults

%% defining variables

twin       =    [-1.8 0.3];
inputDir1  =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/Preprocessing/Data0.1cutoffPREP';
inputDir2  =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/Preprocessing/Data0.1cutoffPRETG';
outDir     =    '/Users/francesco/Documents/MATLAB/MPI/Music&Poetry/OSCanalysis/ResultsPRESTIMULUS2';
subjnum    =    30; % subject 6, 9, 24 were excluded

%% loading pre-processed files before PRIME onset

cd(inputDir1)

for s = 1:subjnum
    
    load(['preprolangOSCPREP_SUBJ' int2str(s) '.mat']);
    
    % exclude peripherical occipital electrodes C22 C26 C54
    cfg                      = []; 
    cfg.channel              = setdiff(1:59, [22 26 53]);       
    data_artifacts_free  = ft_selectdata(cfg, data_no_artifacts);
    
    % Time-frequency analysis post-stimulus for each condition High Frequencies (25-80 Hz)
    
%     cfg                = [];
%     cfg.output         = 'pow';
%     cfg.channel        = 'C*';
%     cfg.method         = 'mtmconvol';
%     cfg.trials         = 'all';   
%     cfg.foi            = 25:2.5:80;                         % analysis 30 to 80 Hz in steps of 1 Hz 
%     cfg.t_ftimwin      = ones(length(cfg.foi),1).*0.4;   % length of time smoothing window = 0.5 sec
%     cfg.pad            = ceil(max(cellfun(@numel, data_artifacts_free.time)/data_artifacts_free.fsample));
%     %cfg.tapsmofrq      = 5;   % fixed 5 Hz frequency smoothing window 
%     cfg.tapsmofrq      = 0.2*cfg.foi; % frequency smoothing window increase with frequency
%     cfg.toi            = -1.5:0.01:0;   % time window "slides" from -0.5 to 1 sec in steps of 0.01 sec (10 ms)
%     PREPhighTFRmultNEW = ft_freqanalysis(cfg,data_artifacts_free);
    
    dataset.prime.data_artifacts_free(s) = {data_artifacts_free};
    dataset.prime.TFALF(s) = {PREPlowTFRhann};
%     dataset.prime.TFAHF(s) = {PREPhighTFRmultNEW};
    dataset.prime.trlxcond(s) = {trlXcond};
    
    disp(['###',num2str(s), '###']);
    clearvars -except twin inputDir1 inputDir2 outDir subjnum s dataset
end

%% loading pre-processed files before PRIME onset

cd(inputDir2)

for s = 1:subjnum
    
    load(['preprolangOSCPRETG_SUBJ' int2str(s) '.mat']);
    
    % exclude peripherical occipital electrodes C22 C26 C54
    cfg                      = []; 
    cfg.channel              = setdiff(1:59, [22 26 53]);      
    data_artifacts_free  = ft_selectdata(cfg, data_no_artifacts);
    
    % Time-frequency analysis post-stimulus for each condition High Frequencies (25-80 Hz)
    
%         cfg                = [];
%         cfg.output         = 'pow';
%         cfg.channel        = 'C*';
%         cfg.method         = 'mtmconvol';
%         cfg.trials         = 'all';   
%         cfg.foi            = 25:2.5:80;                         % analysis 30 to 80 Hz in steps of 1 Hz 
%         cfg.t_ftimwin      = ones(length(cfg.foi),1).*0.4;   % length of time smoothing window = 0.5 sec
%         cfg.pad            = ceil(max(cellfun(@numel, data_artifacts_free.time)/data_artifacts_free.fsample));
%         %cfg.tapsmofrq      = 5;   % fixed 5 Hz frequency smoothing window 
%         cfg.tapsmofrq      = 0.2*cfg.foi; % frequency smoothing window increase with frequency
%         cfg.toi            = -1.5:0.01:0;   % time window "slides" from -0.5 to 1 sec in steps of 0.01 sec (10 ms)
%         PREThighTFRmultNEW = ft_freqanalysis(cfg,data_artifacts_free);
    
    dataset.target.data_artifacts_free(s) = {data_no_artifacts};
    dataset.target.TFALF(s) = {PRElowTFRhann};
%     dataset.target.TFAHF(s) = {PREThighTFRmultNEW};
    dataset.target.trlxcond(s) = {trlXcond};
    
    disp(['###',num2str(s), '###']);
    clearvars -except twin inputDir outDir subjnum s dataset
end

%% Estimating data loss

for s = 1:subjnum
    ploss(s,1) = dataset.prime.trlxcond{s}(1,2);
    ploss(s,2) = dataset.prime.trlxcond{s}(2,2);
    ploss(s,3) = dataset.prime.trlxcond{s}(3,2);
    ploss(s,4) = 135-(ploss(s,1) + ploss(s,2) + ploss(s,3));
end
for s = 1:subjnum
    tloss(s,1) = dataset.target.trlxcond{s}(1,2);
    tloss(s,2) = dataset.target.trlxcond{s}(2,2);
    tloss(s,3) = dataset.target.trlxcond{s}(3,2);
    tloss(s,4) = 135-(tloss(s,1) + tloss(s,2) + tloss(s,3));
end

for i = 1:2
    if i == 1
        loss = ploss;
    else
        loss = tloss;
    end
    
    lossperc = sum(loss(:,4))*100/(135*29);
    freecond1 = mean(loss(:,1));
    freecond2 = mean(loss(:,2));
    freecond3 = mean(loss(:,3));
    stdloss1 = std(loss(:,1));
    stdloss2 = std(loss(:,2));
    stdloss3 = std(loss(:,3));

    condperc = loss(:,1:3);
    for idt = 1:28; ranova_input(idt,:) = [condperc(idt,1), 1, idt]; ranova_input(idt+28,:) = [condperc(idt, 2), 2, idt]; ranova_input(idt+56,:) = [condperc(idt, 3), 3, idt]; end
    [F1, P1, RMAOV1out] = RMAOV1(ranova_input, 0.05);
    
    if i == 1
        p_overview = struct('losspercentage', lossperc,'m_congruent', freecond1, 'm_intermediate', freecond2, 'm_incongruent', freecond3, 's_congruent', stdloss1, 's_intermediate', stdloss2, 's_incongruent', stdloss3, 'fvalue', F1, 'pvalue', P1, 'ranova1', {RMAOV1out});
        fprintf('In the time window preceding the ##PRIME##: \nThe AVERAGE LOSS is : %f \nThe averaged trial number in the ##CONGRUENT## CONDITION is: %f ; STD : %f \nThe averaged trial number in the ##INTERMEDIATE## CONDITION is: %f ; STD: %f \nThe averaged trial number in the ##INCONGRUENT## CONDITION is: %f ; STD : %f \nThe One-Way repeated measures ANOVA for the differences in trial number \nbetween conditions resulted in the following ##F-VALUE##: %f and ##P-VALUE## : %f \n' , p_overview.losspercentage,  p_overview.m_congruent, p_overview.s_congruent, p_overview.m_intermediate, p_overview.s_intermediate, p_overview.m_incongruent, p_overview.s_incongruent, p_overview.fvalue, p_overview.pvalue); 
    else
        t_overview = struct('losspercentage', lossperc,'m_congruent', freecond1, 'm_intermediate', freecond2, 'm_incongruent', freecond3, 's_congruent', stdloss1, 's_intermediate', stdloss2, 's_incongruent', stdloss3, 'fvalue', F1, 'pvalue', P1, 'ranova1', {RMAOV1out});
        fprintf('In the time window following the ##TARGET##: \nThe AVERAGE LOSS is : %f \nThe averaged trial number in the ##CONGRUENT## CONDITION is: %f ; STD : %f \nThe averaged trial number in the ##INTERMEDIATE## CONDITION is: %f ; STD: %f \nThe averaged trial number in the ##INCONGRUENT## CONDITION is: %f ; STD : %f \nThe One-Way repeated measures ANOVA for the differences in trial number \nbetween conditions resulted in the following ##F-VALUE##: %f and ##P-VALUE## : %f \n' , t_overview.losspercentage,  t_overview.m_congruent, t_overview.s_congruent, t_overview.m_intermediate, t_overview.s_intermediate, t_overview.m_incongruent, t_overview.s_incongruent, t_overview.fvalue, t_overview.pvalue);
    end
    
    clearvars -except twin inputDir1 inputDir2 outDir subjnum s dataset p_overview t_overview ploss tloss
end

clear ploss tloss

%% Grand averaged TFAs 2-30 Hz before PRIME and TARGET onset

cfg = [];
cfg.channel        = 'all';
cfg.parameter      = 'powspctrm'; 
cfg.foilim         = [2 30];
cfg.toilim         = [-1.5 0];
GAlow_P            = ft_freqgrandaverage(cfg, dataset.prime.TFALF{:});
GAlow_T            = ft_freqgrandaverage(cfg, dataset.target.TFALF{:});

%% Grand averaged TFAs 25-80 Hz before PRIME and TARGET onset

% cfg = [];
% cfg.channel        = 'all';
% cfg.parameter      = 'powspctrm'; 
% cfg.foilim         = [25 80];
% cfg.toilim         = [-1.5 0];
% GAhigh_P          = ft_freqgrandaverage(cfg, dataset.prime.TFAHF{:});
% GAhigh_T          = ft_freqgrandaverage(cfg, dataset.target.TFAHF{:});

%% Contrast between pre-stimulus prime vs target in low frequency range (2-30 Hz)

cfg = [];
cfg.parameter    = 'powspctrm';
cfg.operation    = '(x1-x2)/(x1+x2)';
TvsPlow       = ft_math(cfg, GAlow_T, GAlow_P);

%% Contrast between pre-stimulus prime vs target in high frequency range (25-80 Hz)

% cfg = [];
% cfg.parameter    = 'powspctrm';
% cfg.operation    = '(x1-x2)/(x1+x2)';
% TvsPhigh       = ft_math(cfg, GAhigh_T, GAhigh_P);

%% plot the difference between TARGET and PRIME in the two frequency of interest range (i.e. low 2-30 Hz, high 30-80 Hz)

for ifreq = 1
   figure;
   cfg = [];
   cfg.marker  = 'on';  
%    cfg.maskstyle    = 'opacity';
%    cfg.colormap     = jet(148);
   ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
   colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
   cfg.fontsize     = 25;
   cfg.xlim         = [-1.5 0];
   if ifreq ==1
       cfg.ylim         = [2 30];
   else
       cfg.ylim         = [30 80];
   end
%    cfg.zlim         = 'maxmin';
   cfg.channel      = {'C37', 'C43','C44', 'C38', 'C39', 'C45'};%{'C29', 'C30','C33', 'C35', 'C36', 'C1'};
   if ifreq == 1
       cfg.zlim = [-1.5e-1 1.5e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
       ft_singleplotTFR(cfg, TvsPlow);
       set(gca,'YTick', 5:5:30);
       b = get(gca,'YTickLabel');
       set(gca,'YTickLabel',b,'fontsize',20)
       set(gca,'XTick', -1.5:0.5:0);
       a = get(gca,'XTickLabel');
       set(gca,'XTickLabel',a,'fontsize',20)
       c = colorbar;
       ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
       colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
       c.Label.String = 'rel. power change';
       set(c.Label,'Rotation',270);
       c.Label.FontSize = 20;
       set(c,'YTick', -0.15:0.3:0.15);
       title({'TFR TARGET vs PRIME'});
   elseif ifreq == 2
       cfg.zlim = [-0.8e-1 1.2e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
       ft_singleplotTFR(cfg, TvsPhigh);
       set(gca,'YTick', 30:10:80);
       c = colorbar;
       c.Label.String = 'rel. power change';
       set(c.Label,'Rotation',270);
       c.Label.FontSize = 15;
       set(c,'YTick', -0.08:0.20:0.12);
       title({'TFR TvsP 30-80 Hz C36'});
   end
    hold on;
    yL = get(gca, 'YLim');
    plot([-0.5 -0.5], yL, 'k-');
end


%% change directory to print results

cd(outDir)

% print('TvsPlowC36.png', '-dpng');
% print('TvsPhighC36.png', '-dpng');

%% cluster based permtation test

for ifrq = 1
    
%     % specifies with which other sensors can form clusters

    % creating clusters based on neighbours DISTANCE in millimeters
    
    cfg_neighb               = [];
    cfg_neighb.neighbourdist = 59; %56/1000 || 58/2000 || 59-61/1000
    cfg_neighb.channel       = 'all';
    cfg_neighb.enableedit    = 'yes';
    cfg_neighb.layout        = 'mpi_customized_acticap64_fran.mat';
    cfg_neighb.method        = 'distance';
    cfg.feedback             = 'yes';

    % creating clusters based on pre-defined TEMPLATES
    
%     cfg_neighb               = [];
%     cfg_neighb.method        = 'template';
%     cfg_neighb.channel       =  'all';
%     cfg_neighb.template      = 'francesco2_59_neighb.mat'; % 2 or 4
% %     cfg_neighb.feedback      = 'yes';
%     cfg_neighb.layout        = 'mpi_customized_acticap64_fran.mat';

    if ifrq == 1
        neighbours = ft_prepare_neighbours(cfg_neighb, GAlow_T);
        ft_neighbourplot(cfg_neighb, GAlow_T)
    else
        neighbours = ft_prepare_neighbours(cfg_neighb, GAhigh_T);
        ft_neighbourplot(cfg_neighb, GAhigh_T)
    end
    
    %estimating neighborhood size
    
    for i=1:length(neighbours)
        neighborhood(i,:) = size(neighbours(i).neighblabel,1);
    end

    mean(neighborhood)

    cfg                  =  [];
    cfg.channel          =  'all';
    cfg.latency          =  [-1.5 0];
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
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization =  1000; 
    cfg.minnbchan        =  2;
    cfg.clustertail      =  0; % two-sided test
    %cfg.clustercritval   =  2.048; % for parametric thresholding
    cfg.clusterthreshold =  'parametric'; %'nonparametric_individual'; 'nonparametric_common';

    cfg.computecritval   =  'yes';
    cfg.alpha            =  0.025;
    cfg.tail             =  0;           % two-sided test
    cfg.correcttail      =  'no';

    % experimental design: WHITHIN SUBJECTS experiment
    cfg.design(1,1:2*subjnum)  = [ones(1,subjnum) 2*ones(1,subjnum)];
    cfg.design(2,1:2*subjnum)  = [1:subjnum 1:subjnum];
    cfg.uvar     = 2;
    cfg.ivar     = 1;

    if     ifrq == 1
%            [LstatPvsT] = ft_freqstatistics(cfg, dataset.prime.TFALF{:}, dataset.target.TFALF{:}); %
           [LstatTvsP] = ft_freqstatistics(cfg, dataset.target.TFALF{:}, dataset.prime.TFALF{:}); %
    elseif ifrq == 2
%            [HstatPvsT] = ft_freqstatistics(cfg, dataset.prime.TFAHF{:}, dataset.target.TFAHF{:}); %
           [HstatTvsP] = ft_freqstatistics(cfg, dataset.target.TFAHF{:}, dataset.prime.TFAHF{:}); %
    end
end

for cid = 1
    if cid ==1
        condlabel = 'LstatTvsP';
    elseif cid == 2
        condlabel = 'HstatTvsP';
    end
    p_min = min(min(min(eval([condlabel, '.prob']))));
    fprintf('In the following contrast: %s the lowest pvalue observed was: %f \n', condlabel(6:end), p_min);
    if p_min < 0.025
        fprintf('There was ##SOME## significant effect in %s \n', condlabel)
    else
        fprintf('There was ##NO## significant effect in %s \n', condlabel)
    end
end

%{
In the following contrast: TvsP the lowest pvalue observed was: 0.000999 
There was ##SOME## significant effect in LstatTvsP 
In the following contrast: TvsP the lowest pvalue observed was: 0.009990 
There was ##SOME## significant effect in HstatTvsP 
%}

%% plot significant clusters

for ifrq = 1:1
       figure;
       cfg = [];
       cfg.marker        = 'on';
       cfg.parameter     = 'stat';  % plot the t-value 
       cfg.mask          = 'yes';
       cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
       cfg.maskstyle     = 'saturation';
%        cfg.colormap      = jet(148);
       ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
       colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
       cfg.xlim          = [-1.5 0];
       if ifrq ==1
           cfg.ylim          = [2 30];
       else 
           cfg.ylim          = [30 80];
       end
       %cfg.zlim         = [-2e-1 3e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
%        cfg.zlim          = 'maxmin';   
       cfg.fontsize      = 12;
       cfg.channel       = 'C*';
       cfg.showlabels    = 'yes';
       cfg.interactive   = 'yes';
       cfg.comment       = ' ';
       cfg.layout        = 'mpi_customized_acticap64_fran.mat';
       condlabel         = 'TvsP';
       if ifrq == 1 
           ft_multiplotTFR(cfg, eval(['Lstat', condlabel]));
           title(['TFR multiplot 2-30 Hz ', condlabel]);
       else
           ft_multiplotTFR(cfg, eval(['Hstat', condlabel]));
           title(['TFR multiplot 30-80 Hz ', condlabel]);
       end
end

for ifreq = 1
   figure;
   cfg               = [];
   cfg.marker        = 'on';
   cfg.parameter     = 'stat';  % plot the t-value 
   cfg.mask          = 'yes';
   cfg.fontsize      = 25;
   cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
   cfg.maskstyle     = 'saturation';
   cfg.colormap      = jet(148);
%    cfg.colormap      = brewermap(64,'RdBu');
   cfg.xlim          = [-1.5 0];
   if ifreq ==1
       cfg.ylim         = [2 30];
   else
       cfg.ylim         = [30 80];
   end
%    cfg.zlim         = 'maxmin';
%    cfg.channel      = {'C37', 'C43','C44', 'C38', 'C39', 'C45'};{'C29', 'C30','C33', 'C35', 'C36', 'C1'};
   cfg.channel      = 'C45';
   if ifreq == 1
%        cfg.zlim = [-1.5e-1 1.5e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
       cfg.zlim = [-3 3];
       cfg.colorbar = 'no';
       ft_singleplotTFR(cfg, LstatTvsP);
       set(gca,'YTick', 2:28:30);
       b = get(gca,'YTickLabel');
       set(gca,'YTickLabel',b,'fontsize',20)
       set(gca,'XTick', -1.5:1.5:0);
       a = get(gca,'XTickLabel');
       set(gca,'XTickLabel',a,'fontsize',20)
%        c = colorbar;
%        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%        colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
%        c.Label.String = 't-value';
%        set(c.Label,'Rotation',270);
%        c.Label.FontSize = 15;
%        set(c,'YTick', -0.15:0.15:0.15);
%        set(c,'YTick', -3:6:3);
       title({'C45'});
   elseif ifreq == 2
%        cfg.zlim = [-1.5e-1 1.5e-1]; % asymmetric color scales can be used to highlight increases or decreases in activity.
       cfg.zlim = [-3 3];
       cfg.renderer = 'painters';
       ft_singleplotTFR(cfg, HstatTvsP);
       set(gca,'YTick', 30:10:80);
       c = colorbar;
       c.Label.String = 't-value';
       set(c.Label,'Rotation',270);
       c.Label.FontSize = 15;
%        set(c,'YTick', -0.15:0.15:0.15);
       set(c,'YTick', -3:6:3);
       title({'TFR TARGETvsPRIME 30-80 Hz C36'});
   end
end

%% topoplot contrasts

% surface Laplacian (FT_SCALPCURRENTDENSITY)
% will increase topographical localization and highlight local spatial features.

toi1 = [-0.6 -0.4];
toi2 = [-1.5 0];
toi3 = [-1.3 -0.5];
foi1 = [2 5];
foi2 = [15 25];
foi3 = [30 35];

for fig=2
    
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
%     cfg.colormap               = jet(148);
    cfg.colormap      = brewermap(64,'RdBu');
    set(gca,'FontSize',20);
    contrast             = 'TvsP';    
    if fig == 1
       cfg.xlim             = toi1;
       cfg.ylim             = foi1;
       cfg.zlim             = [-1e-1 1e-1]; 
       filename             = 'TvsPdeltathetaeffect2-5Hz.png';
       for i = 1:56
       testch = find(abs(LstatTvsP.stat(i,1:4,91:111))>4);
       if ~isempty(testch); ch(i,1) = i; else ch(i,1) = 0; end
       end
       cfg.highlightchannel = ch(find(ch~=0));
    elseif fig == 2
       cfg.xlim             = toi2;
       cfg.ylim             = foi2;
       cfg.zlim             = [-1e-1 1e-1];
       filename             = 'TvsPbetaeffect15-25Hz.png';
       for i = 1:56
       testch = find(abs(LstatTvsP.stat(i,14:24,1:151))>3);
       if ~isempty(testch); ch(i,1) = i; else ch(i,1) = 0; end
       end
       cfg.highlightchannel = ch(find(ch~=0));
    elseif fig == 3
       cfg.xlim             = toi3;
       cfg.ylim             = foi3;
       cfg.zlim             = [-1e-1 1e-1];
       filename             = 'TvsPgammaeffect30-35Hz.png';
       for i = 1:56
       testch = find(abs(HstatTvsP.stat(i,1:3,21:101))>4);
       if ~isempty(testch); ch(i,1) = i; else ch(i,1) = 0; end
       end
       cfg.highlightchannel = ch(find(ch~=0));
    end
    timewin = ([num2str(cfg.xlim(1)*1000), '-', num2str(cfg.xlim(2)*1000),' ms']);
    freqwin = ([num2str(cfg.ylim(1)), '-', num2str(cfg.ylim(2)), 'Hz']);
    if fig == 1
        figure; ft_topoplotTFR(cfg, TvsPlow);
        set(gca,'FontSize',17);
        title({'DELTA-THETA EFFECT', 'TARGET-PRIME'})
%         text(5, 1.5, {freqwin, timewin},'FontSize',15)
        c = colorbar;
        set(c,'YTick', -0.1:0.1:0.1);
    elseif fig == 2  
        figure; ft_topoplotTFR(cfg, TvsPlow);
        set(gca,'FontSize',17);
        title({'-1500 - 0 ms'})
%         text(5, 1.5, {freqwin, timewin},'FontSize',15)
        c = colorbar;
        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(64,'RdBu'))) % change the colormap red blue
        set(c,'YTick', -0.1:0.2:0.1);
    elseif fig == 3
        figure; ft_topoplotTFR(cfg, TvsPhigh);
        set(gca,'FontSize',17);
        title({'GAMMA EFFECT', 'TARGET-PRIME'})
        text(5, 1.5, {freqwin, timewin},'FontSize',15)
        c = colorbar;
        set(c,'YTick', -0.1:0.1:0.1);
    end
    c.Label.String = 'rel. power change';
    set(c.Label,'Rotation',270);
    c.Label.FontSize = 18;
    c.Location = 'EastOutside';
    disp('###Press a key to go on with the next plot###')  % Press a key here
    pause;
    print(filename, '-dpng')
    clf
    close
end

cfg = [];
cfg.layout                 = 'mpi_customized_acticap64_fran.mat';
cfg.highlight              = 'on';
cfg.highlightsymbol        = '.';
cfg.highlightcolor         = 'r';
cfg.highlightsize          = 40;
cfg.comment                = ' ';
cfg.style                  = 'blank';
cfg.gridscale              = 250; 
cfg.parameter              = 'powspctrm';
cfg.highlightchannel       =  {'C37', 'C43','C44', 'C38', 'C39', 'C45'};
ft_topoplotTFR(cfg, TvsPlow);
