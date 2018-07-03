%{
    file_name : francesco_trialfun_ERPmusic
    author : Francesco Mantegna
    project : Music&Poetry
    date : 27/12/2017
%}

function [trl, event] = francesco_trialfun_ERPmusic(cfg);

    % read the header information and the events from the data
    hdr   = ft_read_header(cfg.headerfile);
    event = ft_read_event(cfg.datafile);

    % search for "trigger" events
    value  = {event(find(strcmp('Stimulus', {event.type}))).value}';
    sample = [event(find(strcmp('Stimulus', {event.type}))).sample]';

    % determine the number of samples before and after the trigger
    pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
    posttrig =  round(cfg.trialdef.post * hdr.Fs);

    % look for the combination of a trigger "7" followed by a trigger "64" 
    % for each trigger except the last one
    counter1 = 1;
    counter2 = 49;
    counter3 = 97;
    for j = 1:(length(value))
        trg1 = char(value(j));
        try
        trg5 = char(value(j+5));
        catch
            warning('Problem using function. Matrix exceeds dimensions.');
        end
        if strcmp(trg1,'S  7') & strcmp(trg5,'S 27')
            trlbegin = sample(j+5) + pretrig;       
            trlend   = sample(j+5) + posttrig;       
            offset   = pretrig;
            condnum = 1;
            congtarg   = [trlbegin trlend offset condnum];
            trl(counter1,1:4)   = [congtarg];
            counter1 = counter1 + 1;
        elseif strcmp(trg1,'S  8') & strcmp(trg5,'S 27')
            trlbegin = sample(j+5) + pretrig;       
            trlend   = sample(j+5) + posttrig;       
            offset   = pretrig;
            condnum = 2;
            intertarg   = [trlbegin trlend offset condnum];
            trl(counter2,1:4)   = [intertarg];
            counter2 = counter2 + 1;
        elseif strcmp(trg1,'S  9') & strcmp(trg5,'S 27')
            trlbegin = sample(j+5) + pretrig;       
            trlend   = sample(j+5) + posttrig;       
            offset   = pretrig;
            condnum = 3;
            incongtarg   = [trlbegin trlend offset condnum];
            trl(counter3,1:4)   = [incongtarg];
            counter3 = counter3 + 1;
        end
    end
end