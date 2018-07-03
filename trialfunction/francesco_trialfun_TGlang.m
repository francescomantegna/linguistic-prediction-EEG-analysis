%{
    file_name : francesco_trialfun_TGlang
    author : Francesco Mantegna
    project : Music&Poetry
    date : 27/12/2017
%}

function [trl, event] = francesco_trialfun_TGlang(cfg)

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
    counter2 = 46;
    counter3 = 91;
    for j = 1:(length(value))
        trg1 = char(value(j));
        try
        trg5 = char(value(j+4));
        catch
            warning('Problem using function. Matrix exceeds dimensions.');
        end
        if strcmp(trg1,'S  7') && strcmp(trg5,'S 27')
            trlbegin = sample(j+4) + pretrig;       
            trlend   = sample(j+4) + posttrig;       
            offset   = pretrig;
            condnum = 1;
            trlnumb = char(value(j-1));
            trlnumb = regexp(trlnumb,'\d*','Match');
            trlnumb = str2double(trlnumb);
            congtarg   = [trlbegin trlend offset condnum trlnumb];
            trl(counter1,1:5)   = [congtarg];
            counter1 = counter1 + 1;
        elseif strcmp(trg1,'S  8') && strcmp(trg5,'S 27')
            trlbegin = sample(j+4) + pretrig;       
            trlend   = sample(j+4) + posttrig;       
            offset   = pretrig;
            condnum = 2;
            trlnumb = char(value(j-1));
            trlnumb = regexp(trlnumb,'\d*','Match');
            trlnumb = str2double(trlnumb);
            intertarg   = [trlbegin trlend offset condnum trlnumb];
            trl(counter2,1:5)   = [intertarg];
            counter2 = counter2 + 1;
        elseif strcmp(trg1,'S  9') && strcmp(trg5,'S 27')
            trlbegin = sample(j+4) + pretrig;       
            trlend   = sample(j+4) + posttrig;       
            offset   = pretrig;
            condnum = 3;
            trlnumb = char(value(j-1));
            trlnumb = regexp(trlnumb,'\d*','Match');
            trlnumb = str2double(trlnumb);
            incongtarg   = [trlbegin trlend offset condnum trlnumb];
            trl(counter3,1:5)   = [incongtarg];
            counter3 = counter3 + 1;
        end
    end
end