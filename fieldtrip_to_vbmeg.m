function DATA = fieldtrip_to_vbmeg(fieldtrip_file, vbmeg_file)
% Convert FieldTrip file into VBMEG file
% [Usage]
%    DATA = fieldtrip_to_vbmeg(fieldtrip_file, vbmeg_file)
% [INPUT]
%    fieldtrip_file  (from) FieldTrip file(.mat)
%    vbmeg_file      (to)   VBMEG file(.meg.mat or .eeg.mat)
% [OUTPUT]
%    EEG             EEGLAB data structure
% [NOTE]
%

[~,tmp,~] = fileparts(vbmeg_file);
[~,~,type] = fileparts(tmp); % 'meg' or 'eeg'

% Load FieldTrip data
header = ft_read_header(fieldtrip_file);
data   = ft_read_data(fieldtrip_file);
event  = ft_read_event(fieldtrip_file);

% Check if the data is continuous or epoched
is_epoched = ndims(data)==3;

% Prepare index
% See ft_chantype.m
%--EEG
ix_eeg       = strcmp(header.chantype, 'eeg');
%--MEG
ix_megmag    = strcmp(header.chantype, 'megmag');
ix_meggrad   = strcmp(header.chantype, 'meggrad');
ix_megplanar = strcmp(header.chantype, 'megplanar');
ix_meg       = ix_megmag | ix_meggrad | ix_megplanar;
%--REFERENCE MEG
ix_refmag    = strcmp(header.chantype, 'refmag');
ix_refgrad   = strcmp(header.chantype, 'refgrad');
ix_refplanar = strcmp(header.chantype, 'refplanar');
ix_ref       = ix_refmag | ix_refgrad | ix_refplanar;
%--EXTRA
ix_trigger   = contains(header.chantype, 'trig'); % {'trig','trigger','dtrig'}
ix_ecg       = strcmp(header.chantype, 'ecg');
ix_eog       = contains(header.chantype, 'eog'); % {'eog','veog','heog'}
ix_etc       = strcmp(header.chantype, 'etc') | strcmp(header.chantype, 'misc');
ix_ext       = ix_trigger | ix_ecg | ix_eog | ix_etc;

switch type
    case '.meg'
        % Extract MEG, REFERENCE and others (EEG is ignored)
        DATA.bexp     = data(ix_meg,:,:);
        DATA.refmg    = data(ix_ref,:,:);
        DATA.bexp_ext = data(ix_ext,:,:);

        % Coordinate transformation
        header.grad.coordsys;

        % Unit coefficient for correcting to [m]
        switch header.grad.unit
            case 'mm'
                coef = 1/1000;
            case 'cm'
                coef = 1/100;
            case 'm'
                coef = 1;
        end
        DATA.pick  = header.grad.coilpos .* coef;
        DATA.Qpick = header.grad.coilori;
        DATA.CoordType = 'SPM_Right_m';
        DATA.Measurement = 'MEG';

        % Which device is used?
        dev = header.grad.type;
        if contains(upper(dev), 'NEUROMAG')
            device = 'NEUROMAG';
            is_opm = false;
        else
            error(['Currently ' dev ' is not supported yet'])
        end

        % MEGinfo
        MEGinfo = [];
        MEGinfo.sensor_weight = header.grad.tra;
        MEGinfo.Nchannel = size(DATA.bexp,1);
        MEGinfo.Nsample  = size(DATA.bexp,2);
        MEGinfo.Nrepeat  = size(DATA.bexp,3);
        MEGinfo.Pretrigger = 0; % header.nSamplesPre is ignored
        MEGinfo.SampleFreq = header.Fs;
        MEGinfo.MEGch_id = [1:MEGinfo.Nchannel]';
        MEGinfo.MEGch_name = header.grad.label;
        MEGinfo.MEG_ID = [];
        MEGinfo.MRI_ID = [];
        MEGinfo.device = device;
        MEGinfo.ActiveChannel = ones(MEGinfo.Nchannel,1); % All active
        MEGinfo.ActiveTrial   = ones(MEGinfo.Nrepeat,1); % All active
        MEGinfo.saveman = []; % Data will be saved in the same file
        MEGinfo.cond_id = ones(1,MEGinfo.Nrepeat); % Modify after imporing data using orig.event

        % Make ChannelInfo
        ChannelInfo = [];
        ChannelInfo.Name   = MEGinfo.MEGch_name;
        ChannelInfo.ID     = MEGinfo.MEGch_id;
        ChannelInfo.Active = MEGinfo.ActiveChannel;
        temp = zeros(size(data,1),1);
        if is_opm
            error('Fix here')
        else
            temp(ix_megmag)    = 1;
            temp(ix_meggrad)   = 2;
            temp(ix_megplanar) = 3;
            channel_type = temp(ix_meg);
        end
        ChannelInfo.Type = channel_type;
        MEGinfo.ChannelInfo = ChannelInfo;

        % Make ExtraChannelInfo
        Nchannel_ext = sum(ix_ref | ix_ext);
        ExtraChannelInfo.Channel_name   = header.label(ix_ref | ix_ext);
        ExtraChannelInfo.Channel_id     = [MEGinfo.MEGch_id(end)+1:MEGinfo.MEGch_id(end)+Nchannel_ext]';
        ExtraChannelInfo.Channel_active = ones(Nchannel_ext,1); % All active
        temp = zeros(size(data,1),1);
        temp(ix_refmag)    = 257;
        temp(ix_refgrad)   = 258;
        temp(ix_refplanar) = 259;
        temp(ix_trigger) = -1;
        temp(ix_ecg) = -3;
        temp(ix_eog) = -5;
        temp(ix_etc) = -4;
        ExtraChannelInfo.Channel_type = temp(ix_ref | ix_ext);
        MEGinfo.ExtraChannelInfo = ExtraChannelInfo;

        % Make Trial struct
        if is_epoched
            Trial = [];
            for tt=1:MEGinfo.Nrepeat
                % Trial(tt).sample = [1+MEGinfo.Nsample*(tt-1):MEGinfo.Nsample*tt];
                Trial(tt).sample = [event.sample(tt):event.sample(tt)+MEGinfo.Nsample];
                Trial(tt).number = tt;
                Trial(tt).Active = 1; % All active
            end
        else
            Trial.sample = [1:MEGinfo.Nsample];
            Trial.number = 1;
            Trial.Active = 1;
        end
        MEGinfo.Trial = Trial;

        DATA.MEGinfo = MEGinfo;

        % Keep original header and event
        orig = [];
        orig.format = 'fieldtrip';
        orig.event  = event;
        orig.header = header;
        DATA.orig = orig;

        % Save in VBMEG format
        save(vbmeg_file, '-struct', 'DATA');

    case '.eeg'
        error('Currently eeg.mat is not supported')

    otherwise
        error('Output vbmeg file must be ''meg.mat'' or ''eeg.mat''')
end

