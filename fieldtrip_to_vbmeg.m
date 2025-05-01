function DATA = fieldtrip_to_vbmeg(fieldtrip_file, vbmeg_file, transmat)
% Convert FieldTrip file into VBMEG file
% [Usage]
%    DATA = fieldtrip_to_vbmeg(fieldtrip_file, vbmeg_file)
% [INPUT]
%    fieldtrip_file  (from) FieldTrip file(.mat)
%    vbmeg_file      (to)   VBMEG file(.meg.mat or .eeg.mat)
%    transmat               This matrix is required to transform input coordinate to VBMEG coordinate system (SPM_Right_m)
%                               Neuromag: transform matrix from MEG head coordinate to MRI RAS coordinate
%
% [OUTPUT]
%    DATA             VBMEG data structure
% [NOTE]
%

[~,tmp,~] = fileparts(vbmeg_file);
[~,~,output_type] = fileparts(tmp); % 'meg' or 'eeg'

flag_transe = false;
if exist('transmat','var') && ~isempty(transmat)
    flag_transe = true;
end

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

switch output_type
    case '.meg'
        % Extract MEG, REFERENCE and others (EEG is ignored)
        DATA.bexp     = data(ix_meg,:,:);
        DATA.refmg    = data(ix_ref,:,:);
        DATA.bexp_ext = data(ix_ext,:,:);

        % pick and Qpick in original coordinate system
        pick  = header.grad.coilpos;
        Qpick = header.grad.coilori;

        % Convert unit to [m]
        coef = inner_get_unitcoef(header.grad.unit);
        pick = pick .* coef;

        % Coordinate transformation to SPM_Right_m
        if flag_transe
            switch upper(header.grad.coordsys)
                case 'NEUROMAG'
                    % Transform original coordinate to MRI RAS coordinate system
                    pick  = trans_coord(pick, transmat);
                    Qpick = trans_coord(Qpick, transmat, 1);% Only rotate and normalize
                    coord_system = 'SPM_Right_m';
                otherwise
                    warning(['Currently ' upper(header.grad.coordsys) ' coordinate cannot be transformed to VBMEG coordinate'])
                    flag_transe = false;
            end
        end

        if ~flag_transe
            warning('Original coordinate system is kept (but unit is [m])')
            coord_system = [upper(header.grad.coordsys) '_m']; % Keep original coord [m]
        end

        DATA.pick  = pick ;
        DATA.Qpick = Qpick;
        DATA.CoordType = coord_system;
        DATA.Measurement = 'MEG';

        % Which device is used?
        dev_type = header.grad.type;
        if contains(upper(dev_type), 'NEUROMAG')
            device = 'NEUROMAG';
            is_opm = false;
        else
            error(['Currently ' dev_type ' is not supported yet'])
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
            error('Append here')
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
        % Extract EEG and others (MEGs are ignored)
        eeg_data = data(ix_eeg,:,:);
        ext_data = data(ix_ext,:,:);
        DATA.eeg_data = cat(1, eeg_data, ext_data);
        Nchannel     = size(eeg_data, 1);
        Nchannel_ext = size(ext_data, 1);

        % Set measurement
        DATA.Measurement = 'EEG';

        % Coord(pos) in original coordinate system
        Coord = header.elec.chanpos;

        % Convert unit to [m]
        coef  = inner_get_unitcoef(header.elec.unit);
        Coord = Coord .* coef;

        % Check 'coordsys' info in elec
        warning('Now checking coordinate system of elec object...')
        if ~isfield(header.elec, 'coordsys') || isempty(header.elec.coordsys)
            warning('The coordinate system of elec object is not specified');
            % If elec has no 'coordsys', then refer to grad
            use_grad_coordsys = false;
            warning('Reffering to grad object...')
            if isfield(header, 'grad') && ~isempty(header.grad)
                if isfield(header.grad, 'coordsys') && ~isempty(header.grad.coordsys)
                    use_grad_coordsys = true;
                    warning('Use the coordinate system of grad object')
                else
                    warning('No coordinate system in grad object')
                end
            else
                warning('No grad object in the header')
            end
            if use_grad_coordsys
                % Use same coordinate system as MEG
                elec_coordsys = header.grad.coordsys;
                warning(['Specified coordinate system: ' elec_coordsys])
            else
                % No coordinate info in the header
                flag_transe = false;
                elec_coordsys = [];
                warning('No coordinate transformation')
            end
        else
            elec_coordsys = header.elec.coordsys;
            warning(['Specified coordinate system: ' elec_coordsys])
       end

        % Coordinate transformation to SPM_Right_m
        if flag_transe
            switch upper(elec_coordsys)
                case 'NEUROMAG'
                    % Transform original coordinate to MRI RAS coordinate system
                    Coord  = trans_coord(Coord, transmat);
                    coord_system = 'SPM_Right_m';
                otherwise
                    warning(['Currently ' upper(header.grad.coordsys) ' coordinate cannot be transformed to VBMEG coordinate'])
                    flag_transe = false;
            end
        end

        if ~flag_transe
            warning('Original coordinate system is kept (but unit is [m])')
            coord_system = [upper(header.grad.coordsys) '_m']; % Keep original coord [m]
        end

        % Which device is used?
        if isfield(header.elec, 'type') && ~isempty(header.elec.type)
            dev_type = header.elec.type;
        else
            dev_type = header.grad.type;
        end

        if contains(upper(dev_type), 'NEUROMAG')
            % device = 'NEUROMAG';
            device = 'BRAINAMP'; % Use data format for BRAINAMP
            warning(['Device is registered as ' device ' for convenience...'])
        else
            error(['Currently ' dev_type ' is not supported yet'])
        end

        % Init EEGinfo
        EEGinfo = [];
        EEGinfo.Measurement = 'EEG';

        % Set Coord and CoordType
        EEGinfo.Coord  = Coord;
        EEGinfo.CoordType = coord_system;

        % EEGinfo
        EEGinfo.Nchannel = Nchannel; % Num of eeg channels
        EEGinfo.Nsample  = size(eeg_data,2);
        EEGinfo.Nrepeat  = size(eeg_data,3);
        EEGinfo.Pretrigger = 0; % header.nSamplesPre is ignored
        EEGinfo.SampleFrequency = header.Fs;
        EEGinfo.ChannelID   = [1:EEGinfo.Nchannel]';
        EEGinfo.ChannelName = header.elec.label;
        EEGinfo.MRI_ID = [];
        EEGinfo.Device = device;
        EEGinfo.ActiveChannel = ones(EEGinfo.Nchannel,1); % All active
        EEGinfo.ActiveTrial   = ones(EEGinfo.Nrepeat,1); % All active
        EEGinfo.cond_id  = ones(1,EEGinfo.Nrepeat); % Modify after imporing data using orig.event
        EEGinfo.DataType = repmat({'float32'},Nchannel+Nchannel_ext,1); % Set float32 to all(eeg and extra)

        % Make ChannelInfo
        ChannelInfo = [];
        ChannelInfo.Name   = EEGinfo.ChannelName;
        ChannelInfo.ID     = EEGinfo.ChannelID;
        ChannelInfo.Active = EEGinfo.ActiveChannel;
        ChannelInfo.Type   = ones(EEGinfo.Nchannel,1); % All 1 (electrode)
        ChannelInfo.PhysicalUnit = header.elec.chanunit;
 
        % Make ExtraChannelInfo
        ExtraChannelInfo.Channel_name   = header.label(ix_ext);
        ExtraChannelInfo.Channel_id     = [EEGinfo.ChannelID(end)+1:EEGinfo.ChannelID(end)+Nchannel_ext]';
        ExtraChannelInfo.Channel_active = ones(Nchannel_ext,1); % All active
        temp = zeros(size(data,1),1);
        temp(ix_trigger) = -1;
        temp(ix_ecg) = -3;
        temp(ix_eog) = -5;
        temp(ix_etc) = -4;
        ExtraChannelInfo.Channel_type = temp(ix_ext);
        chanunit = header.chanunit(ix_ext);
        ExtraChannelInfo.PhysicalUnit = chanunit;
        EEGinfo.ExtraChannelInfo = ExtraChannelInfo;

        % Make Trial struct
        if is_epoched
            Trial = [];
            for tt=1:EEGinfo.Nrepeat
                Trial(tt).sample = [event.sample(tt):event.sample(tt)+EEGinfo.Nsample];
                Trial(tt).number = tt;
                Trial(tt).Active = 1; % All active
            end
        else
            Trial.sample = [1:EEGinfo.Nsample];
            Trial.number = 1;
            Trial.Active = 1;
        end
        EEGinfo.Trial = Trial;

        DATA.EEGinfo = EEGinfo;

        % Keep original header and event
        orig = [];
        orig.format = 'fieldtrip';
        orig.event  = event;
        orig.header = header;
        DATA.orig = orig;

        % Save in VBMEG format
        save(vbmeg_file, '-struct', 'DATA');

    otherwise
        error('Output vbmeg file must be ''meg.mat'' or ''eeg.mat''')
end

end

function coef = inner_get_unitcoef(unit)
% Unit coefficient for correcting to [m]
switch unit
    case 'mm'
        coef = 1/1000;
    case 'cm'
        coef = 1/100;
    case 'm'
        coef = 1;
end
end
