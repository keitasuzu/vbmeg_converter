function [lpa, nasion, rpa] = read_fids_from_json(file)
% Read coordinates of fiducials from .json file
%
% - Input
% file : Json file describing position of fiducials (LPA, Nasion, RPA)
%
% - Output
% lpa : Position of LPA
% nasion : Position of nasion
% rpa : Position of RPA
%
% Y. Takeda 2018-09-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

fid = fopen(file, 'rt');

nasion = [];
lpa    = [];
rpa    = [];
while(1)
    s = fgets(fid);
    if s == -1, break; end
    [coord, name] = read_fid(s);
    if ~isempty(name)
        switch(name)
            case 'LPA'
                lpa    = coord;
            case 'Nasion'
                nasion = coord;
            case 'RPA'
                rpa    = coord;
        end
    end
end
fclose(fid);

function [coord, name] = read_fid(s)

if strfind(s, 'LPA')
    start_ix = strfind(s, '[');
    end_ix   = strfind(s, ']');
    coord    = eval(s(start_ix:end_ix));
    name     = 'LPA';
elseif strfind(s, 'Nasion')
    start_ix = strfind(s, '[');
    end_ix   = strfind(s, ']');
    coord    = eval(s(start_ix:end_ix));
    name     = 'Nasion';
elseif strfind(s, 'RPA')
    start_ix = strfind(s, '[');
    end_ix   = strfind(s, ']');
    coord    = eval(s(start_ix:end_ix));
    name     = 'RPA';
else
    coord = []; name = [];
end

