function headpoint = read_head_pos(pos_file)
% Read coordinates of fiducials from .pos file
% 
% - Input
%  pos_file : Position file (.pos)
%
% Y. Takeda 2018-09-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

fid = fopen(pos_file);
c = textscan(fid,'%s %f %f %f');
fclose(fid);
pos = [c{2} c{3} c{4}];
headpoint = pos(pos(:,3)>0,:);% Exclude points below ears