function V1 = trans_coord(V0, trans_matrix, rot_flg)
% Transform coordinates using a transformation matrix
%
% - Input
%  V0 : Coordinates (Npoints x 3)
%  trans_matrix : Transformation matrix
%  rot_fle : Rotation flag. If 1, rotation and normalization will be done
%
% - Output
%  V1 : Transformed coordinates (Npoints x 3)
%
% Y. Takeda 2018-09-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

if ~exist('rot_flg','var')
    rot_flg = 0;
end

if rot_flg == 0
    V1 = [V0  ones(size(V0,1),1)]*trans_matrix;
elseif rot_flg ==1% For Qpick
    V1 = V0 * trans_matrix(1:3, 1:3);
    
    % normalization
    V1 = vb_repmultiply(V1, 1./sqrt(sum(V1.^2, 2)) );
end