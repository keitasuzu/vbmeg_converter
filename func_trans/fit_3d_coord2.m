function trans_opt = fit_3d_coord2(fids_meg, headpoint_meg, fids_ras, V_headsurface, trans_ini)
% Transform MEG coordinate to MRI RAS coordinate using fiducials, head points, and head surface
%
% - Input
%  fids_meg : Fiducial coordinates in MEG head coordinate system
%  headpoint_meg : Headpoint coordinates in MEG head coordinate system
%  fids_ras : Fiducial coordinates in MRI RAS coordinate system
%  V_headsurface : Coordinates of head surface in MRI RAS coordinate system
%  trans_ini : Initial transformation matrix from MEG to MRI RAS coordinate system
%
% - Output
%  trans_opt : Optimal transformation matrix from MEG to MRI RAS coordinate system
%
% Y. Takeda 2018-09-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

if ~exist('trans_ini','var')
    trans_ini = eye(4,3);
end

% Convert fiducials and headpoints into MRI RAS coordinate using trans_ini 
V0 = trans_coord([fids_meg; headpoint_meg], trans_ini);

Vref = [fids_ras; V_headsurface];

% Find optimal transformation matrix 
OPTIONS = optimset( ...
                    'MaxIter', 200, ...
                    'MaxFunEvals', 2*200,...
				    'TolX',    1e-6, ...
				    'TolFun',  1.e-10, ...
					'display', 'iter', ...
					'GradObj', 'off' ...
				   ); 
optparam= ...
	    fminsearch(@calc_matching_to_fids_surf, zeros(6, 1), OPTIONS, V0, Vref);

% Make output
trans_opt = vb_rigid_trans_matrix(optparam);
trans_opt = [trans_ini [zeros(3,1); 1]] * trans_opt;