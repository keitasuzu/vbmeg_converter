function trans_head2ras = make_trans_head2ras(fids_meg, headpoint_meg, fids_t1, V_headsurface)
% Make transform matrix from MEG head coordinate to MRI RAS coordinate
%
% - Input
%  fids_meg : Fiducial coordinates in MEG head coordinate system
%  headpoint_meg : Headpoint coordinates in MEG head coordinate system
%  fids_t1 : Fiducial coordinates in MRI RAS coordinate system
%  V_headsurface : Coordinates of head surface in MRI RAS coordinate system
%
% - Output
%  trans_head2ras : Transformation matrix from MEG head to MRI RAS coordinate system
%
% Y. Takeda 2018-09-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

% Make initial transformation matrix (MEG -> MRI RAS) using fids
trans_ini = vb_fit_3d_coord(fids_meg, fids_t1);

% Update the transformation matrix (MEG -> MRI RAS) using fids, head points, and head surface
trans_head2ras = fit_3d_coord2(fids_meg, headpoint_meg, fids_t1, V_headsurface, trans_ini);

