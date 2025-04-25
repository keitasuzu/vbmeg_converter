function	dmin = calc_matching_to_fids_surf(para, V0, V)
% Calculate distance for transformed points
%
% --- Input
% para   : [ dX ; th ] rigid body transformation parameter (6x1 vector)
% V0     : Target points (Npoint x 3) (1-3: fiducials, 4-: headpoints)
% V      : Reference points (Npoint x 3) (1-3: fiducials, 4-: headsurface)
%
% --- Output
% dmin = Min distance from rigid-transformed 'V0' to reference point set 'Vref'
%
% Y. Takeda 2018-09-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see POSITIONING_LICENSE.txt)

% Calculate distances between reference fiducials and transformed fiducials of target
V0_fids = V0(1:3, :);
V_fids = V(1:3, :);
dmin_fids = sum( sum(( vb_rigid_transform(V0_fids, para) - V_fids ).^2 , 2));

% Calculate distances between head surface and transformed head points
V0_headpoint = V0(4:end, :);
V_surf = V(4:end, :);

V1_headpoint = vb_rigid_transform(V0_headpoint, para);
Npoint = size(V1_headpoint, 1);
Nv_surf = size(V_surf, 1);

dmin_headpoint = 0;
for point = 1:Npoint
    d = sum((V_surf-repmat(V1_headpoint(point, :), Nv_surf, 1)).^2, 2);
    dmin_headpoint = dmin_headpoint + min(d);
end

% Calculate all the distances
dmin = dmin_fids + dmin_headpoint;



