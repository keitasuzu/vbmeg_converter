function   [vx, m] = convert_mri_fids_to_ras(original_t1, mri_fids)
% Convert mri_fids voxel coordinate to RAS (origin = center of image)
% --- Usage
%   convert_mri_fids_to_ras(fname,mri_fids);
%
% --- Input
% original_t1 : Original T1 NIFTI file name
% mri_fids : Voxel coordinate of fiducials in the original T1 NIFTI file
%
% --- Output
% vx : Voxel coordinates of fiducials in RAS (origin = center of image)
% m : Meter coordinates of fiducials in RAS (origin = center of image)
%
% Y. Takeda 2018-09-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

% Load NIFTI header & image
avw = load_nii_cbi(original_t1);

%  Check NIFTI format
if ~isfield(avw,'filetype') || avw.filetype == 0
    error('File is not Nifti-format');
end

% Orient : axis dim to get RAS coordinate
orient = get_orient_info(avw);

% Change axis of image to RAS according to orient
vx = mri_fids(:, abs(orient));

% Flip coordinate according to orient
img = permute(avw.img, abs(orient));
dim = size(img);
for d = 1:3
    if orient(d)<0
        vx(:,d) = dim(d)-vx(:,d)+1;
    end
end

% Transform voxel to meter coordinate
vsize = avw.hdr.dime.pixdim(2:4);% Voxel size of original image
R = get_rot_from_orient(orient);% Rotation matrix from orient
vsize = abs(R * vsize(:));% Voxel size of transformed RAS image
trans_x = [1 0 0 -dim(1)/2]'*vsize(1);
trans_y = [0 1 0 -dim(2)/2]'*vsize(2);
trans_z = [0 0 1 -dim(3)/2]'*vsize(3);
trans = [trans_x,trans_y,trans_z]*1e-3;
m = [vx ones(size(vx,1),1)]*trans;
