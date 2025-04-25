function transinfo_file = make_transinfo_bids(bids_dir, subid, taskname, output_dir)
% Make transinfo_file for bids data
% 
% Y. Takeda 2018-09-21
% K. Suzuki 2025-04-26
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

sid = ['sub-' subid];

% Unzip and copy T1 image
input_file = fullfile(bids_dir, sid, 'ses-mri', 'anat', [sid '_ses-mri_acq-mprage_T1w.nii.gz']);
gunzip(input_file, output_dir);
t1_file = fullfile(output_dir, [sid '_ses-mri_acq-mprage_T1w.nii']);

% Convert fiducial coordinates to RAS coordinate (origin: ceter of the image)
json_file = fullfile(bids_dir, sid, 'ses-mri', 'anat', [sid '_ses-mri_acq-mprage_T1w.json']);
[lpa, nasion, rpa] = read_fids_from_json(json_file);
fids_orig = [lpa; nasion; rpa];
[~, fids_t1] = convert_mri_fids_to_ras(t1_file, fids_orig);
fids_file = fullfile(output_dir,'fiducial.mat');
save(fids_file, 'fids_t1')
disp([fids_file ' was saved.'])

% Convert T1 image to RAS coordinate (origin: ceter of the image)
convert_nifti_to_ras(t1_file, t1_file);
disp([t1_file ' was converted to RAS coordinate (origin is the center of the image).'])

% Import fiducial coordinates in MEG head coordinates system
json_file = fullfile(bids_dir, sid, 'ses-meg', 'meg', [sid '_ses-meg_task-' taskname '_coordsystem.json']);
[lpa, nasion, rpa] = read_fids_from_json(json_file);
fids_meg = [lpa; nasion; rpa]*1e-3;% [m]

% Import headpoint coordinates in MEG head coordinates system
pos_file = fullfile(bids_dir, sid, 'ses-meg', 'meg', [sid '_ses-meg_headshape.pos']);
headpoint_meg = read_head_pos(pos_file)*1e-3;% [m]

% Load fiducial coordinates in MRI RAS coordinate system
load(fids_file, 'fids_t1')

% Extract head surface from T1 image
[F_headsurface, V_headsurface] = vb_face_extract(t1_file);

% Make transform matrix from MEG head coordinate to MRI RAS coordinate
% using fiducials, head points, and head surface
trans_head2ras = make_trans_head2ras(fids_meg, headpoint_meg, fids_t1, V_headsurface);

% Save tranform information
transinfo_file = fullfile(output_dir, 'transinfo.mat');
save(transinfo_file, 'fids_meg', 'headpoint_meg', 'fids_t1',...
    'F_headsurface', 'V_headsurface', 'trans_head2ras');
disp([transinfo_file ' was saved.'])
