% Authors: Moreau clara Urchs Sebastian
% 2019 Montreal Pierre Bellec lab CRIUGM, Jacquemont Sebastien lan Ste Justine hospital
% More information about NIAK are available here https://niak.simexp-lab.org/build/html/CONNECTOME.html
% More information about MIST parcellation https://simexp.github.io/multiscale_dashboard/index.html

clear;

root_path = 'CNV_project/';
preproc_path = [root_path 'preproc/'];
path_out = [root_path 'glm/sample'] ;
model_path = [root_path 'pheno/pheno.csv'];

files_in.networks.cambridge64 = [root_path 'template/template_cambridge_basc_multiscale_sym_scale064.nii.gz'];
files_in.networks.cambridge12 = [root_path 'template/template_cambridge_basc_multiscale_sym_scale012.nii.gz'];

opt_g.exclude_subject = {'s14725xx46xFCAP1','s14785xx5xFCAP1', 's14871xx1xFCAP1', 's14927xx1xFCAP1', 's14928xx1xFCAP1', 's14880xx1xFCAP1', 's14952xx5xFCAP1'};
opt_g.min_nb_vol = 40; % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.type_files = 'glm_connectome'; % Specify to the grabber to prepare the files for the glm_connectome pipeline
tmp =  niak_grab_fmri_preprocess(preproc_path, opt_g);
files_in.fmri = tmp.fmri;
files_in.model.group = model_path;

opt.fdr = 0.05;
opt.folder_out = path_out; % Where to store the results

%%%%%%%%%%%%
%Tests
%%%%%%%%%%%%

%%%% mean con  %%%%%
opt.test.mean_con.group.contrast.g1 = 0; %del
opt.test.mean_con.group.contrast.g2 = 1; %con
opt.test.mean_con.group.contrast.g3 = 0; %dup
opt.test.mean_con.group.contrast.sex_dummy = 0;
opt.test.mean_con.group.contrast.FD_scrubbed_both_norm = 0;
opt.test.mean_con.group.contrast.Site_dummy = 0;
opt.test.mean_con.group.flag_intercept = false;
opt.test.mean_con.group.normalize_x = false;
opt.test.mean_con.group.normalize_y = false;

opt.flag_test = false; 
[pipeline,opt] = niak_pipeline_glm_connectome(files_in,opt);
