#!/bin/bash

#Written by Noor Z. Al Dahhan

#Running script for each subject in folder
for subject in ####
do 
    SUBJ=$subject
	cd $subject
	echo $subject

function []= Cortex_Dyslexia_ADHD_EF()

%onsets
%instr	[0,42,64,86,128,170,212,234,256,278,320,362,404,426,448,490]			
%fix	[22,108,150,192,300,342,384,470]			
%rhyme	[2, 66, 130, 214,  258, 322, 428, 492]			
%faces	[44, 88, 172,236,  280, 364, 406,  450]			


%Set the paths for main directory and spm12
maindir='___';
addpath('...spm12/');

%Specify subject's functional data:

FUNCTIONAL_FILE={''};


%Specify data output directory:
datadir = '';


%STEPS TO BE RUN. 1 = RUN, 0=DON'T RUN
firstlevel = 1;
firstlevel_cons = 1;


%% L1
if firstlevel == 1
mask = ',1'; %define whole brain mask

%Phono Task
for i=1:length(FUNCTIONAL_FILE)

    curSubj = FUNCTIONAL_FILE{i};
    fprintf('\nWorking on %s...\n',curSubj); %tells you what participant is being run
    subdatadir = strcat(datadir,curSubj,'/phono/'); %PATH TO THE SMOOTHED DATA

    rpfile = strcat(subdatadir,curSubj,'_task-phono_fd.txt'); %PATH TO OUTLIERS (Actual file)
        if exist(rpfile)
        rpfile = rpfile;
        else 
        rpfile = '';
        end
    outdir = strcat(datadir,curSubj,'/phono/results/'); %PATH TO WHERE SPM.MAT WILL BE SAVED    %

    spm_jobman('initcfg')
    spm('defaults', 'fmri');

    clear matlabbatch

matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {subdatadir}; %change to task-specific directory for smoothed files
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^s_'; %change filter to s_ - the smoothed file.
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{2}.spm.util.exp_frames.files(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^s_)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.util.exp_frames.frames = Inf; %Expand everything
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {subdatadir}; %is ok.
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^txt';
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{4}.spm.stats.fmri_spec.dir = {outdir}; %change to task-specific output directory
matlabbatch{4}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{4}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t = 16; %was 32 before;verified (number of slices per volume)
matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t0 = 8; %was 1 before; Acquisition is interleaved ascending. For even # it is 2, 4, 6, ..., 1, 3, 5. middle is 1 then?
matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(1).name = 'rhyme'; %Cut/paste from timing file
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(1).onset = [2
                                                         66
                                                         130
                                                         214
                                                         258
                                                         322
                                                         428
                                                         492];

matlabbatch{4}.spm.stats.fmri_spec.sess.cond(1).duration = [20
                                                            20
                                                            20
                                                            20
                                                            20
                                                            20
                                                            20
                                                            20]; %duration of block(s)
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(2).name = 'faces'; %same as above
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(2).onset = [44
                                                         88
                                                         172
                                                         236
                                                         280
                                                         364
                                                         406
                                                         450];
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(2).duration = [20
                                                            20
                                                            20
                                                            20
                                                            20
                                                            20
                                                            20
                                                            20];
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{4}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
matlabbatch{4}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{4}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg = {rpfile}; %outlier file as regressor (1s and 0s file) %for no outliers, put ''. if outliers, put rpfile with no quotes
matlabbatch{4}.spm.stats.fmri_spec.sess.hpf = 128; %leave this and all below
matlabbatch{4}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{4}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{4}.spm.stats.fmri_spec.volt = 1;
matlabbatch{4}.spm.stats.fmri_spec.global = 'None';
matlabbatch{4}.spm.stats.fmri_spec.mthresh = -Inf;
matlabbatch{4}.spm.stats.fmri_spec.mask = {mask};
matlabbatch{4}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{5}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{5}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch);

clear matlabbatch
clear subdatadir
clear outdir

end 
end

%% L1 Contrasts
if firstlevel_cons == 1

for i=1:length(FUNCTIONAL_FILE)
   
    cursubj = FUNCTIONAL_FILE{i};
    datadir = ''; %root directory
    outdir = strcat(datadir,cursubj,'/phono/results/');
    phono_spmfile = strcat(outdir,'SPM.mat');
    spm_jobman('initcfg')
    spm('defaults', 'fmri');

%-----------------------------------------------------------------------
% PHONO CONTRASTS & T-TESTS
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.con.spmmat = {phono_spmfile};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'WordRhyme';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'FaceMatch';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Words>Faces';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run', matlabbatch);
clear matlabbatch


end


%%T-tests for each group

spm_jobman('initcfg')
spm('defaults', 'fmri');


%% For each group:
% Phono - Words
outdir = '';
matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {
'AllTDIndividuals'};

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'mask,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);

% Phono - Faces
outdir = '';
matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {
'AllTDIndividuals'};
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'mask,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);


% Phono - Words v Faces
outdir = '';
matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {
'AllTDIndividuals'};
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'mask,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);

%



