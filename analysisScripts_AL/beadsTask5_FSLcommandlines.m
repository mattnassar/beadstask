subjects = {'CC0528','FL1136','FT3594','IL3520','LL3555',...
            'NS0950','OF3592','TK3556','TQ3543','TQ3600'};
        
betsofmags = [];
erodemasks = [];
prepfield = [];
medfilt = [];

% commands for B0maps
for i = 1:length(subjects)
    b0mag = [subjects{i},'_B0map_mag.nii.gz'];
    b0magMask = [subjects{i},'_B0map_Magmask.nii.gz'];
    b0Mask = [subjects{i},'_B0map_Magmask_mask.nii.gz'];
    b0mag_brain = [subjects{i},'_B0map_mag_brain.nii.gz'];
    b0phase = [subjects{i},'_B0map_phase.nii.gz'];
    b0prep = [subjects{i},'_B0map_prep.nii.gz'];
    b0prep_med = [subjects{i},'_B0map_prep_med.nii.gz'];
    
    betsofmags = [betsofmags, 'bet ', b0mag,' ', b0magMask, ' -n -m & '];
    erodemasks = [erodemasks, 'fslmaths ', b0Mask, ' -kernel boxv 5 -ero -mul ', b0mag, ' ', b0mag_brain, ' & '];
    prepfield = [prepfield, 'fsl_prepare_fieldmap SIEMENS ', b0phase, ' ', b0mag_brain, ' ', b0prep, ' 2.46 & '];
    medfilt = [medfilt, 'fugue --loadfmap=',b0prep,' --savefmap=',b0prep_med,' --mask=',b0Mask,' --median & '];

end