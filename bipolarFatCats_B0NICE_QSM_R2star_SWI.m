%
clear all;
close all;
clc;
flag_invivo = 0; % 0 for phantom; 1 for invivo
step = 1;
if step == 0
% for Siemens
    if flag_invivo == 0; 
    load_Siemens3Tdicom_JL_new('515.mat','0001','0002');
    else
    load_Siemens3Tdicom_JL_new('603.mat','0011','0012');    
    end
end
%
for index_case = 515:515 % please check
    %
    case_str = num2str(index_case);
    if index_case < 10
    input_file = ['0',case_str,'.mat'];    
    %
    bc_file = ['0',case_str,'_bc.mat'];
    fm_file = ['0',case_str,'_fm.mat'];
    FF_file = ['0',case_str,'_FF.mat'];
    %
    QSM_file = ['0',case_str,'_QSM.mat'];
    unwrap_file = ['0',case_str,'_unwrap.mat'];
    SWI_file = ['0',case_str,'_SWI.mat'];
    R2star_file = ['0',case_str,'_R2star.mat'];
    R2star_file5echo = ['0',case_str,'_R2star_5echo.mat'];
    else
    input_file = [case_str,'.mat']; 
    %
    bc_file = [case_str,'_bc.mat'];
    fm_file = [case_str,'_fm.mat'];
    FF_file = [case_str,'_FF.mat'];
    %
    QSM_file = [case_str,'_QSM.mat'];
    unwrap_file = [case_str,'_unwrap.mat'];
    SWI_file = [case_str,'_SWI.mat'];
    R2star_file = [case_str,'_R2star.mat'];
    R2star_file5echo = [case_str,'_R2star_5echo.mat'];
    end
    %
%-------------------------------------------------------------------
    %---------------------------
    %!!!!! very important!!!!please check your lab book
    sample_type = 2; % 2 for Siemens; 1 for invivo; 0 for ex-vivo GE
    data_mode = 'bipolar';
    %---------------------------
load(input_file);
%
imDataParams.images = double(imDataParams.images);
matrix_size = size(imDataParams.images);
%
ori_imDataParams.images = imDataParams.images;
ori_imDataParams.TE = imDataParams.TE;
%
    if strcmp(data_mode,'bipolar')   
    [ imDataParams ] = bipolar_PhaseCorrection( imDataParams );
        save(bc_file, 'imDataParams');    
    end
%
if matrix_size(5) > 5
imDataParams.images(:,:,:,:,6:matrix_size(5)) = [];
imDataParams.TE(6:matrix_size(5)) = [];
end
%-------------------------------------------------------------------
algoOptions.debug = 0;
algoOptions.do_FF = 1;
%-------------------------------------------------------------------
if length(imDataParams.TE) > 2
%------------------------------
[algoParams] = FatCats_VerifyInputData(imDataParams);

%---------------------------

algoParams.debug = 2; % 0 off all the figure; 1 on all; 2 only on the final
    matrix_size = size(imDataParams.images);
    if matrix_size(4) > 1
        mag = abs(imDataParams.images).^2;
        sum_mag = sqrt(sum(mag,4));
        algoParams.mag4magfitting = sum_mag;
    else
        algoParams.mag4magfitting = abs(imDataParams.images);
    end    
        algoParams.TE4magfitting = imDataParams.TE;        
%------------------------------
%pixel-by-pixel magnitude fitting
% generate R2* maps and mag-based fat and water masks
[algoParams] = FatCats_SetAlgoParams4MagFitting(algoParams);

[algoParams] = FatCats_MagFitting(algoParams);

%------------------------------
algoParams.FiltSize = 7;
% initial B0 Phase mapping
algoParams.echo_selection = 2; % 0 auto, 1 manual, 2 fixed
if algoParams.echo_selection == 2   
    algoParams.index_B0 = [1 1 4]; % alternatively, [1 2 5];
end
[algoParams] = FatCats_initialB0PhaseMapping_SelectEcho(algoParams);
%
[algoParams] = FatCats_SetAlgoParams4PUROR(algoParams);
algoParams.upSampling = 2;
[algoParams] = FatCats_initialB0PhaseMapping_RawB0(algoParams);
%
[algoParams] = FatCats_SetAlgoParams4PhaseErrorCorrection(algoParams);
%------------------------------
% phase error correction
[algoParams] = FatCats_PhaseErrorCorrection(algoParams);
%#############################################
%bias removal
[outParams] = FatCats_BiasRemoval( algoParams );
%
%-------------------------
%
FF = outParams.FF;
FF(isnan(FF)) = 0;
fm = outParams.fm;
fm(isnan(fm)) = 0;
R2star = outParams.R2star;
R2star(isnan(R2star)) = 0;

save(FF_file, 'FF'); 
save(fm_file, 'fm');
save(R2star_file5echo, 'R2star');
%--------------------------------------------
    matrix_size = size(ori_imDataParams.images);
    flag_usingFF = 1;
    if matrix_size(5) < 6
       R2star = outParams.R2star;
    else
%-----------------------
    % re-calculating R2star with all 10 echo
    matrix_size = size(ori_imDataParams.images);
        if matrix_size(4) > 1
        mag = abs(ori_imDataParams.images).^2;
        sum_mag = sqrt(sum(mag,4));
        algoParams.mag4magfitting = sum_mag;
        else
        algoParams.mag4magfitting = abs(ori_imDataParams.images);
        end    
        algoParams.TE4magfitting = ori_imDataParams.TE;        
    %------------------------------
    %pixel-by-pixel magnitude fitting
    % generate R2* maps with all echoes
    if flag_usingFF == 1
        
    end
    [algoParams] = FatCats_MagFitting(algoParams);
    R2star = algoParams.R2star_map_raw;
    end
    save(R2star_file, 'R2star');
%
elseif length(imDataParams.TE) == 2
[ outParams ] = FatCats_main_2pt( imDataParams,algoOptions );    
else
[ outParams ] = FatCats_main_1pt( imDataParams,algoOptions );     
end
%

if strcmp(data_mode,'bipolar')
    flag_inphase = 1;
    if flag_invivo == 1
    [chi, RDF, Mask, iFreq] = test_MEDI_head_usingB0NICEfm(bc_file, fm_file, flag_inphase);
    else
    [chi, RDF, Mask, iFreq] = test_MEDI_phantom_usingB0NICEfm(bc_file, fm_file, flag_inphase);   
    end
    save(QSM_file, 'chi', 'RDF', 'Mask', 'iFreq');
    %stop
    %
    flag_HP = 1;
    [unwrapphase] = EchoPhaseUnwrapping(bc_file, fm_file,flag_HP);
    save(unwrap_file,'unwrapphase');
else
    flag_inphase = 1;
    if flag_invivo == 1
    [chi, RDF, Mask, iFreq] = test_MEDI_head_usingB0NICEfm(input_file, fm_file, flag_inphase);
    else
    [chi, RDF, Mask, iFreq] = test_MEDI_phantom_usingB0NICEfm(input_file, fm_file, flag_inphase);   
    end
    save(QSM_file, 'chi', 'RDF', 'Mask', 'iFreq');
    %
    flag_HP = 1;
    [unwrapphase] = EchoPhaseUnwrapping(input_file, fm_file,flag_HP);
    save(unwrap_file,'unwrapphase');
end
%
if flag_invivo == 1
[SWI, Phase_filt, PhaseMask] = SWI_echoBYecho(input_file, unwrap_file);
save(SWI_file,'SWI', 'Phase_filt', 'PhaseMask');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
index_case
%
end
