clear all; close all; clc;

% ====================user's modify=====================
radarFolder = 'K:\ourDataset\20210715\radar\mode3Group1';
calibFileName = 'input/20210715.mat';
pathGenParaFolder = 'input';
PARAM_FILE_GEN_ON = 1;

%% get radarStartTime
radarTimeFile = dir(fullfile(radarFolder, '*.startTime.txt'));
% radar start time from PC
radarStartTime = getRadarStartTime(radarTimeFile);

%% create mmwave.json parameters file
radarInfoFile = dir(fullfile(radarFolder, '*.mmwave.json'));
if PARAM_FILE_GEN_ON == 1
    parameter_files_path = parameter_file_gen_json(fullfile(radarInfoFile.folder, radarInfoFile.name), calibFileName, pathGenParaFolder);
else
    parameter_files_path = struct();
    parameter_files_path(1).path = 'input\subFrame0_param.m';
    parameter_files_path(2).path = 'input\subFrame1_param.m';
end

%% radar process
% read .bin files from radar data folder
radarBinFile_list = dir(fullfile(radarFolder, '*.bin'));
% get sliceId based path
sliceIdBasedPath = getSliceIdBasedPath(radarBinFile_list);

% create objects
% for subFrameId = 0:length(parameter_files_path)-1
%     simTopObj = simTopCascade('pfile', pathGenParaFile); %内含如何读取adc数据
%     calibrationObj = calibrationCascade('pfile', pathGenParaFile, 'calibrationfilePath', dataFolder_calib);
%     rangeFFTObj = rangeProcCascade('pfile', pathGenParaFile);%内含如何对数据进行rangeFFT处理
%     DopplerFFTObj = DopplerProcClutterRemove('pfile', pathGenParaFile);%内含如何对数据进行DopplerFFT处理
%     detectionObj = CFAR_OS('pfile', pathGenParaFile);%CFAR_OS算法
%     DOAObj = DOACascade('pfile', pathGenParaFile);
% end

% initialization
subFramesInfo = struct();
cnt_globalSubFrames = -1;
cnt_globalFrames = -1;

% get parameter from 
totNumFrames = getPara(parameter_files_path(1).path, 'frameCount');
NumSubFramesPerFrame = getPara(parameter_files_path(1).path, 'NumSubFrames');
totNumSubFrames = totNumFrames * NumSubFramesPerFrame;


for sliceId = 1:length(sliceIdBasedPath)
    num_validSubFrames = getValidFrames(sliceIdBasedPath, sliceId);
    num_validFrames = num_validSubFrames / NumSubFramesPerFrame;
    for frameId = 0 : num_validFrames-1 
        cnt_globalFrames = cnt_globalFrames +1;
        disp('===========================================================');
        fprintf('正在访问第 %s 片中第 %d/%d 帧（全局的第 %d/%d 帧）\n',  sliceIdBasedPath(sliceId).sliceId, frameId, num_validFrames-1, cnt_globalFrames, totNumFrames);
        
        for i_subFrame = 0 : NumSubFramesPerFrame -1
            cnt_globalSubFrames = cnt_globalSubFrames + 1;
            subFrameId = i_subFrame + frameId * NumSubFramesPerFrame;
            fprintf('正在读取第 %d/%d 子帧（全局的第 %d/%d 子帧）\n',  subFrameId, num_validSubFrames-1, cnt_globalSubFrames, totNumSubFrames);
            
            % record frame information
            subFramesInfo(cnt_globalSubFrames+1).globalFrameId = cnt_globalSubFrames;
            curSubFrameInfo = getFrameInfo(sliceIdBasedPath, sliceId, subFrameId);
            fieldnames_cell = fieldnames(curSubFrameInfo);
            for i_field = 1: length(fieldnames_cell)
                fieldname = fieldnames_cell{i_field};
                eval(['subFramesInfo(cnt_globalSubFrames+1).', fieldname, ' = curSubFrameInfo.', fieldname, ';']);            
            end
        end
        
        % read raw data
        rawAdcData = readAdcData(subFramesInfo, frameId, parameter_files_path);
        % calibrate raw data
        adcData = calibAdcData(rawAdcData, calibFileName, parameter_files_path);
         
         
        
        
        
        
        
        
    end
end



