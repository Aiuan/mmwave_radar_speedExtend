clear all; close all; clc;

% ====================user's modify=====================
radarFolder = 'F:\mmwaveData\20211008\20211008_mode3_group1';

calibFileName = 'input/20210715.mat';
pathGenParaFolder = 'input';
saveFileName = 'runs/20211008_mode3_group1.mat';

%% check startTime.txt
radarTimeFile = dir(fullfile(radarFolder, '*.startTime.txt'));
% radar start time from PC
pcStartTime = getPCStartTime(radarTimeFile);

%% check mmwave.json
radarInfoFile = dir(fullfile(radarFolder, '*.mmwave.json'));
parameter_files_path = parameter_file_gen_json(fullfile(radarInfoFile.folder, radarInfoFile.name), calibFileName, pathGenParaFolder);

%% check data.bin and idx.bin
radarBinFile_list = dir(fullfile(radarFolder, '*.bin'));
% get sliceId based path
sliceIdBasedPath = getSliceIdBasedPath(radarBinFile_list);

% initialization
framesInfo = struct();
cnt_globalFrames = -1;

% get parameter from 
totNumFrames = getPara(parameter_files_path(1).path, 'frameCount') * length(parameter_files_path);

% check data and save framesInfo
for sliceId = 1:length(sliceIdBasedPath)
    % check every slice file information
    sliceIdxInfo = getSliceIdxInfo(sliceIdBasedPath, sliceId);
    num_validFrames = getValidFrames(sliceIdBasedPath, sliceId);
    for frameId = 0 : num_validFrames-1 
        cnt_globalFrames = cnt_globalFrames + 1;
        disp('===========================================================');
        fprintf('正在访问第 %s 片中第 %d/%d 帧（全局的第 %d/%d 帧）\n',  sliceIdBasedPath(sliceId).sliceId, frameId, num_validFrames-1, cnt_globalFrames, totNumFrames);
        
        framesInfo(cnt_globalFrames+1).globalFrameId = cnt_globalFrames;
        curFrameInfo = getFrameInfo(sliceIdBasedPath, sliceId, frameId);
        fieldnames_cell = fieldnames(curFrameInfo);
        for i_field = 1: length(fieldnames_cell)
            fieldname = fieldnames_cell{i_field};
            eval(['framesInfo(cnt_globalFrames+1).', fieldname, ' = curFrameInfo.', fieldname, ';']);            
        end
    end
end

%%  save check results
save(saveFileName);




