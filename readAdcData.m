%% read raw adc data with MIMO 

function rawAdcData = readAdcData(subFramesInfo, globalSubFrameId_start, parameter_files_path)
    rawAdcData = struct();
    NumSubFramesPerFrame = length(parameter_files_path);
    for i_subFrame = 0 : NumSubFramesPerFrame - 1
        tic;
        subframeIdx = i_subFrame + globalSubFrameId_start;
        % 20210903读取bug修正subframeIdx_inSlice
        subframeIdx_inSlice = subFramesInfo(subframeIdx + 1).sliceFrameId;
        fprintf('>>正在读取第 %d/%d 子帧       ',  subframeIdx, length(subFramesInfo) - 1);
        fileFullPath_master = subFramesInfo(subframeIdx + 1).master_adcDataPath;
        fileFullPath_slave1 = subFramesInfo(subframeIdx + 1).slave1_adcDataPath;
        fileFullPath_slave2 = subFramesInfo(subframeIdx + 1).slave2_adcDataPath;
        fileFullPath_slave3 = subFramesInfo(subframeIdx + 1).slave3_adcDataPath; 
        numSamplePerChirp = getPara(parameter_files_path(i_subFrame+1).path, 'numADCSample');
        numChirpPerLoop = getPara(parameter_files_path(i_subFrame+1).path, 'NumChirp');
        numLoops = getPara(parameter_files_path(i_subFrame+1).path, 'NumChirpLoops');
        numRX = getPara(parameter_files_path(i_subFrame+1).path, 'numRxToEnable');
        numDevices = getPara(parameter_files_path(i_subFrame+1).path, 'NumDevices');
        numRXPerDevice = numRX / numDevices;

        [radar_data_Rxchain_master] = readBinFile(fileFullPath_master, subframeIdx_inSlice,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);
        [radar_data_Rxchain_slave1] = readBinFile(fileFullPath_slave1, subframeIdx_inSlice,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);
        [radar_data_Rxchain_slave2] = readBinFile(fileFullPath_slave2, subframeIdx_inSlice,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);
        [radar_data_Rxchain_slave3] = readBinFile(fileFullPath_slave3, subframeIdx_inSlice,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);

        % Arranged based on Master RxChannels, Slave1 RxChannels, slave2 RxChannels, slave3 RxChannels
        clear radar_data_Rxchain;
        radar_data_Rxchain(:,:,1:4,:) = radar_data_Rxchain_master;
        radar_data_Rxchain(:,:,5:8,:) = radar_data_Rxchain_slave1;
        radar_data_Rxchain(:,:,9:12,:) = radar_data_Rxchain_slave2;
        radar_data_Rxchain(:,:,13:16,:) = radar_data_Rxchain_slave3;
        
        rawAdcData(i_subFrame+1).subframeIdx = subframeIdx;
        rawAdcData(i_subFrame+1).rawAdcData = radar_data_Rxchain;
        
        fprintf('耗时%.3f s\n', toc);
    end
        
end


function [adcData1Complex] = readBinFile(fileFullPath, frameIdx_inSlice,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice)
    Expected_Num_SamplesPerFrame = numSamplePerChirp*numChirpPerLoop*numLoops*numRXPerDevice*2; % 这里的2表示IQ两路复数信号
    fp = fopen(fileFullPath, 'r');
    % 20210903读取bug修正frameIdx_inSlice，不应为（frameIdx_inSlice-1）
    fseek(fp,frameIdx_inSlice*Expected_Num_SamplesPerFrame*2, 'bof');% 这里的2表示信号以16位格式存储，占2字节
    adcData1 = fread(fp,Expected_Num_SamplesPerFrame,'uint16');
    neg             = logical(bitget(adcData1, 16));
    adcData1(neg)    = adcData1(neg) - 2^16;
    %% 
    adcData1 = adcData1(1:2:end) + sqrt(-1)*adcData1(2:2:end);
    adcData1Complex = reshape(adcData1, numRXPerDevice, numSamplePerChirp, numChirpPerLoop, numLoops);
    adcData1Complex = permute(adcData1Complex, [2 4 1 3]);
    fclose(fp);
end