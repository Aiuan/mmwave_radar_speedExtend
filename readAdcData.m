%% read raw adc data with MIMO 

function rawAdcData = readAdcData(subFramesInfo, frameId, parameter_files_path)
    rawAdcData = struct();
    NumSubFramesPerFrame = length(parameter_files_path);
    for i_subFrame = 0 : NumSubFramesPerFrame - 1
        subframeIdx = i_subFrame + frameId * NumSubFramesPerFrame;        
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

        [radar_data_Rxchain_master] = readBinFile(fileFullPath_master, subframeIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);
        [radar_data_Rxchain_slave1] = readBinFile(fileFullPath_slave1, subframeIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);
        [radar_data_Rxchain_slave2] = readBinFile(fileFullPath_slave2, subframeIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);
        [radar_data_Rxchain_slave3] = readBinFile(fileFullPath_slave3, subframeIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice);

        % Arranged based on Master RxChannels, Slave1 RxChannels, slave2 RxChannels, slave3 RxChannels
        clear radar_data_Rxchain;
        radar_data_Rxchain(:,:,1:4,:) = radar_data_Rxchain_master;
        radar_data_Rxchain(:,:,5:8,:) = radar_data_Rxchain_slave1;
        radar_data_Rxchain(:,:,9:12,:) = radar_data_Rxchain_slave2;
        radar_data_Rxchain(:,:,13:16,:) = radar_data_Rxchain_slave3;
        
        rawAdcData(i_subFrame+1).subframeIdx = subframeIdx;
        rawAdcData(i_subFrame+1).rawAdcData = radar_data_Rxchain;
    end
        
end


function [adcData1Complex] = readBinFile(fileFullPath, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice)
    Expected_Num_SamplesPerFrame = numSamplePerChirp*numChirpPerLoop*numLoops*numRXPerDevice*2;
    fp = fopen(fileFullPath, 'r');
    fseek(fp,(frameIdx-1)*Expected_Num_SamplesPerFrame*2, 'bof');
    adcData1 = fread(fp,Expected_Num_SamplesPerFrame,'uint16');
    neg             = logical(bitget(adcData1, 16));
    adcData1(neg)    = adcData1(neg) - 2^16;
    %% 
    adcData1 = adcData1(1:2:end) + sqrt(-1)*adcData1(2:2:end);
    adcData1Complex = reshape(adcData1, numRXPerDevice, numSamplePerChirp, numChirpPerLoop, numLoops);
    adcData1Complex = permute(adcData1Complex, [2 4 1 3]);
    fclose(fp);
end