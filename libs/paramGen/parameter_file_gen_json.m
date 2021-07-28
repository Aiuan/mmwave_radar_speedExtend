% parameter_file_gen_Jason.m
%  
% Function call to generate a complete parameter file using Jason file before signal processing, including chirp parameters and module parameters.
% This function does not support advanced frame config
%input:
%   test_name: test data file name which contains the chirp configuration
%   file in a Jason file format
%   dataFolder_calib: calibration data file path which contains the pre-calibrated .mat file. This .mat file includes the calibration matrix to be used in adc data calibration module
%   module_param_file: contains parameters to initialize each signal processing module
%   pathGenParaFile: the complete parameter file name to be generated
%   dataPlatform: flag for data capture platform


function param_results = parameter_file_gen_json(paramFile, dataFolder_calib, pathGenParaFolder)
    % 预设信息
    module_param_file = 'module_param.m';
    dataPlatform = 'TDA2';
 
    % 解析mmwave.json文件
    params_chirp = JsonParser(paramFile);
    
    numChirpConfig = length(params_chirp.DevConfig(1).Chirp);%device1中设置的chirps数量
    numTXPerDev = 3;%每个芯片发射通道为3
    totTx = numTXPerDev*params_chirp.NumDevices;%3*4=12，一共12发射通道
    TxEnableTable = zeros(numChirpConfig, totTx);%构建每个device的每个chirp的占用Tx情况表，维度本应为 device数量*chirp数量*12发射通道，实际合并了device维度
    for iDev = 1:params_chirp.NumDevices%第iDev个设备
        for iconfig = 1:numChirpConfig%第iconfig个chirp
            TxEnableTable(iconfig,1+(iDev-1)*numTXPerDev) = ...
                params_chirp.DevConfig(iDev).Chirp(iconfig).Tx0Enable;
            TxEnableTable(iconfig, 2+(iDev-1)*numTXPerDev) = ...
                params_chirp.DevConfig(iDev).Chirp(iconfig).Tx1Enable;
            TxEnableTable(iconfig, 3+(iDev-1)*numTXPerDev) = ...
                params_chirp.DevConfig(iDev).Chirp(iconfig).Tx2Enable;
        end
    end    
    TxChannelEnabled = zeros(1,numChirpConfig);%构建按顺序发射的chirp所使用的TX的索引矩阵，维度为 1*chirp数量
    for iconfig = 1:numChirpConfig%第iconfig个chirp
        [channelID] = find(TxEnableTable(iconfig,:)~=0);%寻找第iconfig个chirp使能的天线
         TxChannelEnabled(iconfig) = channelID;     
    end
    
    param_results = struct();
    for subFrameId = 0:params_chirp.DevConfig(1).AdvFrame.NumSubFrames-1
        pathGenParaFile = fullfile(pathGenParaFolder, ['subFrame', num2str(subFrameId), '_param.m']);
        param_results(subFrameId+1).path = pathGenParaFile;
        
        % open or create param file
        fidParam = fopen(pathGenParaFile, 'w');            
        fprintf(fidParam, 'subFrameId = %d; \n', subFrameId);
        fprintf(fidParam, 'dataPlatform = ''%s''; \n', dataPlatform);
        
        % 将mmwave.json中的配置信息全部写入将要新生成的参数文件
        fprintf(fidParam, '%%%% pass the devices parameters associated with mmwave.json \n');
        fprintf(fidParam, 'NumDevices = %d; \n', params_chirp.NumDevices);
        fprintf(fidParam, '\n');
        % chirp parameters
        fprintf(fidParam, '%%master chirp parameters: \n');
        fprintf(fidParam, 'numADCSample = %e; \n', params_chirp.DevConfig(1).Profile.NumSamples);
        fprintf(fidParam, 'adcSampleRate = %e; %%Hz/s \n', params_chirp.DevConfig(1).Profile.SamplingRate*1e3);
        fprintf(fidParam, 'startFreqConst = %e; %%Hz \n', params_chirp.DevConfig(1).Profile.StartFreq*1e9);
        fprintf(fidParam, 'chirpSlope = %e; %%Hz/s \n',params_chirp.DevConfig(1).Profile.FreqSlope*1e12);

        num_chirpsInLoop = params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).NumChirp;
        if params_chirp.DevConfig(1).Chirp(subFrameId * num_chirpsInLoop + 1).IdleTimeVar == 0
            fprintf(fidParam, 'chirpIdleTime = %e; %%s \n', params_chirp.DevConfig(1).Profile.IdleTime*1e-6);
        else
            fprintf(fidParam, 'chirpIdleTime = %e; %%s \n', params_chirp.DevConfig(1).Chirp(subFrameId * num_chirpsInLoop + 1).IdleTimeVar *1e-6);
        end
        
        fprintf(fidParam, 'adcStartTimeConst = %e; %%s \n', params_chirp.DevConfig(1).Profile.AdcStartTime*1e-6);
        fprintf(fidParam, 'chirpRampEndTime = %d; %%s \n', params_chirp.DevConfig(1).Profile.RampEndTime*1e-6);
        fprintf(fidParam, '\n');

        % 'advancedFrameChirp'
        fprintf(fidParam, '%%master advancedFrameChirp parameters: \n');
        fprintf(fidParam, 'frameCount = %d; \n', params_chirp.DevConfig(1).AdvFrame.NumFrames);
        fprintf(fidParam, 'NumSubFrames = %d; \n', params_chirp.DevConfig(1).AdvFrame.NumSubFrames);
        fprintf(fidParam, '%%current subFrame advancedFrameChirp parameters: \n');
        fprintf(fidParam, 'ChirpStartIdx = %d; \n', params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).ChirpStartIdx);
        fprintf(fidParam, 'NumChirp = %d; \n', params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).NumChirp);
        fprintf(fidParam, 'NumChirpLoops = %d; \n', params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).NumChirpLoops);
        fprintf(fidParam, 'BurstPeriod = %d; %%ms \n', params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).BurstPeriod);
        fprintf(fidParam, 'NumBurst = %d; \n', params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).NumBurst);
        fprintf(fidParam, 'NumBurstLoops = %d; \n', params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).NumBurstLoops);
        fprintf(fidParam, 'SubFramePeriod = %d; %%ms \n', params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).SubFramePeriod);
        
        fprintf(fidParam, '\n');    

        fprintf(fidParam, 'numTxToEnable = %d; \n', length(params_chirp.TxToEnable));
        fprintf(fidParam, ['TxToEnable = [' num2str(params_chirp.TxToEnable) ']' ';\n']);
        fprintf(fidParam, 'numRxToEnable = %d; \n', length(params_chirp.RxToEnable));
        fprintf(fidParam, ['RxToEnable = [' num2str(params_chirp.RxToEnable) ']' ';\n']);
        
        totTransferOrder = [];
        for i = 0:params_chirp.DevConfig(1).AdvFrame.NumSubFrames-1
            temp_chirpStartIdx = params_chirp.DevConfig(1).AdvFrame.SubFrame(i+1).ChirpStartIdx;
            temp_num_chirpsInLoop = params_chirp.DevConfig(1).AdvFrame.SubFrame(i+1).NumChirp;
            totTransferOrder = [totTransferOrder TxChannelEnabled(temp_chirpStartIdx + 1 : temp_chirpStartIdx + temp_num_chirpsInLoop)];
        end
        fprintf(fidParam, ['totTransferOrder = [' num2str(totTransferOrder) ']' ';\n']);
        
        chirpStartIdx = params_chirp.DevConfig(1).AdvFrame.SubFrame(subFrameId+1).ChirpStartIdx;
        curTransferOrder = TxChannelEnabled(chirpStartIdx + 1 : chirpStartIdx + num_chirpsInLoop);
        fprintf(fidParam, ['curTransferOrder = [' num2str(curTransferOrder) ']' ';\n']);

        Start_Freq_Hz = params_chirp.DevConfig(1).Profile.StartFreq * 1e9;
        AdcStartTime_s = params_chirp.DevConfig(1).Profile.AdcStartTime * 1e-6;
        Samples_per_Chirp = params_chirp.DevConfig(1).Profile.NumSamples;
        Sampling_Rate_sps = params_chirp.DevConfig(1).Profile.SamplingRate * 1e3;
        Slope_Hzpers = params_chirp.DevConfig(1).Profile.FreqSlope * 1e12;
        CenterFreq = Start_Freq_Hz + (AdcStartTime_s + Samples_per_Chirp / Sampling_Rate_sps / 2) * Slope_Hzpers;
        fprintf(fidParam, 'centerFreq = %d; \n', (CenterFreq / 1e9));
        fprintf(fidParam, '\n');   


        %% 加载校准相关参数
        paramsCalib = load(dataFolder_calib);
        fprintf(fidParam, '%%%% pass the slope used for calibration \n');
        fprintf(fidParam, 'Slope_calib = %d; %%Hz \n', paramsCalib.params.Slope_MHzperus*1e12);
        fprintf(fidParam, 'fs_calib = %d; %%sps \n\n', paramsCalib.params.Sampling_Rate_sps);

        %% 利用module_param_file载入额外的配置信息，包括rangeFFT、DopplerFFT、CFAR、DOA过程中的超参数
        fprintf(fidParam, '%%%% pass all other parameters \n');
        fidCommon = fopen(module_param_file);
        %跳过module_param_file中前32行注释
        for ii = 1:32
            tline = fgets(fidCommon);
        end
        %将module_param_file后续部分写入pathGenParaFile
        tline = fgets(fidCommon);
        while ischar(tline)
            fwrite(fidParam, tline);
            tline = fgets(fidCommon);
        end    
        fclose(fidCommon);


        fclose(fidParam);
        
    end
    
    
    
    
end


