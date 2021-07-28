function adcData = calibAdcData(rawAdcData, calibFileName, parameter_files_path)
    adcData = struct();
    calibData = load(calibFileName);
    NumSubFramesPerFrame = length(rawAdcData);
    for i_subFrame = 0 : NumSubFramesPerFrame - 1
        adcData(i_subFrame+1).subframeIdx = rawAdcData(i_subFrame+1).subframeIdx;
        
        radar_data_Rxchain = rawAdcData(i_subFrame+1).rawAdcData;
        
        curTransferOrder = getPara(parameter_files_path(i_subFrame+1).path, 'curTransferOrder');        
        numTX = length(curTransferOrder);
        Sampling_Rate_sps = getPara(parameter_files_path(i_subFrame+1).path, 'adcSampleRate');
        chirpSlope = getPara(parameter_files_path(i_subFrame+1).path, 'chirpSlope');
        numSamplePerChirp = getPara(parameter_files_path(i_subFrame+1).path, 'numADCSample');
        nchirp_loops = getPara(parameter_files_path(i_subFrame+1).path, 'NumChirpLoops');
        
        RangeMat = calibData.calibResult.RangeMat;
        fs_calib = calibData.params.Sampling_Rate_sps;
        Slope_calib = calibData.params.Slope_MHzperus * 1e12;
        calibrationInterp = getPara(parameter_files_path(i_subFrame+1).path, 'calibrationInterp');
        PeakValMat = calibData.calibResult.PeakValMat;
        phaseCalibOnly = getPara(parameter_files_path(i_subFrame+1).path, 'calibrationCascade_phaseCalibOnly');
        
        clear outData;
        TX_ref = curTransferOrder(1);
        for iTX = 1: numTX
            TXind = curTransferOrder(iTX);
            
            %construct the frequency compensation matrix  
            freq_calib = (RangeMat(TXind,:)-RangeMat(TX_ref,1))*fs_calib/Sampling_Rate_sps *chirpSlope/Slope_calib;       
            freq_calib = 2*pi*(freq_calib)/(numSamplePerChirp * calibrationInterp);
            correction_vec = (exp(1i*((0:numSamplePerChirp-1)'*freq_calib))');
            
            freq_correction_mat = repmat(correction_vec, 1, 1, nchirp_loops);
            freq_correction_mat = permute(freq_correction_mat, [2 3 1]);
            outData1TX = radar_data_Rxchain(:,:,:,iTX).*freq_correction_mat;
            
            %construct the phase compensation matrix
            phase_calib = PeakValMat(TX_ref,1)./PeakValMat(TXind,:);
            %remove amplitude calibration
            if phaseCalibOnly == 1
                phase_calib = phase_calib./abs(phase_calib);
            end
            phase_correction_mat = repmat(phase_calib.', 1,numSamplePerChirp, nchirp_loops);
            phase_correction_mat = permute(phase_correction_mat, [2 3 1]);
            outData(:,:,:,iTX) = outData1TX.*phase_correction_mat;
            
        end
        
        adcData(i_subFrame+1).adcData = outData;
    end
end