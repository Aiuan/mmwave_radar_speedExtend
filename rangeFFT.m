function rangeFFTOut = rangeFFT(adcData, parameter_files_path, rangeFFTFilter_ON)
    rangeFFTSize = getPara(parameter_files_path, 'rangeProcCascade_rangeFFTSize');
    rangeFFTOut = zeros(rangeFFTSize, size(adcData, 2), size(adcData, 3), size(adcData, 4));
    
    for i_tx = 1 : size(adcData, 4)
        for i_rx = 1 : size(adcData, 3)
            % vectorized version
            inputMat = adcData(:, :, i_rx, i_tx);
            
            if rangeFFTFilter_ON
                % DC offset compensation
                inputMat = bsxfun(@minus, inputMat, mean(inputMat));
                % apply range-domain windowing
                rangeWindowCoeff = getPara(parameter_files_path, 'rangeProcCascade_rangeWindowCoeff');
                numAdcSamplePerChirp = getPara(parameter_files_path, 'rangeProcCascade_numAdcSamplePerChirp');
                rangeWindowCoeffVec = ones(numAdcSamplePerChirp, 1);
                rangeWindowCoeffVec(1 : length(rangeWindowCoeff)) = rangeWindowCoeff;
                rangeWindowCoeffVec(numAdcSamplePerChirp-length(rangeWindowCoeff)+1 : numAdcSamplePerChirp) = rangeWindowCoeffVec(length(rangeWindowCoeff):-1:1);
                inputMat = bsxfun(@times, inputMat, rangeWindowCoeffVec);
            end
            
            % Range FFT
            fftOutput   = fft(inputMat, rangeFFTSize);
            
            if rangeFFTFilter_ON
                FFTOutScaleOn = getPara(parameter_files_path, 'rangeProcCascade_FFTOutScaleOn');
                if  FFTOutScaleOn == 1
                    scaleFactorRange = getPara(parameter_files_path, 'rangeProcCascade_scaleFactorRange');
                    fftOutput   = fftOutput * scaleFactorRange;
                end
            end
            
            % populate in the data cube
            rangeFFTOut(:, :, i_rx, i_tx) = fftOutput;
        end
    end

end