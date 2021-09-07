function dopplerFFTOut = dopplerFFT(rangeFFTOut, parameter_files_path, dopplerFFTFilter_ON)
    dopplerFFTSize = getPara(parameter_files_path, 'DopplerProcClutterRemove_dopplerFFTSize');
    dopplerFFT = zeros(size(rangeFFTOut, 1), dopplerFFTSize, size(rangeFFTOut, 3), size(rangeFFTOut, 4));
    
     for i_tx = 1 : size(rangeFFTOut, 4)
        for i_rx = 1 : size(rangeFFTOut, 3)
            % vectorized version
            inputMat = rangeFFTOut(:, :, i_rx, i_tx);
            
            if dopplerFFTFilter_ON
                % apply doppler-domain windowing
                dopplerWindowCoeff = getPara(parameter_files_path, 'DopplerProcClutterRemove_dopplerWindowCoeff');
                numChirpsPerVirAnt = getPara(parameter_files_path, 'DopplerProcClutterRemove_numChirpsPerVirAnt');
                dopplerWindowCoeffVec = ones(numChirpsPerVirAnt, 1);
                dopplerWindowCoeffVec(1 : length(dopplerWindowCoeff)) = dopplerWindowCoeff;
                dopplerWindowCoeffVec(numChirpsPerVirAnt-length(dopplerWindowCoeff)+1:numChirpsPerVirAnt) = dopplerWindowCoeffVec(length(dopplerWindowCoeff):-1:1);
                inputMat    = bsxfun(@times, inputMat, dopplerWindowCoeffVec.');

                clutterRemove = getPara(parameter_files_path, 'DopplerProcClutterRemove_clutterRemove');
                if clutterRemove ==1
                    inputMat = inputMat - (repmat(mean(inputMat'),size(inputMat,2),1))';
                end
            end
            
            % Doppler FFT
            fftOutput   = fft(inputMat, dopplerFFTSize, 2);                         
            
            if dopplerFFTFilter_ON
                FFTOutScaleOn = getPara(parameter_files_path, 'DopplerProcClutterRemove_FFTOutScaleOn');
                if FFTOutScaleOn ==1
                    scaleFactorDoppler = getPara(parameter_files_path, 'DopplerProcClutterRemove_scaleFactorDoppler');
                    fftOutput   = fftshift(fftOutput, 2) * scaleFactorDoppler;
                else
                    fftOutput   = fftshift(fftOutput, 2);
                end
            end
            
            % populate in the data cube
            dopplerFFTOut(:, :, i_rx, i_tx) = fftOutput;
            
        end
     end
  
end