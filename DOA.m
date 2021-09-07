function out = DOA(detected_obj, parameter_file_path)
    method = getPara(parameter_file_path, 'DOACascade_method');
    
    numObj = length(detected_obj);
    out = detected_obj;
    numAoAObjCnt = 0;
    
    % extended detection_obj to include the angles information
    for i_obj = 1:numObj
        current_obj = detected_obj(i_obj);
%         estSNR = 10*log10(sum(abs(current_obj.bin_val).^2)/sum(current_obj.noise_var));
        X = current_obj.bin_val; 
%         R = X*X';
        switch method
            case 1
                %2D beamforming angle estimation, azimuth is estimated based on 1D FFT output                        
                [DOA_angles, angle_sepc_2D_fft ]= DOA_beamformingFFT_2D(parameter_file_path, X);
                if (numAoAObjCnt == 0)
                    out = [];
                end
                for ii_obj = 1 : size(DOA_angles, 2)
                    numAoAObjCnt = numAoAObjCnt+1;
                    out(numAoAObjCnt).rangeInd = current_obj.rangeInd;
                    out(numAoAObjCnt).dopplerInd = current_obj.dopplerInd;
                    out(numAoAObjCnt).range = current_obj.range;
                    out(numAoAObjCnt).doppler_corr = current_obj.doppler_corr;
                    out(numAoAObjCnt).dopplerInd_org = current_obj.dopplerInd_org;

                    out(numAoAObjCnt).noise_var = current_obj.noise_var;
                    out(numAoAObjCnt).bin_val = current_obj.bin_val;
                    out(numAoAObjCnt).estSNR = current_obj.estSNR;
                    out(numAoAObjCnt).doppler_corr_overlap = current_obj.doppler_corr_overlap;
                    out(numAoAObjCnt).doppler_corr_FFT = current_obj.doppler_corr_FFT;

                    out(numAoAObjCnt).angles = DOA_angles(:, ii_obj);
                    out(numAoAObjCnt).spectrum = angle_sepc_2D_fft;
                end
                
            case 2
                %2D beamforming, angle estimated after 2D FFT jointly
                [DOA_angles, angle_sepc_2D_fft]= DOA_beamformingFFT_2D_joint(parameter_file_path, X);
                if (numAoAObjCnt == 0)
                    out = [];
                end
                for ii_obj = 1:size(DOA_angles,2)
                    numAoAObjCnt = numAoAObjCnt+1;
                    out(numAoAObjCnt).rangeInd = current_obj.rangeInd;
                    out(numAoAObjCnt).dopplerInd = current_obj.dopplerInd;
                    out(numAoAObjCnt).range = current_obj.range;
                    out(numAoAObjCnt).doppler_corr = current_obj.doppler_corr;
                    out(numAoAObjCnt).dopplerInd_org = current_obj.dopplerInd_org;

                    out(numAoAObjCnt).noise_var = current_obj.noise_var;
                    out(numAoAObjCnt).bin_val = current_obj.bin_val;
                    out(numAoAObjCnt).estSNR = current_obj.estSNR;
                    out(numAoAObjCnt).doppler_corr_overlap = current_obj.doppler_corr_overlap;
                    out(numAoAObjCnt).doppler_corr_FFT = current_obj.doppler_corr_FFT;

                    out(numAoAObjCnt).angles = DOA_angles(:,ii_obj);
                    out(numAoAObjCnt).spectrum = angle_sepc_2D_fft;
                end

            otherwise
                error('Not specified DOA method')

        end
    end
end



%DOA_beamformingFFT_2D.m
% DOA_beamformingFFT_2D function perform 2D angle estimation based on FFT beamforming, the azimuth peak selection
%is done in 1D FFT domain, the elevation peak selection is done after 2D FFT
%input:
%   obj: object instance
%   sig: complex signal vector, with each value corresponding to each
%   antenna. The length of this vector equals to numTX x numRX enabled.
%   There can be overlapped antennas. this signal needs to be re-arranged
%   based on D value to form the virtual antenna array
%output:
%   angleObj_est: angle estimation results
%   angle_sepc_2D_fft: angle 2D fft spectrum
function [angleObj_est, angle_sepc_2D_fft]= DOA_beamformingFFT_2D(parameter_file_path, sig)
    sidelobeLevel_dB_azim = getPara(parameter_file_path, 'DOACascade_sidelobeLevel_dB_azim');
    sidelobeLevel_dB_elev = getPara(parameter_file_path, 'DOACascade_sidelobeLevel_dB_elev');

    %field of view to do beamforming
    angles_DOA_az = getPara(parameter_file_path, 'DOACascade_angles_DOA_az');
    angles_DOA_ele = getPara(parameter_file_path, 'DOACascade_angles_DOA_ele');

    %distance unit in terms of wavelength
    d = getPara(parameter_file_path, 'DOACascade_antDis');
    %2D matrix providing antenna coordinates
    D = getPara(parameter_file_path, 'DOACascade_D');
    angleFFTSize = getPara(parameter_file_path, 'DOACascade_DOAFFTSize');
    

    %FFT based implementation
    %first form a 2D matrix based on the antenna coordinates
    D = D + 1;
    apertureLen_azim = max(D(:,1));
    apertureLen_elev = max(D(:,2));
    sig_2D = zeros(apertureLen_azim,apertureLen_elev);
    for i_line = 1:apertureLen_elev
        ind = find(D(:,2) == i_line);
        D_sel = D(ind,1);%azi_position, when ele_position == i_line
        sig_sel = sig(ind);%value, when ele_position == i_line
        [~, indU] = unique(D_sel);%val is nonredundant azi_position; indU is index
        sig_2D(D_sel(indU),i_line) = sig_sel(indU);
    end

    %run FFT on azimuth and elevation
    angle_sepc_1D_fft=fftshift(fft(sig_2D,angleFFTSize,1),1); 
    angle_sepc_2D_fft=fftshift(fft(angle_sepc_1D_fft,angleFFTSize,2),2); 

    wx_vec=-pi : 2*pi/angleFFTSize : pi;
    wz_vec=-pi : 2*pi/angleFFTSize : pi;
    wx_vec = wx_vec(1:end-1);
    wz_vec = wz_vec(1:end-1);
    %use one row with complete azimuth antenna of 1D FFT output for azimuth
    %estimation
    spec_azim = abs(angle_sepc_1D_fft(:,1));
    [~, peakLoc_azim] = DOA_BF_PeakDet_loc(parameter_file_path, spec_azim, sidelobeLevel_dB_azim);

    if apertureLen_elev ==1
        %azimuth array only, no elevation antennas
        obj_cnt = 1;
        angleObj_est= [];
        for i_obj = 1:length(peakLoc_azim)
            ind = peakLoc_azim(i_obj);
            azim_est = asind(wx_vec(ind)/(2*pi*d));
            if (azim_est >= angles_DOA_az(1) && azim_est <= angles_DOA_az(2))
                angleObj_est(1,obj_cnt) = azim_est;
                angleObj_est(2,obj_cnt) = 0;
                angleObj_est(3,obj_cnt) = ind;
                angleObj_est(4,obj_cnt) = 0;
                obj_cnt = obj_cnt+1;
            else
                continue;
            end
        end
    else
        %azimuth and elevation angle estimation

        % figure(1);plot(spec_azim); hold on; grid on
        % plot(peakLoc_azim, spec_azim(peakLoc_azim),'ro');hold on

        %for each detected azimuth, estimate its elevation
        % figure(2)
        obj_cnt = 1;
        angleObj_est= [];
        for i_obj = 1:length(peakLoc_azim)
            ind = peakLoc_azim(i_obj);
            spec_elev = abs(angle_sepc_2D_fft(ind,:));
            [peakVal_elev, peakLoc_elev] = DOA_BF_PeakDet_loc(parameter_file_path, spec_elev, sidelobeLevel_dB_elev);
            %calcualte the angle values
            for j_elev = 1:length(peakVal_elev)
                azim_est = asind(wx_vec(ind)/(2*pi*d));
                elev_est = asind(wz_vec(peakLoc_elev(j_elev))/(2*pi*d));
                if (azim_est >= angles_DOA_az(1) && azim_est <= angles_DOA_az(2) ...
                        &&elev_est >= angles_DOA_ele(1) && elev_est <= angles_DOA_ele(2))
                    angleObj_est(1,obj_cnt) = azim_est;
                    angleObj_est(2,obj_cnt) = elev_est;              
                    angleObj_est(3,obj_cnt) = ind;
                    angleObj_est(4,obj_cnt) = peakLoc_elev(j_elev);
                    %plot(angleObj_est(4,obj_cnt),angleObj_est(3,obj_cnt) ,'x','MarkerSize',12, 'LineWidth',2);
                    %hold on
                    obj_cnt = obj_cnt+1;
                else
                    continue;
                end
            end        
        end    
        %hold off
    end
end


%DOA_beamformingFFT_2D_joint.m
%DOA_beamformingFFT_2D_joint function perform 2D angle estimation based on FFT beamforming, both azimuth/elevation peak selection
%is done in 2D FFT domain
%input:
%   obj: object instance
%   sig: complex signal vector, with each value corresponding to each
%   antenna. The length of this vector equals to numTX x numRX enabled.
%   There can be overlapped antennas. this signal needs to be re-arranged
%   based on D value to form the virtual antenna array
%output:
%   angleObj_est: angle estimation results
%   angle_sepc_2D_fft: angle 2D fft spectrum
function [angleObj_est, angle_sepc_2D_fft ]= DOA_beamformingFFT_2D_joint(parameter_file_path, sig)

    sidelobeLevel_dB_azim = getPara(parameter_file_path, 'DOACascade_sidelobeLevel_dB_azim');
    %field of view to do beamforming
    angles_DOA_az = getPara(parameter_file_path, 'DOACascade_angles_DOA_az');
    angles_DOA_ele = getPara(parameter_file_path, 'DOACascade_angles_DOA_ele');
    %distance unit in terms of wavelength
    d = getPara(parameter_file_path, 'DOACascade_antDis');
    %2D matrix providing antenna coordinates
    D = getPara(parameter_file_path, 'DOACascade_D');
    angleFFTSize = getPara(parameter_file_path, 'DOACascade_DOAFFTSize');

    %FFT based implementation
    %first form a 2D matrix based on the antenna coordinates
    D = D + 1;
    apertureLen_azim = max(D(:,1));
    apertureLen_elev = max(D(:,2));
    sig_2D = zeros(apertureLen_azim,apertureLen_elev);
    for i_line = 1:apertureLen_elev
        ind = find(D(:,2) == i_line);
        D_sel = D(ind,1);
        sig_sel = sig(ind);
        [~, indU] = unique(D_sel);
        sig_2D(D_sel(indU),i_line) = sig_sel(indU);
    end

    %run FFT on azimuth and elevation
    angle_sepc_1D_fft=fftshift(fft(sig_2D,angleFFTSize,1),1); %corr1 is simply a 2D-fft of s1
    angle_sepc_2D_fft=fftshift(fft(angle_sepc_1D_fft,angleFFTSize,2),2); %corr1 is simply a 2D-fft of s1

    wx_vec=-pi : 2*pi/angleFFTSize : pi;
    wz_vec=-pi : 2*pi/angleFFTSize : pi;
    wx_vec = wx_vec(1:end-1);
    wz_vec = wz_vec(1:end-1);
    %convert to degree
    elev_vec = asind(wz_vec/(2*pi*d));

    angleObj_est= [];
    cnt = 1;
    if apertureLen_elev ==1
        %azimuth array only, no elevation antennas
        spec_azim = abs(angle_sepc_1D_fft(:,1));
        [~, peakLoc_azim] = DOA_BF_PeakDet_loc(parameter_file_path, spec_azim, sidelobeLevel_dB_azim);

        for i_obj = 1:length(peakLoc_azim)
            ind = peakLoc_azim(i_obj);
            azim_est = asind(wx_vec(ind)/(2*pi*d));
            if (azim_est >= angles_DOA_az(1) && azim_est <= angles_DOA_az(2))
                angleObj_est(1,cnt) = azim_est;
                angleObj_est(2,cnt) = 0;
                angleObj_est(3,cnt) = ind;
                angleObj_est(4,cnt) = 0;
                cnt = cnt+1;
            else
                continue;
            end
        end
    else
        %for each azimuth angle, find the max value in elevation direction, only
        %single target detection in elevation direction
        [spec_azim, peak_ele_ind]= max(abs(angle_sepc_2D_fft'));

        %detect multiple peaks in azimuth direction
        [peakVal_azim, peakLoc_azim] = DOA_BF_PeakDet_loc(parameter_file_path, spec_azim, sidelobeLevel_dB_azim);

        for i_obj = 1:length(peakVal_azim)
            elev_est = elev_vec(peak_ele_ind(peakLoc_azim(i_obj)));
            azim_sind = (wx_vec(peakLoc_azim(i_obj))/(2*pi*d))/cosd(elev_est);
            if (abs(azim_sind)<1)
                azim_est = asind(azim_sind);
                if ( azim_est>= angles_DOA_az(1) && azim_est <= angles_DOA_az(2) ...
                        &&elev_est >= angles_DOA_ele(1) && elev_est <= angles_DOA_ele(2))
                    angleObj_est(1,cnt) = azim_est;
                    angleObj_est(2,cnt) = elev_est;
                    angleObj_est(3,cnt) = peakLoc_azim(i_obj);
                    angleObj_est(4,cnt) = peak_ele_ind(peakLoc_azim(i_obj));
                    %plot(angleObj_est(4,cnt),angleObj_est(3,cnt) ,'x','MarkerSize',12, 'LineWidth',2);
                    %hold on
                    cnt = cnt+1;
                end        
            end
        end

    end

end


%DOA_BF_PeakDet_loc.m
% DOA_BF_PeakDet_loc function perform peak detection based on the input
% angle spectrum
%input:
%   obj: object instance
%   inData: angle spectrum
%output:
%   peakVal: value of detected peaks
%   peakLoc: index of detected peaks
function [peakVal, peakLoc] = DOA_BF_PeakDet_loc(parameter_file_path, inData, sidelobeLevel_dB)
    gamma = getPara(parameter_file_path, 'DOACascade_gamma');

    inData = inData(:);
    minVal = Inf;
    maxVal = 0;
    maxLoc = 0;
    maxData = [];
    locateMax = 0;  % at beginning, not ready for peak detection遇见满足条件的峰值启动开关，与参数gamma有关
    numMax = 0;
    extendLoc = 0;
    initStage = 1;
    absMaxValue = 0;
    i = 0;
    N = length(inData);
    while (i < (N + extendLoc - 1))
        i = i+1;
        i_loc = rem(i-1, N)+1;
        currentVal = inData(i_loc);
        % record the maximum value
        if currentVal > absMaxValue
            absMaxValue = currentVal;
        end
        % record the current max value and location
        if currentVal > maxVal
            maxVal = currentVal;
            maxLoc = i_loc;
            maxLoc_r = i;
        end

        % record for the current min value and location
        if currentVal < minVal
            minVal = currentVal;
        end

        if locateMax
            if currentVal < maxVal/gamma
                numMax = numMax + 1;
                bwidth = i - maxLoc_r;
                % Assign maximum value only if the value has fallen below the max by
                % gamma, thereby declaring that the max value was a peak
                maxData = [maxData(1:numMax-1,:) ; maxLoc maxVal bwidth maxLoc_r];

                minVal = currentVal;
                locateMax = 0;
            end
        else
            if currentVal > minVal*gamma
                % Assign minimum value if the value has risen above the min by
                % gamma, thereby declaring that the min value was a valley
                locateMax = 1;
                maxVal = currentVal;

                % aifuyuan add 20210608
                maxLoc = i_loc;
                maxLoc_r = i;

                if (initStage == 1)
                    extendLoc = i;
                    initStage = 0;
                end
            end
        end
    end


    % make sure the max value needs to be cetain dB higher than the side lobes to declare any detection

    % [v ind] = max(maxData(:,2));
    % peakMean = mean(maxData([1:(ind-1) ind+1:end],2));
    % SNR_DOA = v/peakMean;
    % if v>peakMean*maxPeakThre
    % if the max is different by more than sidelobeLevel_dB dB from the
    % peak, then removed it as a sidelobe
    maxData_ = [];
    numMax_ = 0;
    totPower = 0;
    for i = 1:numMax
        if maxData(i, 2) >= absMaxValue * (10^(-sidelobeLevel_dB/10))
            numMax_ = numMax_ + 1;
            maxData_(numMax_,:) = maxData(i, :);
            totPower = totPower + maxData(i, 2);
        end
    end
    maxData = maxData_;
    numMax = numMax_;

    peakVal = zeros(numMax, 1);
    peakLoc = zeros(numMax, 1);
    for ind = 1:numMax
        peakVal(ind) = maxData(ind,2);
        peakLoc(ind) = rem(maxData(ind,1)-1, N)+1;
    end
    
end

