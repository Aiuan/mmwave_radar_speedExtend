function [detection_results] = CFAR(input, parameter_file_path)
    detectMethod = getPara(parameter_file_path, 'CFAR_OS_detectMethod');
    rangeBinSize = getPara(parameter_file_path, 'CFAR_OS_rangeBinSize');
    dopplerFFTSize = getPara(parameter_file_path, 'CFAR_OS_dopplerFFTSize');
    velocityBinSize = getPara(parameter_file_path, 'CFAR_OS_velocityBinSize');
    numAntenna = getPara(parameter_file_path, 'CFAR_OS_numAntenna');
    minDisApplyVmaxExtend = getPara(parameter_file_path, 'CFAR_OS_minDisApplyVmaxExtend');
    overlapAntenna_ID = getPara(parameter_file_path, 'CFAR_OS_overlapAntenna_ID');
    TDM_MIMO_numTX = getPara(parameter_file_path, 'CFAR_OS_TDM_MIMO_numTX');
    numRxAnt = getPara(parameter_file_path, 'CFAR_OS_numRxAnt');
    antenna_azimuthonly = getPara(parameter_file_path, 'CFAR_OS_antenna_azimuthonly');
    applyVmaxExtend = getPara(parameter_file_path, 'CFAR_OS_applyVmaxExtend');
    

    sig_integrate = sum((abs(input)).^2,3) + 1; %沿天线阵列组合非相干信号变成2维数据

    angleFFTSize = 128; 
    angleBinSkipLeft = 4;
    angleBinSkipRight = 4;

    if  detectMethod
        [N_obj_Rag, Ind_obj_Rag, noise_obj, CFAR_SNR] = CFAR_OS_Range(parameter_file_path, sig_integrate); %进行距离维上目标检测
        %距离维检测返回检测目标数目和检测目标索引以及门限水平和SNR
        N_obj = 0; %多普勒维目标数目初始化
        Ind_obj = []; %多普勒维目标索引初始化
        detection_results = {}; %检测结果初始化
        if (N_obj_Rag>0) %如果距离维检测到目标
            [N_obj, Ind_obj] = CFAR_OS_Doppler_overlap(parameter_file_path, Ind_obj_Rag, input, sig_integrate); %则进行多普勒维检测
            detection_results = [];

            %使用首次噪声估计并应用于第二次检测
            noise_obj_agg = [];
            for i_obj = 1:N_obj %对多普勒维的每个目标
                indx1R = Ind_obj(i_obj,1); %返回多普勒维目标索引
                indx1D = Ind_obj(i_obj,2); %返回多普勒维循环次数
                ind2R = find(Ind_obj_Rag(:,1) == indx1R); 
                ind2D = find(Ind_obj_Rag(ind2R,2) == indx1D);
                noiseInd = ind2R(ind2D);
                noise_obj_agg(i_obj) = noise_obj(noiseInd);
            end

            for i_obj = 1:N_obj
                xind = (Ind_obj(i_obj,1)-1) +1;
                detection_results(i_obj).rangeInd = Ind_obj(i_obj, 1) - 1;  %range index
                detection_results(i_obj).range = (detection_results(i_obj).rangeInd) * rangeBinSize;  %range estimation
                dopplerInd  = Ind_obj(i_obj, 2) - 1;  %Doppler index
                detection_results(i_obj).dopplerInd_org = dopplerInd;
                detection_results(i_obj).dopplerInd = dopplerInd;

                %velocity estimation
                detection_results(i_obj).doppler = (dopplerInd-dopplerFFTSize/2)*velocityBinSize;
                detection_results(i_obj).doppler_corr = detection_results (i_obj).doppler;
                detection_results(i_obj).noise_var = noise_obj_agg(i_obj);       %noise variance
                detection_results(i_obj).bin_val  = reshape(input(xind, Ind_obj(i_obj,2),:),numAntenna,1);  %2d FFT value for the 4 antennas
                %detection_results(i_obj).estSNR  = 10*log10(sum(abs(detection_results (i_obj).bin_val).^2)/sum(detection_results (i_obj).noise_var));  %2d FFT value for the 4 antennas
                detection_results(i_obj).estSNR  = (sum(abs(detection_results (i_obj).bin_val).^2)/sum(detection_results (i_obj).noise_var));  

                sig_bin = [];
                %only apply max velocity extention if it is enabled and distance is larger
                %than minDisApplyVmaxExtend
                if (applyVmaxExtend == 1 && (detection_results(i_obj).range > minDisApplyVmaxExtend) && (~isempty(overlapAntenna_ID)))
                    velocityObj_est = detection_results(i_obj).doppler;
                    if mod(TDM_MIMO_numTX,2)==1
                        %odd number
                        dopplerInd_unwrap = dopplerInd + ((1:TDM_MIMO_numTX)-ceil(TDM_MIMO_numTX/2))*dopplerFFTSize;

                    else
                        %even number
                        if velocityObj_est>0
                            dopplerInd_unwrap = dopplerInd + ((1:TDM_MIMO_numTX)-(TDM_MIMO_numTX/2+1))*dopplerFFTSize;

                        else
                            dopplerInd_unwrap = dopplerInd + ((1:TDM_MIMO_numTX)-TDM_MIMO_numTX/2)*dopplerFFTSize;

                        end
                    end
                    sig_bin_org = detection_results (i_obj).bin_val;                
                    %Doppler phase correction due to TDM MIMO             
                    deltaPhi = 2*pi*(dopplerInd_unwrap-dopplerFFTSize/2)/( TDM_MIMO_numTX*dopplerFFTSize);

                    % construct all possible signal vectors based on the number
                    % of possible hypothesis
                    for i_TX = 1:TDM_MIMO_numTX
                        RX_ID = (i_TX-1)*numRxAnt+1 : i_TX*numRxAnt;
                        sig_bin(RX_ID,: )= sig_bin_org(RX_ID )* exp(-1j*(i_TX-1)*deltaPhi);
                    end

                    % use overlap antenna to do max velocity unwrap                
                    signal_overlap = sig_bin_org(overlapAntenna_ID(:,1:2));                

                    %check the phase difference of each overlap antenna pair
                    %for each hypothesis
                    angle_sum_test = [];                
                    for i_sig = 1:size(signal_overlap,1)
                        for i_test = 1:length(deltaPhi)
                            signal2 = signal_overlap(1:i_sig,2)*exp(-j*deltaPhi(i_test));
                            angle_sum_test(i_sig,i_test) = angle(sum(signal_overlap(1:i_sig,1).*conj(signal2)));
                        end
                    end

                    %chosee the hypothesis with minimum phase difference to
                    %estimate the unwrap factor
                    [val_doppler_unwrap_integ_overlap doppler_unwrap_integ_overlap] = min(abs(angle_sum_test),[],2);


                    %test the angle FFT SNR
                    sig_bin_row1 = sig_bin(antenna_azimuthonly,:);
                    sig_bin_row1_fft = fftshift(fft(sig_bin_row1,angleFFTSize),1);
                    sig_bin_row1_fft_cut = abs(sig_bin_row1_fft(angleBinSkipLeft+1:(angleFFTSize-angleBinSkipRight),:));
                    [val doppler_unwrap_integ_FFT] = max(max(sig_bin_row1_fft_cut));


                    b = unique(doppler_unwrap_integ_overlap);
                    c = histc(doppler_unwrap_integ_overlap(:),b);
                    [val ind] = max(c);
                    doppler_unwrap_integ_overlap_sel = b(ind);
                    doppler_unwrap_integ = doppler_unwrap_integ_overlap_sel;           


                    %overlap antenna method is applied by default
                    detection_results(i_obj).bin_val = sig_bin(:,doppler_unwrap_integ);                  

                    %corret velocity after applying the integer value
                    dopplerInd = dopplerInd_unwrap(doppler_unwrap_integ);
                    dopplerInd_FFT = dopplerInd_unwrap(doppler_unwrap_integ_FFT);
                    dopplerInd_overlap = dopplerInd_unwrap(doppler_unwrap_integ_overlap_sel);
                    detection_results (i_obj).dopplerInd = dopplerInd;
                    %velocity estimation
                    detection_results (i_obj).doppler_corr = (dopplerInd-dopplerFFTSize/2)*velocityBinSize;                
                    %both overlap antenna and FFT results are reported 
                    detection_results(i_obj).doppler_corr_overlap = (dopplerInd_overlap-dopplerFFTSize/2)*velocityBinSize;
                    detection_results(i_obj).doppler_corr_FFT = (dopplerInd_FFT-dopplerFFTSize/2)*velocityBinSize;
                    detection_results(i_obj).overlapTests = doppler_unwrap_integ_overlap;
                    detection_results(i_obj).overlapTestsVal = val_doppler_unwrap_integ_overlap;
                else
                    %Doppler phase correction due to TDM MIMO without apply
                    %Vmax extention algorithm

                    deltaPhi = 2*pi*(dopplerInd-dopplerFFTSize/2)/( TDM_MIMO_numTX*dopplerFFTSize);
                    sig_bin_org = detection_results (i_obj).bin_val;
                    for i_TX = 1:TDM_MIMO_numTX
                        RX_ID = (i_TX-1)*numRxAnt+1 : i_TX*numRxAnt;
                        sig_bin(RX_ID,: )= sig_bin_org(RX_ID )* exp(-1j*(i_TX-1)*deltaPhi);
                    end
                    detection_results(i_obj).bin_val = sig_bin;
                    detection_results(i_obj).doppler_corr_overlap = detection_results(i_obj).doppler_corr;
                    detection_results(i_obj).doppler_corr_FFT = detection_results(i_obj).doppler_corr;

                end

            end
        end
        
    end

end


function [N_obj, Ind_obj, noise_obj, CFAR_SNR] = CFAR_OS_Range(parameter_file_path, sig)
    % 本函数为OS CFAR距离维CFAR
    %输入的sig为2D数据，分别表示每个chrip采样点数，Loop次数
    cellNum = getPara(parameter_file_path, 'CFAR_OS_refWinSize');%读取对象的参考窗
    gapNum = getPara(parameter_file_path, 'CFAR_OS_guardWinSize'); %读取对象的保护单元
    cellNum = cellNum(1); %距离维CFAR的参考窗长度
    gapNum = gapNum(1); %距离维CFAR的保护单元
    K0 = getPara(parameter_file_path, 'CFAR_OS_K0'); %距离维的门限系数
    K0 = K0(1);

    M_samp = size(sig, 1); %获取采样点数
    N_pul = size(sig, 2); %获取Loop次数

    gaptot = gapNum + cellNum; %参考单元+保护单元（相对半参考窗而言）
    N_obj = 0; %目标数目
    Ind_obj = []; %目标在sig中的索引
    noise_obj = []; %噪声存储
    CFAR_SNR = []; %CFAR信噪比存储

    discardCellLeft = getPara(parameter_file_path, 'CFAR_OS_discardCellLeft'); %丢弃左侧单元的数目
    discardCellRight = getPara(parameter_file_path, 'CFAR_OS_discardCellRight'); %丢弃右侧单元的数目
    maxEnable = getPara(parameter_file_path, 'CFAR_OS_maxEnable'); %是否检测峰值最大状态

    for k = 1:N_pul %每次循环
        sigv = (sig(:,k))'; %提取每次循环的采样点数并转换为行向量
        vec = sigv(discardCellLeft+1:M_samp-discardCellRight); %提取去除近场和远场后的数据
        vecLeft = vec(1:(gaptot)); %提取左侧近场背景水平填补左边界
        vecRight = vec(end-(gaptot)+1:end);  %提取右侧远场背景水平填补右边界
        vec = [vecLeft vec vecRight]; %构建闭环检测单元
        for j = 1:(M_samp-discardCellLeft-discardCellRight)
            cellInd = [j-gaptot:j-gapNum-1 j+gapNum+1:j+gaptot]; %获取j检测单元的参考窗
            cellInd = cellInd + gaptot; %加上填补的边界长度

            cellInda = [j-gaptot:j-gapNum-1]; %左侧半参考窗
            cellInda = cellInda + gaptot; %加上填补的左边界长度
            cellIndb = [j+gapNum+1:j+gaptot]; %右侧参考窗
            cellIndb = cellIndb + gaptot; %加上填补的右边界长度
    %         cellave1a =sum(vec(cellInda))/(cellNum); %对左侧参考窗采样值取平均
    %         cellave1b =sum(vec(cellIndb))/(cellNum); %对右侧参考窗蚕养殖取平均
    %         cellave1 = min(cellave1a,cellave1b); %选择平均值最小的一侧参考窗值
            cell_k = sort(vec([cellInda,cellIndb])); %进行排序
            cellave1 = cell_k(ceil(3/4*cellNum*2)); %选择第3/4参考窗的单元作为估计

            %if((j > discardCellLeft) && (j < (M_samp-discardCellRight)))
            if maxEnable == 1 %进行峰值判断
                maxInCell = max(vec([cellInd(1):cellInd(end)])); %检查是否为局部最大值并输出
                if (vec(j+gaptot)>K0*cellave1 && ( vec(j+gaptot) >= maxInCell)) 
                    %j+gaptot是指原始vec有效单元中的目标，判断检测单元是否大于门限水平
                    N_obj = N_obj+1; %检测单元大于门限水平则目标数目+1
                    Ind_obj(N_obj,:) = [j+discardCellLeft, k]; %记录下目标的位置（这里要加上删除的近场数据）
                    noise_obj(N_obj) = cellave1; %保存杂波功率水平
                    CFAR_SNR(N_obj) = vec(j+gaptot)/cellave1; %记录目标信噪比
                end
            else %不进行峰值判断，其余步骤与上相同
                if vec(j+gaptot)>K0*cellave1
                    N_obj=N_obj+1;
                    Ind_obj(N_obj,:)=[j+discardCellLeft, k];
                    noise_obj(N_obj) = cellave1; %保存噪声水平
                    CFAR_SNR(N_obj) = vec(j+gaptot)/cellave1;
                end
            end        
        end
    end


    %获取每个阵列的噪声方差系数
    for i_obj = 1:N_obj %每个目标
        ind_range = Ind_obj(i_obj,1); %目标在sig的索引
        ind_Dop = Ind_obj(i_obj,2); %第几个循环
        if ind_range <= gaptot %左边界参考窗处理
            cellInd = [ind_range+gapNum+1:ind_range+gaptot ind_range+gapNum+1:ind_range+gaptot];
        elseif ind_range >= M_samp-gaptot+1 %右边界参考窗处理
            cellInd = [ind_range-gaptot:ind_range-gapNum-1 ind_range-gaptot:ind_range-gapNum-1];
        else%中间边界参考窗处理
            cellInd = [ind_range-gaptot:ind_range-gapNum-1 ind_range+gapNum+1:ind_range+gaptot];   
        end
    end
    
end


function [N_obj, Ind_obj, noise_obj_an] = CFAR_OS_Doppler_overlap(parameter_file_path, Ind_obj_Rag, sigCpml, sig_integ)
    % 本函数为OS CFAR多普勒维CFAR
    maxEnable = getPara(parameter_file_path, 'CFAR_OS_maxEnable'); %是否检测峰值最大状态
    cellNum0 = getPara(parameter_file_path, 'CFAR_OS_refWinSize'); %读取对象的参考窗长度
    gapNum0 = getPara(parameter_file_path, 'CFAR_OS_guardWinSize'); %读取对象的保护单元长度
    cellNum = cellNum0(2); %读取多普勒维的参考窗长度
    gapNum = gapNum0(2); %读取多普勒维的保护窗长度
    K0 = getPara(parameter_file_path, 'CFAR_OS_K0'); %多普勒维的门限系数
    K0 = K0(2);

    rangeNumBins = size(sig_integ,1); %输入的sig为2D数据，分别表示每个chrip采样点数（距离维概念），Loop次数（阵列概念）

    detected_Rag_Cell = unique(Ind_obj_Rag(:,1)); %获取不重复的目标距离索引
    sig = sig_integ(detected_Rag_Cell,:); %统计所有目标的循环情况即radarcube的列向量

    M_samp = size(sig, 1); %获取采样点数，及距离维检测到目标的点数
    N_pul = size(sig, 2); %获取Loop次数

    gaptot = gapNum + cellNum; %参考窗

    N_obj = 0; %检测目标数初始化
    Ind_obj = []; %目标索引初始化
    noise_obj_an = []; 
    vec = zeros(1,N_pul+gaptot*2);
    for k = 1:M_samp %对每个目标点
        detected_Rag_Cell_i = detected_Rag_Cell(k); %获取第k个目标点的距离
        ind1 = find(Ind_obj_Rag(:,1) == detected_Rag_Cell_i);%第一次循环中目标位置与目标点的距离进行匹配
        indR = Ind_obj_Rag(ind1, 2); %返回索引

        sigv=(sig(k,:)); %获取每个目标距离对应的多普勒维信息
        vec(1:gaptot) = sigv(end-gaptot+1:end); %将右边界信息拷贝到左边界
        vec(gaptot+1: N_pul+gaptot) = sigv; %中间信息用sigv拷贝
        vec(N_pul+gaptot+1:end) = sigv(1:gaptot); %将左边界信息拷贝到右边界

        %CFAR处理
        ind_loc_all = []; %距离维索引矩阵
        ind_loc_Dop = []; %多普勒维索引矩阵
        ind_obj_0 = 0;
        noiseEst = zeros(1,N_pul); %噪声估计，长度为阵列数
        for j = 1+gaptot:N_pul+gaptot %去中间原始的数据进行处理
            cellInd = [j-gaptot:j-gapNum-1 j+gapNum+1:j+gaptot]; %获取参考窗
            noiseEst(j-gaptot) = sum(vec(cellInd)); %对参考窗内数据求和
        end
        for j = 1+gaptot:N_pul+gaptot
            j0 = j - gaptot;
            cellInd = [j-gaptot:j-gapNum-1 j+gapNum+1:j+gaptot];
            cellInda = [j-gaptot: j-gapNum-1]; %取左侧参考窗
            cellIndb =[j+gapNum+1:j+gaptot]; %取右侧参考窗

    %         cellave1a = sum(vec(cellInda))/(cellNum); %左侧参考窗多普勒维数据平均
    %         cellave1b = sum(vec(cellIndb))/(cellNum); %右侧参考窗多普勒维数据平均
    %         cellave1 = min(cellave1a,cellave1b); %两侧参考窗平均取最小        
            cell_k = sort(vec([cellInda,cellIndb]));
            cellave1 = cell_k(ceil(3/4*cellNum*2));
            maxInCell = max(vec(cellInd)); %获取峰值数据
            if maxEnable == 1 %检测该单元是否为参考窗内的峰值
                condition = ((vec(j)>K0*cellave1)) && ((vec(j)>maxInCell));
            else
                condition = vec(j)>K0*cellave1;
            end

            if condition == 1 %如果j单元是参考窗的峰值
                %检查该检测是否与多普勒检测重叠
                if(find(indR == j0))
                    %如果重叠则声明检测结果
                    ind_win = detected_Rag_Cell_i; %返回j0目标的多普勒维信息
                    ind_loc_all = [ind_loc_all ind_win]; 
                    ind_loc_Dop = [ind_loc_Dop j0]; %返回j0目标的索引，对应于距离维距离
                end
            end
        end
        ind_obj_0 = [];

        if (length(ind_loc_all)>0)
            ind_obj_0(:,1) = ((ind_loc_all)); %indobj第一列存放数据
            ind_obj_0(:,2) = ind_loc_Dop; %indobj第二列存放距离索引
            if size(Ind_obj,1) == 0 %如果没有数据则赋值给IndObj
                Ind_obj = ind_obj_0;
            else%否则
                %以下过程是为了避免重复的检测点
                ind_obj_0_sum = ind_loc_all + 10000 * ind_loc_Dop;
                Ind_obj_sum = Ind_obj(:,1) + 10000 * Ind_obj(:,2);
                for ii= 1: length(ind_loc_all)
                    if (length(find(Ind_obj_sum == ind_obj_0_sum(ii)))==0)
                        Ind_obj = [Ind_obj ; ind_obj_0(ii,:)];
                    end
                end
            end
        end
    end

    N_obj = size(Ind_obj,1);

    %重置参考窗口
    cellNum = cellNum0(1); %距离维窗口长度
    gapNum = gapNum0(1); %距离维保护单元
    gaptot = gapNum + cellNum;

    %获取每个阵列的噪声系数
    N_obj_valid = 0;
    Ind_obj_valid = [];
    powerThre = getPara(parameter_file_path, 'CFAR_OS_powerThre');
    numAntenna = getPara(parameter_file_path, 'CFAR_OS_numAntenna');
    for i_obj = 1:N_obj %对每个目标    
        ind_range = Ind_obj(i_obj,1);
        ind_Dop = Ind_obj(i_obj,2);
        %跳过信号功率小于powerThre的检测点
        if (min(abs(sigCpml(ind_range, ind_Dop,:)).^2) < powerThre)
            continue;
        end
        if ind_range <= gaptot%右边界
            cellInd = [ind_range+gapNum+1:ind_range+gaptot ind_range+gapNum+1:ind_range+gaptot];
        elseif ind_range >= rangeNumBins-gaptot+1%左边界
            cellInd = [ind_range-gaptot:ind_range-gapNum-1 ind_range-gaptot:ind_range-gapNum-1];
        else %中间边界
            cellInd = [ind_range-gaptot: ind_range-gapNum-1 ind_range+gapNum+1:ind_range+gaptot];
        end

        N_obj_valid = N_obj_valid +1;
        noise_obj_an(:, i_obj) = reshape((mean(abs(sigCpml(cellInd, ind_Dop, :)).^2, 1)), numAntenna, 1, 1);
        Ind_obj_valid(N_obj_valid,:) = Ind_obj(i_obj,:);    

    end

    N_obj = N_obj_valid;
    Ind_obj = Ind_obj_valid;
end