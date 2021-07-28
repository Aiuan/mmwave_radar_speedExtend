function [sliceIdBasedPath] = getSliceIdBasedPath(radarBinFile_list)
    
    if mod(length(radarBinFile_list), 8) ~= 0
        disp('####Error：.bin文件数量不是8的整数倍');
    end
    num_slices = ceil(length(radarBinFile_list)/8);
    sliceIdBasedPath = struct();
    for i =1: num_slices
        sliceIdBasedPath(i).sliceId = sprintf('%04d', i-1);
    end
    
    for i = 1:length(radarBinFile_list)
        path = fullfile(radarBinFile_list(i).folder, radarBinFile_list(i).name);
        nameSplitResults = strsplit(radarBinFile_list(i).name, '_');
        device = nameSplitResults{1};
        sliceId = nameSplitResults{2};
        type = nameSplitResults{3};
        type = strsplit(type, '.');
        type = type{1};        
        eval(['sliceIdBasedPath(str2num(sliceId)+1).', device, '_', type, ' = ', 'path;']);
        
    end
    
    
end