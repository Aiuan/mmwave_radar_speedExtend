function startTimestamp_image_fixed = radarTimestampMatchCameraTimestamp(imageFolder, imageType, timestamp_radar)
    imagesInfo = dir(fullfile(imageFolder,['*',imageType]));
    imagesTimestamp = [];
    
    for i = 1: length(imagesInfo)
        temp = split(imagesInfo(i).name, '.');
        imagesTimestamp = [imagesTimestamp; str2num(temp{1})];  
    end
    
    timestamp_diffs = imagesTimestamp - double(timestamp_radar);
    [~, index] = min(abs(timestamp_diffs));
    
    startTimestamp_image_fixed = uint64(imagesTimestamp(index));
    
end