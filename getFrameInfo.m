function curFrameInfo = getFrameInfo(sliceIdBasedPath, sliceId, frameId)
    curFrameInfo = struct();
    curFrameInfo.sliceFrameId = frameId;

    masterData_path = sliceIdBasedPath(sliceId).master_data;
    curFrameInfo.master_adcDataPath = masterData_path;
    masterIdx_path = sliceIdBasedPath(sliceId).master_idx;
    [~, ~,  ~, ~, ~, ~, curFrameInfo.master_timestamp, curFrameInfo.master_offset] = getFrameInfoFromIdx(masterIdx_path, frameId);
    
    slave1Data_path = sliceIdBasedPath(sliceId).slave1_data;
    curFrameInfo.slave1_adcDataPath = slave1Data_path;
    slave1Idx_path = sliceIdBasedPath(sliceId).slave1_idx;
    [~, ~,  ~, ~, ~, ~, curFrameInfo.slave1_timestamp, curFrameInfo.slave1_offset] = getFrameInfoFromIdx(slave1Idx_path, frameId);
    
    slave2Data_path = sliceIdBasedPath(sliceId).slave2_data;
    curFrameInfo.slave2_adcDataPath = slave2Data_path;
    slave2Idx_path = sliceIdBasedPath(sliceId).slave2_idx;
    [~, ~,  ~, ~, ~, ~, curFrameInfo.slave2_timestamp, curFrameInfo.slave2_offset] = getFrameInfoFromIdx(slave2Idx_path, frameId);
    
    slave3Data_path = sliceIdBasedPath(sliceId).slave3_data;
    curFrameInfo.slave3_adcDataPath = slave3Data_path;
    slave3Idx_path = sliceIdBasedPath(sliceId).slave3_idx;
    [~, ~,  ~, ~, ~, ~, curFrameInfo.slave3_timestamp, curFrameInfo.slave3_offset] = getFrameInfoFromIdx(slave3Idx_path, frameId);
    
    
end


function [tag, version,  flags, width, height, size, timestamp, offset] = getFrameInfoFromIdx(path, frameId)
    f = fopen(path, 'r');
    % File header in *_idx.bin:
    %     struct Info
    %     {
    %         uint32_t tag;
    %         uint32_t version;
    %         uint32_t flags;
    %         uint32_t numIdx;       // number of frames 
    %         uint64_t dataFileSize; // total data size written into file
    %     };
    % 
    % Index for every frame from each radar:
    %     struct BuffIdx
    %     {
    %         uint16_t tag;
    fseek(f, 24 + frameId * 48 + 0, 'bof');
    tag = uint16(fread(f, 1,'uint16'));
    %         uint16_t version; /*same as Info.version*/
    fseek(f, 24 + frameId * 48 + 2, 'bof');
    version = uint16(fread(f, 1,'uint16'));
    %         uint32_t flags;
    fseek(f, 24 + frameId * 48 + 4, 'bof');
    flags = uint32(fread(f, 1,'uint32'));
    %         uint16_t width;        
    fseek(f, 24 + frameId * 48 + 8, 'bof');
    width = uint16(fread(f, 1,'uint16'));
    %         uint16_t height;
    fseek(f, 24 + frameId * 48 + 10, 'bof');
    height = uint16(fread(f, 1,'uint16'));
    %         uint32_t pitchOrMetaSize[4]; /*For image data, this is pitch.
    %                                                        For raw data, this is size in bytes per metadata plane.*/
    %         uint32_t size; /*total size in bytes of the data in the buffer (sum of all planes)*/
    fseek(f, 24 + frameId * 48 + 28, 'bof');
    size = uint32(fread(f, 1,'uint32'));
    %         uint64_t timestamp;
    fseek(f, 24 + frameId * 48 + 32, 'bof');
    timestamp = uint64(fread(f, 1,'uint64'));
    %         uint64_t offset;
    fseek(f, 24 + frameId * 48 + 40, 'bof');
    offset = uint64(fread(f, 1,'uint64'));
    %     };
    fclose(f);
end