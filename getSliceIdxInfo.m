function sliceIdxInfo = getSliceIdxInfo(sliceIdBasedPath, sliceId)
    masterIdx_path = sliceIdBasedPath(sliceId).master_idx;
    slave1Idx_path = sliceIdBasedPath(sliceId).slave1_idx;
    slave2Idx_path = sliceIdBasedPath(sliceId).slave2_idx;
    slave3Idx_path = sliceIdBasedPath(sliceId).slave3_idx;    
    
    master = struct();
    slave1 = struct();
    slave2 = struct();
    slave3 = struct();
    
    [master.num_frames, master.size_file, master.framesInfo] = getInfoFromIdx(masterIdx_path);
    [slave1.num_frames, slave1.size_file, slave1.framesInfo] = getInfoFromIdx(slave1Idx_path);
    [slave2.num_frames, slave2.size_file, slave2.framesInfo] = getInfoFromIdx(slave2Idx_path);
    [slave3.num_frames, slave3.size_file, slave3.framesInfo] = getInfoFromIdx(slave3Idx_path);
    
    sliceIdxInfo = struct();
    sliceIdxInfo.master = master;
    sliceIdxInfo.slave1 = slave1;
    sliceIdxInfo.slave2 = slave2;
    sliceIdxInfo.slave3 = slave3;
    
    
end



function [num_frames, size_file, framesInfo] = getInfoFromIdx(path)
    f = fopen(path, 'r');
    % File header in *_idx.bin:
    %     struct Info
    %     {
    %         uint32_t tag;
    %         uint32_t version;
    %         uint32_t flags;
    %         uint32_t numIdx;       // number of frames 
    fseek(f, 12, 'bof');
    num_frames = uint32(fread(f, 1,'uint32'));
    %         uint64_t dataFileSize; // total data size written into file
    fseek(f, 16, 'bof');
    size_file = uint64(fread(f, 1,'uint64'));
    %     };
    % 
    framesInfo = struct();
    % Index for every frame from each radar:
    for i = 0:num_frames-1
        framesInfo(i+1).frameId = i;
    %     struct BuffIdx
    %     {
    %         uint16_t tag;
        fseek(f, 24 + i * 48 + 0, 'bof');
        framesInfo(i+1).tag = uint16(fread(f, 1,'uint16'));
    %         uint16_t version; /*same as Info.version*/
        fseek(f, 24 + i * 48 + 2, 'bof');
        framesInfo(i+1).version = uint16(fread(f, 1,'uint16'));
    %         uint32_t flags;
        fseek(f, 24 + i * 48 + 4, 'bof');
        framesInfo(i+1).flags = uint32(fread(f, 1,'uint32'));
    %         uint16_t width;        
        fseek(f, 24 + i * 48 + 8, 'bof');
        framesInfo(i+1).width = uint16(fread(f, 1,'uint16'));
    %         uint16_t height;
        fseek(f, 24 + i * 48 + 10, 'bof');
        framesInfo(i+1).height = uint16(fread(f, 1,'uint16'));
    %         uint32_t pitchOrMetaSize[4]; /*For image data, this is pitch.
    %                                                        For raw data, this is size in bytes per metadata plane.*/
    %         uint32_t size; /*total size in bytes of the data in the buffer (sum of all planes)*/
        fseek(f, 24 + i * 48 + 28, 'bof');
        framesInfo(i+1).size = uint32(fread(f, 1,'uint32'));
    %         uint64_t timestamp;
        fseek(f, 24 + i * 48 + 32, 'bof');
        framesInfo(i+1).timestamp = uint64(fread(f, 1,'uint64'));
    %         uint64_t offset;
        fseek(f, 24 + i * 48 + 40, 'bof');
        framesInfo(i+1).offset = uint64(fread(f, 1,'uint64'));
    
    %     };
    end   
    
    fclose(f);
end