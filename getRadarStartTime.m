function radarStartTime = getRadarStartTime(path)
    f = fopen(path, 'r');
%     if skipOneSubframe_ON
%         subFrameId = 1;
%     else
%         subFrameId = 0;
%     end
    subFrameId = 0;
    %         uint64_t timestamp;
    fseek(f, 24 + subFrameId * 48 + 32, 'bof');
    radarStartTime = uint64(fread(f, 1,'uint64'));
    
    fclose(f);
end