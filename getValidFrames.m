function num_validFrames = getValidFrames(sliceIdBasedPath, sliceId)
    masterIdx_path = sliceIdBasedPath(sliceId).master_idx;
    f = fopen(masterIdx_path, 'r');
    fseek(f, 12, 'bof');
    num_validFrames = uint32(fread(f, 1,'uint32'));
    fclose(f);
end