function radarStartTime = getPCStartTime(radarTimeFile)
    f = fopen(fullfile(radarTimeFile.folder, radarTimeFile.name), 'r');
    radarStartTime = fscanf(f, '%f');
    fclose(f);
    radarStartTime = uint64(radarStartTime*1000000);%16位UNIX时间
end