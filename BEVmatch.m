function [results] = BEVmatch(subFrame_results, epsilon, minpts, maxMatchErr)
    results = struct();

    % 原始
    figure();
    subplot(2,2,[1,2]);
    xyz = subFrame_results(1).xyz;
    scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 10, (xyz(:, 4)),'filled');%x,y,z,v
    hold on;
    xyz = subFrame_results(2).xyz;
    scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 10, (xyz(:, 4)),'o');%x,y,z,v
    hold off;    
    c = colorbar;
    c.Label.String = 'velocity (m/s)'; 
    grid on;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    colormap('jet');
    set(gca, 'Color', [0.8,0.8,0.8]);
    axis(gca, 'equal');
    title('原始 子帧点云效果聚合图');
    legend('subFrame0', 'subFrame1');
    
    %% 原始点云聚类结果    
    %subframe 0
    subFrame0ObjInfo = struct();
    cnt_subFrame0obj = 0;
    xyz = subFrame_results(1).xyz;
    idx = dbscan(xyz(:, 1:2),epsilon,minpts);
    subplot(2,2,3);
    classIds = unique(idx);
    legend_list = {};
    for i = 1:length(classIds)
        legend_list{i} = num2str(classIds(i));
        mask = (idx == classIds(i));
        if classIds(i) == -1
            scatter3(xyz(mask, 1), xyz(mask, 2), xyz(mask, 3), 10, '*');
        else
            scatter3(xyz(mask, 1), xyz(mask, 2), xyz(mask, 3), 10, 'filled');
            cnt_subFrame0obj = cnt_subFrame0obj + 1;
            subFrame0ObjInfo(cnt_subFrame0obj).xyz = xyz(mask, :);
            subFrame0ObjInfo(cnt_subFrame0obj).objX = mean(xyz(mask, 1));
            subFrame0ObjInfo(cnt_subFrame0obj).objY = mean(xyz(mask, 2));
            subFrame0ObjInfo(cnt_subFrame0obj).objZ = mean(xyz(mask, 3));
        end
        hold on;
    end
    hold off;
    set(gca, 'Color', [0.8,0.8,0.8]);
    axis(gca, 'equal');
    title('原始 subframe0点云聚类结果');
    legend(legend_list);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    results.subFrame0ObjInfo = subFrame0ObjInfo;
    
    
    %subframe 1
    subFrame1ObjInfo = struct();
    cnt_subFrame1obj = 0;
    xyz = subFrame_results(2).xyz;
    idx = dbscan(xyz(:, 1:2),epsilon,minpts);
    subplot(2,2,4);
    classIds = unique(idx);
    legend_list = {};
    for i = 1:length(classIds)
        legend_list{i} = num2str(classIds(i));
        mask = (idx == classIds(i));
        if classIds(i) == -1
            scatter3(xyz(mask, 1), xyz(mask, 2), xyz(mask, 3), 10, '*');
        else
            scatter3(xyz(mask, 1), xyz(mask, 2), xyz(mask, 3), 10, 'filled');
            cnt_subFrame1obj = cnt_subFrame1obj + 1;
            subFrame1ObjInfo(cnt_subFrame1obj).xyz = xyz(mask, :);
            subFrame1ObjInfo(cnt_subFrame1obj).objX = mean(xyz(mask, 1));
            subFrame1ObjInfo(cnt_subFrame1obj).objY = mean(xyz(mask, 2));
            subFrame1ObjInfo(cnt_subFrame1obj).objZ = mean(xyz(mask, 3));
        end
        hold on;
    end
    hold off;
    set(gca, 'Color', [0.8,0.8,0.8]);
    axis(gca, 'equal');
    title('原始 subframe1点云聚类结果');
    legend(legend_list);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    results.subFrame1ObjInfo = subFrame1ObjInfo;
    
    %% XYZ空间内匹配    
    matchInfo = struct();
    cnt_matchObj = 0;
    for sub0ObjId = 1:cnt_subFrame0obj
        errorDis_list = zeros(1, cnt_subFrame1obj);
        for sub1ObjId = 1:cnt_subFrame1obj
            errorDis_list(sub1ObjId) = sqrt((subFrame0ObjInfo(sub0ObjId).objX - subFrame1ObjInfo(sub1ObjId).objX).^2 + ...
                (subFrame0ObjInfo(sub0ObjId).objY - subFrame1ObjInfo(sub1ObjId).objY).^2);
        end
        
        if min(errorDis_list) < maxMatchErr
            cnt_matchObj = cnt_matchObj + 1;
            matchInfo(cnt_matchObj).subFrame0ObjId = sub0ObjId;
            matchInfo(cnt_matchObj).subFrame0_xyz = subFrame0ObjInfo(sub0ObjId).xyz;
            matchInfo(cnt_matchObj).subFrame0_objX = subFrame0ObjInfo(sub0ObjId).objX;
            matchInfo(cnt_matchObj).subFrame0_objY = subFrame0ObjInfo(sub0ObjId).objY;
            matchInfo(cnt_matchObj).subFrame0_objZ = subFrame0ObjInfo(sub0ObjId).objZ;
            [errorDis, sub1ObjId] = min(errorDis_list);
            matchInfo(cnt_matchObj).errorDis = errorDis;
            matchInfo(cnt_matchObj).subFrame1ObjId = sub1ObjId; 
            matchInfo(cnt_matchObj).subFrame1_xyz = subFrame1ObjInfo(sub1ObjId).xyz;
            matchInfo(cnt_matchObj).subFrame1_objX = subFrame1ObjInfo(sub1ObjId).objX;
            matchInfo(cnt_matchObj).subFrame1_objY = subFrame1ObjInfo(sub1ObjId).objY;
            matchInfo(cnt_matchObj).subFrame1_objZ = subFrame1ObjInfo(sub1ObjId).objZ;            
        end
    end
    results.matchInfo = matchInfo;
    
    % 匹配效果展示
    figure();
    cnt_object = 0;
    legend_list = {};
    subFrame0_matchedId = [];
    subFrame1_matchedId = [];
    for i = 1:length(matchInfo)
        cnt_object = cnt_object + 1;
        legend_list{cnt_object} = ['matchedObj', num2str(i)];
        subFrame0_matchedId = [subFrame0_matchedId, matchInfo(i).subFrame0ObjId];
        subFrame1_matchedId = [subFrame1_matchedId, matchInfo(i).subFrame1ObjId];
        xyz = [matchInfo(i).subFrame0_xyz; matchInfo(i).subFrame1_xyz];
        scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 10, 'filled');
        hold on;        
    end
    results.subFrame0_matchedId = subFrame0_matchedId;
    results.subFrame1_matchedId = subFrame1_matchedId;
    
    % 查询未匹配到的object
    for i = 1: length(subFrame0ObjInfo)
        if ~ismember(i, subFrame0_matchedId)
            cnt_object = cnt_object + 1;
            legend_list{cnt_object} = ['subFrame0NotMatchedObj', num2str(i)];
            xyz = subFrame0ObjInfo(i).xyz;
            scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 10, 'o');
            hold on;
        end
    end
    for i = 1: length(subFrame1ObjInfo)
        if ~ismember(i, subFrame1_matchedId)
            cnt_object = cnt_object + 1;
            legend_list{cnt_object} = ['subFrame1NotMatchedObj', num2str(i)];
            xyz = subFrame1ObjInfo(i).xyz;
            scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 10, 'o');
            hold on;
        end
    end
    hold off;
    set(gca, 'Color', [0.8,0.8,0.8]);
    axis(gca, 'equal');
    title('匹配结果');
    legend(legend_list);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    
    
end