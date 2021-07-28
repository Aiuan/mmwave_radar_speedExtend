% 从/input/xxx_param.m文件中获取雷达参数
function par = getPara(paraFile, paraName)
    run(paraFile);
    par = eval(paraName);
end