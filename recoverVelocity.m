function  v_recover = recoverVelocity(v, time, velocityBin_list)
%     if time > 0
%         % 频谱向+方向搬移
%         
%     elseif time < 0
%         % 频谱向-方向搬移
%         
%     else
%         v_recover = v;
%     end
    
    v_recover = v + time * (velocityBin_list(end) * 2);
        

end