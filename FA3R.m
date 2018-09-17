% Fast Analytical 3D Registration
% Submitted to IEEE Transactions on Robotics
% Authors: Jin Wu, Ming Liu et al.
% Copytight (c) 2018

function [C, T, time, loss, expn] = FA3R( MM, r_base, b_base, compute_expn, compute_metric )

    nn = length(b_base(1, :));
    
    mean_X = zeros(3, 1);
    mean_Y = zeros(3, 1);
    
    for i = 1 : nn
        mean_X = mean_X + 1 / nn * r_base(:, i);
        mean_Y = mean_Y + 1 / nn * b_base(:, i);
    end
    
    if(MM == 0)
        MM = zeros(3, 3);
   
        for i = 1 : length(b_base(1, :))
            MM = MM + 1 / nn * (b_base(:, i) - mean_Y) * (r_base(:, i) - mean_X)';
        end
    end
    
    if(compute_expn)
        HX = [];
        HY = [];
        HZ = [];
    end
    
    tic;
    
    hx = MM(:, 1)';
    hy = MM(:, 2)';
    hz = MM(:, 3)';
    
    for iter = 1 : 10
        k = 2.0 / (hx(1) * hx(1) + hx(2) * hx(2) + hx(3) * hx(3) + ...
                   hy(1) * hy(1) + hy(2) * hy(2) + hy(3) * hy(3) + ...
                   hz(1) * hz(1) + hz(2) * hz(2) + hz(3) * hz(3) + 1);
             
        hx_ = hx;  hy_ = hy; hz_ = hz;
        hz = (hz_ + cross(hx_, hy_)) * k;
        hy = (hy_ + cross(hz_, hx_)) * k;
        hx = (hx_ + cross(hy_, hz_)) * k;
        
        if(compute_metric)
            HX = [HX hx'];
            HY = [HY hy'];
            HZ = [HZ hz'];
        end
    end
    
    C = [hx; hy; hz];
    T = mean_X - C * mean_Y;
    time = toc;
    
    
    if(compute_expn)
        expn = zeros(iter, 1);
        HX = HX';
        HY = HY';
        HZ = HZ';
        for i = 1 : iter
            C = [HX(i, :); HY(i, :); HZ(i, :)];
            expn(i) = 0;
            for j = 1 : nn
                expn(i) = expn(i) + 1 / nn * norm(r_base(:, j) - C * b_base(:, j) - T)^2;
            end
        end
    end

    if(compute_metric)
        loss = 0;
        for i = 1 : nn
            loss = loss + 1 / nn * norm(r_base(:, i) - C * b_base(:, i) - T)^2;
        end
    end
end

