function ratio = vent_block_ratio(hz,vent_H,vent_D,method)

% 開いている面積の割合を算出

vH = vent_H;
vD = vent_D;
r = vent_D/2;
S0 = pi*r^2;

% method
% 1 ... 一方向に進む
% 2 ... 中心から上下方向に進む

switch method
    case '1'
        if ~(hz >= vH-r && hz <= vH+r)
            error('使い方が違う！！')
        end

        if vH > hz
            theta = 2*acos((vH-hz)/vD);
            S_block = 1/2 * r^2 * theta - 1/2 * r^2 * sin(theta);

        else
            theta = 2*acos((hz-vH)/vD);
            S = 1/2 * r^2 * theta - 1/2 * r^2 * sin(theta);
            S_block = S0-S;

        end

    case '2'
        if ~(vH == 0)
            error('使い方が違う！！')
        end

        theta = 2*acos(hz/vD);

        S_block = 2*(1/2 * r^2 * (pi-theta) + 1/2 * r^2 * sin(theta));

end


ratio = 1 - S_block/S0;

