%% ベントホール半径が縁-中心間距離より"小さい(r<=h0)" 範囲
% "下"から閉じる

function ratio = vent_block_ratio3(hz,h0,D,method)

hz = hz/2;
h0 = h0/2;
r = D/2;
S0 = pi*r^2;

switch method
    case '1'
        if ~(hz >= h0-r && hz <= h0+r)
            error('使い方が違う！！')
        end

        dh = abs(hz-h0);

        if hz <= h0
            theta   = 2*acos(dh/r);            % 中心角θ [rad.]
            S_block = 1/2*r^2*(theta - sin(theta)); % 塞がっている面積 [m2]
            S       = S0 - S_block;
        else
            theta   = 2*acos(dh/r);
            S = 1/2*r^2*(theta - sin(theta)); % 空いている面積 [m2]
        end

        ratio = S / S0;


    case '2'
        dh = hz;
        theta = 2*acos(dh/r);            % 中心角θ [rad.]
        S = 1/2*r^2*(theta - sin(theta)); % 空いている面積 [m2]
        S = 2*S;

        ratio = S / S0;

end
end