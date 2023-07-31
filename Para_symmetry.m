function para_calc = Para_symmetry(para,input)
% 上半分で対称になるように配置

des_num = length(para); %設計変数の数
Hs_half = input.bag_Hs / 2; %エアバッグ上半分の高さ

vent_H = para(des_num/2+1:end);
vent_H = [vent_H,flip(Hs_half-vent_H)];

vent_D = para(1:des_num/2);
vent_D = [vent_D,flip(vent_D)];

para_calc = [vent_D,vent_H];
end