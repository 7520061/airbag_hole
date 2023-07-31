%% Main Straight Line
%
close all
clear all
dbstop if error

%% PATH
% addpath(genpath('util'))
dir0 = pwd;
addpath(fullfile(dir0,'/util'));

global DIRECTORY

% -----------------------------------
% 初期条件
ini.pos_select = 1;     % 前後選択：0 = front, 1 = rear
ini.area_select = 1;    % 場所選択：0 = エスレンジ, 1 = メッペン, 2 = 生花苗沼
ini.T0_cel = 15;       % バッグ内初期温度　セルシウス[c.]
ini.margin = 0;         % ガスマージン[g]
ini.p_bag = 0 * 1000;   % バッグ初期内圧[PaG]
% -----------------------------------
% Constraint
ini.P_limit  = 40 * 10^3; % Maximum airbag inner pressure [Pa]
ini.lf_limit = 7.0; % Maximum load factor [G]
% -----------------------------------
% mode
% 1 ... 上半分を6セクションに分けて各セクションに1つずつ穴を配置　→ 下半分に対称移動,様々なつぶれ方で評価
% 2 ... 上半分で上下対称になるように配置 → 下半分に対称移動，上下からでも中腹からでも同じ結果になる
ini.calc_mode = 0;
% -----------------------------------------------

input = Loadmodel(ini);

input.sim_mode = 'SL';

%% つぶれ方定義
input.vent_close_method = '中腹から順次閉じる_面積比から閉じ量求める';
input.vent_close_num = 3; %どの穴から閉じるか　下から何番目かを指定　1以上vent_num以下
% input.close_height = 0.100;

% 現在の日付と時間を取得
datetime_now = datetime('now');
date_str = datestr(datetime_now, 'yyyymmdd_HHMMSS');

if ~ini.pos_select
    res_str = ['SL_result/front/' date_str];
else
    res_str = ['SL_result/rear/' date_str];
end
    mkdir([res_str '/data_ind'])
    DIRECTORY = (res_str);

input.DIRECTORY = DIRECTORY;

input.hole_dif = 0.025; %穴ふちの最小距離
section_range = input.bag_Hs - 0.06; %穴が配置される範囲

% Vent Hole Diameter [m]
%
D_min = 0.00;
D_max = 0.07;
D_span = 0.001;
D = D_min:D_span:D_max;

% 穴ふち距離は最小にして並べる
for j = 1:1:length(D)

    vent_range = section_range - D(j);
    vent_num = floor(vent_range / (D(j)+input.hole_dif)) + 1;
    input.vent_num = vent_num;
    center_dif = D(j)+input.hole_dif;

    if center_dif <= 0.00001 
        continue
    end

    if rem(vent_num,2) == 1
        vent_num_half = (vent_num-1)/2;
        vent_H_end = vent_num_half*center_dif;
    else
        vent_num_half = vent_num/2;
        vent_H_end = (vent_num_half - 0.5)*center_dif;
    end

    vent_H = (-vent_H_end:center_dif:vent_H_end);
    vent_D = D(j)*ones(1,vent_num);
    para = [vent_D,vent_H];

    [~,landing,~] = ObjFunc_sim(para, input);
    result_SL(j,:) = [D(j),landing.Gmax_nominal,landing.Vmin_nominal];

%     if landing.Vmin_nominal >= 0.5 % 最小降下速度
%         break
%     end

end


% %% ベントホール穴位置座標
% input.hole_dif_min = 0.025; %穴ふちの最小距離
% input.section_range = input.bag_Hs - 0.06; %穴が配置される範囲
% vent_num = 11;% Number of vent hole
% input.vent_num = vent_num;
%
% % Vent Hole Diameter [m]
% %
% D_min = 0.000; % zeros(nVar,1); % Lower bound of design variable
% D_max = 0.070; % ones(nVar,1);  % Upper bound of design variable
% D_span = 0.001; %Dの刻み幅
% D = D_min:D_span:D_max;
%
% center_dif = input.section_range / vent_num;
%
%
% % エアバッグ中心を原点とした座標系での穴位置の座標
% vent_abs_H = linspace(-(input.section_range-center_dif) / 2, (input.section_range-center_dif) / 2, vent_num);
%
% %% 穴径を振って性能評価
%
% for j = 1:1:length(D)
%     vent_abs_D = D(j) .* ones(1,vent_num);
%     vent_abs_para = [vent_abs_D,vent_abs_H];
%
%     if input.hole_dif_min > center_dif - D(j) % 穴ふち最小距離
%         break
%     end
%
%     [~,landing,~] = ObjFunc_sim(vent_abs_para, input);
%     result_D_Gmax(j,:) = [j-1,landing.Gmax_nominal];
%
%     if landing.Vmin_nominal >= 1.00 % 最小降下速度
%         break
%     end
%
%     para(j,:) =vent_abs_para;
% end



%% plot
result_D = 1000*result_SL(1:j-1,1);
result_Gmax = result_SL(1:j-1,2);
result_Vmin = result_SL(1:j-1,3);
figure;
hold on
grid on
scatter(result_D,result_Gmax,[],result_Vmin,'filled')
yline(7,'-r','7G');
c = colorbar;
c.Label.String = '最小降下速度';
xlabel('穴直径[mm]')
ylabel('最大荷重倍数[G]')
fig_result = gca;
saveas(fig_result,[DIRECTORY '/Result_D_Gmax.fig']);

% -----------------------------------------------------------------
% for i = 1:1:length(vent_num)
% center_dif = input.section_range / vent_num(i); %穴の中心間距離
%
% % Vent Hole Location [m]
% %
% vent_H_start = -(input.section_range - center_dif)/2;
% vent_H_end = (input.section_range - center_dif)/2;
% vent_H = linspace(vent_H_start,vent_H_end,vent_num(i));
%
% for j = 1:1:length(D)
%     vent_D = D(j) .* ones(1,vent_num(i));
%     para_temp = [vent_D,vent_H];
%
%
% end
%
% end
% ----------------------------------------------------------------

%   target_ind = coeff(coeff(:,1)==i,:);

