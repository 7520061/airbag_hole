close all
clear all

%% PATH
% addpath(genpath('util'))
dir0 = pwd;
addpath(fullfile(dir0,'/util'));

global DIRECTORY
rng(1); % 乱数生成
flag = 0; % 0:初期解ランダム, 1:初期解の1つに現状の最適解を採用(opt.matがないとエラー)
% -----------------------------------
% 初期条件
ini.pos_select = 1;     % 前後選択：0 = front, 1 = rear
ini.area_select = 1;    % 場所選択：0 = エスレンジ, 1 = メッペン, 2 = 生花苗沼
ini.T0_cel = 15;       % バッグ内初期温度　セルシウス[c.]
ini.margin = 0;         % ガスマージン[g]
ini.p_bag = 0 * 1000;   % バッグ初期内圧[PaG]
ini.close_H = 0.3;       % ベントホール閉じる高さ[m]
% ini.the = -10;       % 初期ピッチ角 [deg.]
% ini.phi = 0;         % 初期ロール角 [deg.]
% ini.AoA = 90;        % 初期迎角 [deg.]
% ini.beta = 0;        % 初期横滑り角 [deg.]
% ini.UE  = 0;         % 地球固定座標X,風速 [m/s]
% ini.VE  = 0;         % 地球固定座標Y,風速 [m/s]
% ini.WE  = 0;         % 地球固定座標Z,風速 [m/s]
% -----------------------------------
% Constraint
ini.P_limit  = 40 * 10^3; % Maximum airbag inner pressure [Pa]
ini.lf_limit = 7.0; % Maximum load factor [G]
% -----------------------------------
ini.calc_mode = 0;

input = Loadmodel(ini);

input.sim_mode = 'sim';

input.vent_close_method = '中腹から順次閉じる_面積比から閉じ量求める';
input.Vc0 = -5;

% 現在の日付と時間を取得
datetime_now = datetime('now');
date_str = datestr(datetime_now, 'yyyymmdd_HHMMSS');

if ~ini.pos_select
mkdir('Simu_result/front/',date_str)
    DIRECTORY = append('Simu_result/front/',date_str);
else
mkdir('Simu_result/rear/',date_str)
    DIRECTORY = append('Simu_result/rear/',date_str);
end

input.DIRECTORY = DIRECTORY;

vent_num = 12;
input.vent_num = 12;
nVar = 12; % input.nVar; % Number of design variable
input.nVar = nVar;

% 設計変数読み込み
% load(append(dir0,'\_Rear Airbag_hd偶数_千鳥_25mm\20230510_132455\opt.mat'),'x_opt');

%       [   D1    D2    D3    D4    D5    D6    h7    h8    h9   h10   h11   h12] 
% x_opt = [0.046 0.052 0.024 0.062 0.056 0.054 0.520 0.508 0.492 0.466 0.400 0.324];
% x_opt = [0.046 0.052 0.024 0.062 0.056 0.054 0.506 0.500 0.492 0.466 0.400 0.324];
% x_opt = [0.033 0.037 0.024 0.062 0.056 0.054 0.506 0.500 0.492 0.466 0.400 0.324];

x_opt = [];
x_opt = [0.042 0.056 0.006 0.070 0.070 0 0.448 0.436 0.426 0.422 0.310 0.0240]; %中腹から閉じるで最適化

[obj_opt,landing_opt,pnlt_opt] = ObjFunc_sim(x_opt, input);

% Plot
fig_function(landing_opt.t_result, landing_opt.result, input);
fig_result = gca;
saveas(fig_result,[DIRECTORY '/Result_opt.fig']);
saveas(fig_result,[DIRECTORY '/Result_opt.png']);

% vent_hole_plot_Zigzag(input, x_opt);
vent_hole_plot_Straight(input,x_opt);
vent_hole_result = gca;
saveas(gca,[DIRECTORY '/vent_hole.fig']);
saveas(gca,[DIRECTORY '/vent_hole.png']);

figure(1);
fig_fv = gca;
fig_fv.Children.MarkerSize = 12;
hold on; box on;
set(gca,'fontname','Times New Roman','fontsize',14);
set(gcf,'Color',[1 1 1]);
set(gca,'xlim',[0 inf],'FontSize',14,'Fontname','Times New Roman')
saveas(fig_fv,[DIRECTORY '/ObjFncV.fig']);
saveas(fig_fv,[DIRECTORY '/ObjFncV.png']);

if obj_opt <= 80
        set(gca,'ylim',[0 80])
    saveas(fig_fv,[DIRECTORY '/ObjFncV_zoom.fig']);
    saveas(fig_fv,[DIRECTORY '/ObjFncV_zoom.png']);
elseif obj_opt <= 1500
    set(gca,'ylim',[0 1500])
    saveas(fig_fv,[DIRECTORY '/ObjFncV_zoom.fig']);
    saveas(fig_fv,[DIRECTORY '/ObjFncV_zoom.png']);
elseif obj_opt <= 3000
    set(gca,'ylim',[0 3000])
    saveas(fig_fv,[DIRECTORY '/ObjFncV_zoom.fig']);
    saveas(fig_fv,[DIRECTORY '/ObjFncV_zoom.png']);
end

% save opt
save([DIRECTORY '/opt'],'obj_opt','landing_opt','pnlt_opt','input','x_opt');
