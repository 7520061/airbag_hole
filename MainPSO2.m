%% MainPSO2
% ベントホールが高さ方向に散らばるように分布させる
% 様々なつぶれ方で評価

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
% -----------------------------------
% Constraint
ini.P_limit  = 40 * 10^3; % Maximum airbag inner pressure [Pa]
ini.lf_limit = 7.0; % Maximum load factor [G]
% -----------------------------------
% mode
% 1 ... 上半分を6セクションに分けて各セクションに1つずつ穴を配置　→ 下半分に対称移動,様々なつぶれ方で評価
% 2 ... 上半分で上下対称になるように配置 → 下半分に対称移動，上下からでも中腹からでも同じ結果になる
ini.calc_mode = 2;
% -----------------------------------------------

input = Loadmodel(ini);


%% つぶれ方定義
input.vent_close_method1 = '順次閉じる';
input.vent_close_method2 = '中腹から順次閉じる';

%% 結果保存用ファイル作成
% 現在の日付と時間を取得
datetime_now = datetime('now');
date_str = datestr(datetime_now, 'yyyymmdd_HHMMSS');

%front or rear
if ~ini.pos_select
mkdir('PSO2_result/front/',date_str)
    DIRECTORY = append('PSO2_result/front/',date_str);
else
mkdir('PSO2_result/rear/',date_str)
    DIRECTORY = append('PSO2_result/rear/',date_str);
end

input.DIRECTORY = DIRECTORY;

if ~flag
    status = mkdir(DIRECTORY);
    % delete('trajs/*')
    if ~isempty(dir([DIRECTORY '/iter*']))
        rmdir([DIRECTORY '/iter*'],'s');
    end
end
mkdir([DIRECTORY '/iter=0'])

%% PSO
%PSO option
input.method = 'PSO';
nGen = 40; % 100; % Number of iteration or generation
nPop = 200; % 100; % Number of population
input.nGen = nGen;
input.nPop = nPop;

options = optimoptions('particleswarm');
options.Display           = 'iter';  % Information is displayed at each iteration
% options.FunctionTolerance = 0;  % Stopping tolerance
options.MaxIterations     = input.nGen;  % Number of generations
options.MaxStallIterations = 15;
options.OutputFcn         = @mypsofun;   %
options.SwarmSize         = input.nPop;  % Number of populations
options.PlotFcn           = 'pswplotbestf'; % Plots the best objective function value against iterations
options.UseParallel       = true;  % Parallel evaluation option

tic

switch input.calc_mode
    case 1
        input.vent_num = 6; % 設計変数で定義しているベントホールの数
        nVar = 2 * input.vent_num; % 設計変数の数

        % Vent Hole Diamiter [m]
        %                a      b      c      d      e      f
        l_bound_D = [0.000; 0.000; 0.000; 0.000; 0.000; 0.000]; % zeros(nVar,1); % Lower bound of design variable
        u_bound_D = [0.070; 0.070; 0.070; 0.070; 0.070; 0.070]; % ones(nVar,1);  % Upper bound of design variable

        % Vent Hole Location [m]
        % 6セクションに分けて各セクションに1つずつ穴を配置
        max_L = 0.5; % ベントホールの中心が取りうる範囲[m]
        vent_hole_L = max_L / input.vent_num; %各セクションの範囲[m]

        l_bound_L = (vent_hole_L * linspace(0,input.vent_num - 1,input.vent_num))'; % Lower bound of design variable
        u_bound_L = (vent_hole_L * linspace(1,input.vent_num,input.vent_num))'; % Upper bound of design variable


        l_bound = [l_bound_D;l_bound_L];
        u_bound = [u_bound_D;u_bound_L];

        ObjFunc = @ObjFunc_PSO2; % 目的関数の計算をする関数
        FigFunc = @fig_function; % 結果表示する関数

    case 2
        input.vent_num = 3; % 設計変数で定義しているベントホールの数
        nVar = 2 * input.vent_num; % 設計変数の数

        % Vent Hole Diamiter [m]
        %                a      b      c      
        l_bound_D = [0.000; 0.000; 0.000]; % zeros(nVar,1); % Lower bound of design variable
        u_bound_D = [0.070; 0.070; 0.070]; % ones(nVar,1);  % Upper bound of design variable

        % Vent Hole Location [m]
        max_L = input.bag_Hs; % ベントホールの中心が取りうる範囲[m]
        vent_hole_L = max_L / 2; 

        l_bound_L = vent_hole_L .* ones(input.vent_num,1); % Lower bound of design variable
        u_bound_L = 2 * vent_hole_L .* ones(input.vent_num,1); % Upper bound of design variable

        l_bound = [l_bound_D;l_bound_L];
        u_bound = [u_bound_D;u_bound_L];

        ObjFunc = @ObjFunc_PSO2; % 目的関数の計算をする関数
        FigFunc = @fig_function; % 結果表示する関数

end

if flag
    load('_Rear Airbag_hd偶数_千鳥_25mm\20230510_132455\opt.mat','x_opt');
    options.InitialSwarmMatrix = x_opt;
end

x_opt = particleswarm(@(TBD)ObjFunc(TBD, input),nVar,l_bound,u_bound,options);


%% best
x_opt = Para_sort(x_opt);
[obj_opt,landing_opt,pnlt_opt] = ObjFunc(x_opt, input);

%% Plot
FigFunc(landing_opt.t_result, landing_opt.result, input);
fig_result = gca;
saveas(fig_result,[DIRECTORY '/Result_opt.fig']);
saveas(fig_result,[DIRECTORY '/Result_opt.png']);

vent_hole_plot_Zigzag(input, x_opt);
% vent_hole_plot_Straight(input,x_opt);
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
xlabel('Iteration','FontSize',16,'Fontname','Times New Roman');
ylabel('\itJ_{\rmEA} ','interpreter','tex','FontSize',16);
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
delete(gcp('nocreate'))

toc
