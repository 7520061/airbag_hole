function input = Loadmodel(ini)

% 時間の設定
input.t0 = 0;
input.te_plot = 0.5;                 % シミュレーション時間 %%
input.dt_plot = 0.2 * input.te_plot; % プロットしていく刻み時間
input.tl_plot = input.te_plot;       % グラフのx軸はどこまで書きますか~？
input.dt = 0.001;                    % シミュレーション計算の刻み時間
input.t_end = input.te_plot;         % 何秒までシミュレーション回しますか~？

input.atmos = loadAtmos_airbag; % Atmosphere model U.S Standard Atomosphere 1976

%% Constraint
input.P_limit  = ini.P_limit;  % Maximum airbag inner pressure [Pa]
input.lf_limit = ini.lf_limit; % Maximum load factor [G]

%% Static Variable
temp = csvread('StaticVariable.csv',0,1);

% #015機体諸元
input.body_l = temp(1);
input.body_S = temp(2);
% input.body_cmac = temp(3);
% input.body_bsp = temp(4);
input.body_Cd = temp(5);

% input.body_ElevonL = 0; % エレボンLの角度 [deg.]
% input.body_ElevonR = 0; % エレボンRの角度 [deg.]
% input.body_Rudder  = 0; % ラダー舵角 [deg.]

% input.WIRES_Aero = AerodynamicCoefficient_Init(); % 空力データ，KHI修正版

% 重心
input.CGx = temp(6);
input.CGmm = input.body_l * 10^3 * input.CGx / 100;
% input.CGy = temp(7);
% input.CGz = temp(8);
% input.CGpos  = [input.body_l*input.CGx/100;
%                 input.body_l*input.CGy/100;
%                 input.body_l*input.CGz/100];    % 重心位置[m] [x y z]

% 慣性モーメント
% input.Ixx = temp(9);
% input.Iyy = temp(10);
% input.Izz = temp(11);
% input.Ixz = temp(12);
% input.Ixy = 0;
% input.Iyz = 0;
% input.I = [ input.Ixx, -input.Ixy, -input.Ixz;
%            -input.Ixy,  input.Iyy, -input.Iyz;
%            -input.Ixz, -input.Iyz,  input.Izz];

% エアバッグ搭載位置
input.front_pos = temp(13);
input.rear_pos  = temp(14);
input.front_ratio = (input.rear_pos - input.CGmm)/(input.rear_pos - input.front_pos);
input.rear_ratio  = (input.CGmm - input.front_pos)/(input.rear_pos - input.front_pos)/2; % リアは2個あるから÷2

% 機体重量
input.m_lf        = temp(15);
input.m_LOX       = temp(16);
input.m_LNG       = temp(17);
input.mp_residual = temp(18);
input.m_margin    = temp(19);

input.m_recovery = input.m_lf - (input.m_LOX + input.m_LNG) * (1 - input.mp_residual / 100) + input.m_margin;
% input.m_recovery = 1050;
input.m0_front = input.m_recovery * input.front_ratio;
input.m0_rear  = input.m_recovery * input.rear_ratio;

% 充填気体諸元
input.R0 = 8.3144598;                         % 気体定数 [J/(molK)]
input.atomicm = 28.966;                       % 乾燥空気平均分子量 [g/mol]
input.R = input.R0 * 1000 / input.atomicm;    % 気体定数　N2 [J/(kgK)]
input.GAM = 1.4;                              % 比熱比　N2[-]
input.cv = input.R / (input.GAM - 1);         % 定積比熱

% input.drogue_mass = 1.9 + 0.507 + 1.01 + 0.222 + 0.579 + 0.224;

% メインシュート諸元
input.main_S = 175;  % メインシュート投影面積[m2]

input.g_1000 = calc1DInterp(input.atmos.htab, input.atmos.gtab, 1000, input.atmos.gpp);        % 重力加速度@1000m[m/s2]
input.g_1500 = calc1DInterp(input.atmos.htab, input.atmos.gtab, 1500, input.atmos.gpp);        % 重力加速度@1500m[m/s2]

input.rho_1000 = calc1DInterp(input.atmos.htab, input.atmos.Rhotab, 1000, input.atmos.Rhopp);  % 大気密度@1000m[kg/m3]
input.rho_1500 = calc1DInterp(input.atmos.htab, input.atmos.Rhotab, 1500, input.atmos.Rhopp);  % 大気密度@1500m[kg/m3]

cd_1 = 2 * 800 * input.g_1000/(input.rho_1000 * input.main_S * 7.1^2);
cd_2 = 2 * 800 * input.g_1500/(input.rho_1500 * input.main_S * 7.3^2);
cd_3 = 2 * 990 * input.g_1000/(input.rho_1000 * input.main_S * 7.9^2);
cd_4 = 2 * 990 * input.g_1500/(input.rho_1500 * input.main_S * 8.1^2);

input.main_cd = (cd_1 + cd_2 + cd_3 + cd_4)/4;

% エアバッグ諸元
input.bag_Hs  = temp(20);             % 側面高さ（ベントホール縫製面[m])
input.bag_H0  = 0.25 + input.bag_Hs;  % 展開時バッグ高さ [m]
% input.bag_Hs = 0.552+0.06;
% input.bag_H0 = 0.802;
input.bag_H   = input.bag_H0;         % リバウンド更新用 [m]

%% Initial state
input.Flag = 0;

% Environment
% input.UE   = ini.UE;   % 地球固定座標X,風速
% input.VE   = ini.VE;   % 地球固定座標Y,風速
% input.WE   = ini.WE;   % 地球固定座標Z,風速

if ini.area_select == 0
    input.h_sl = 341; % 海抜高度[ｍ] (エスレンジ)
elseif ini.area_select == 1
    input.h_sl = 23;  % 海抜高度[ｍ] (メッペン)
elseif ini.area_select == 2
    input.h_sl = 2;   % 海抜高度[ｍ] (生花苗沼)
end

input.h0 = input.bag_H0 + 0.001;% + 2; % Initial Altitude @ground level [m]
input.h0_sl = input.bag_H0 + input.h_sl; % Initial Altitude @sea level [m]

input.g = calc1DInterp(input.atmos.htab, input.atmos.gtab, input.h0_sl, input.atmos.gpp);        % 重力加速度[m/s2]
input.p = calc1DInterp(input.atmos.htab, input.atmos.Prestab, input.h0_sl, input.atmos.Presspp); % 外気圧 [PaA] (U.S. STANDARD ATMOSPHRE,1976 線形補完)
input.rho = calc1DInterp(input.atmos.htab, input.atmos.Rhotab, input.h0_sl, input.atmos.Rhopp);  % 大気密度[kg//m3]

% エアバッグ直径計算
input.loadmax_front = input.m0_front * input.g * input.lf_limit;    % エアバッグ一個当たりにかかる荷重[N]
input.loadmax_rear  = input.m0_rear * input.g * input.lf_limit;     % エアバッグ一個当たりにかかる荷重[N]

input.S_correction_front = (0.93^2 * pi / 4) / 0.664;
input.Sreq_front = input.S_correction_front * (input.loadmax_front / input.P_limit);
input.Sreq_rear  = input.loadmax_rear / input.P_limit;

input.bag_D_front = 2 * sqrt(input.Sreq_front / pi);
input.bag_D_rear  = 2 * sqrt(input.Sreq_rear  / pi);

input.bag_D_front = round(input.bag_D_front,3);
input.bag_D_rear  = round(input.bag_D_rear, 3);
input.pos_select = ini.pos_select;
if ~ini.pos_select
    input.m0 = input.m0_front;    % エアバッグ一個当たりにかかる重量[kg]
    input.bag_s_correction = input.S_correction_front; % 接触面積修正係数 (0.93mのときのCADモデル実接触面積との比より)
    input.bag_D   = input.bag_D_front;
else
    input.m0 = input.m0_rear;     % エアバッグ一個当たりにかかる重量[kg]
    input.bag_s_correction = 1;
    input.bag_D   = input.bag_D_rear;
end

input.bag_wa  = 1.608;                % 楕円表面積係数
input.bag_wet_area_d = 4 * pi * ((2 * (input.bag_D / 2 * (input.bag_H0 - input.bag_Hs) / 2)^input.bag_wa + (input.bag_D / 2)^(2 * input.bag_wa)) / 3)^(1 / input.bag_wa);               % ドーム形状の表面積[m2]
input.bag_wet_area_s = input.bag_D * pi * input.bag_Hs;                            % 側面の表面積[m2]
input.bag_V_reef = (input.bag_wet_area_d + input.bag_wet_area_s) / 0.85 * 0.00128; % 畳んだ時の体積[m3]
input.bag_S = 0; % 接触面積 [m2]
input.bag_V = 4 / 3 * pi * input.bag_D^2 / 4 * (input.bag_H0 - input.bag_Hs) / 2 + input.bag_D^2 * pi / 4 * input.bag_Hs; % バッグ体積(ガス量過不足なし)[m3]
input.bag_V_nominal = input.bag_V ; % ノミナル初期バッグ体積(ガス量過不足なし)[m3]

% エアバッグ内部ガス
input.T0_cel = ini.T0_cel;
input.gas_margin = ini.margin;

input.T0     = 273.15 + input.T0_cel;   % バッグ内初期温度[K]
input.p0     = input.p + ini.p_bag;     % 初期バッグ圧力[PaA]
input.mol    = input.p0 * input.bag_V_nominal /(input.R0 * input.T0); % 充填ガスモル量[mol]
input.gas_m0 = (input.mol * input.atomicm - input.gas_margin) * 10^-3;

input.bag_V0 = input.mol * input.R0 * input.T0 / input.p0; % 初期バッグ体積　[m3], ゲージ圧0, ガス量に合わせて体積は小さくなる計算
input.bag_rho0 = input.gas_m0 / input.bag_V0;  % 初期密度[kg/m3]
input.acc = -input.g;
input.mae = 0;                                 % 初期排気マッハ数[-]
input.ve = 0;                                  % 初期排気速度[m/s]
input.d_mass = 0;

% input.AoA0  = ini.AoA * pi/180;    % 初期迎角 [rad]
% input.beta0 = ini.beta * pi/180;  % 初期横滑り角 [rad]
% 
% input.the0 = ini.the * pi/180;   % 初期ピッチ角 [rad]
% input.phi0 = ini.phi * pi/180;   % 初期ロール角 [rad]
% input.psi0 = 0  *pi/180;         % 初期ヨー角 [rad]
% 
% input.p0   = 0  *pi/180;      % 初期ロール角速度 [rad/s]
% input.q0   = 0  *pi/180;      % 初期ピッチ角速度 [rad/s]
% input.r0   = 0  *pi/180;      % 初期ヨー角速度 [rad/s]

input.Vc0 = -sqrt((2 * input.m_recovery * input.g)/(input.rho * input.main_cd * input.main_S));  % 降下速度[m/s] (機体抗力無視)
% input.Vc0 = -6;
% input.Vc0 = -sqrt((2 * input.m_recovery * input.g)/(input.rho * (input.main_cd * input.main_S + input.body_S * input.body_Cd)));  % 降下速度[m/s] (機体迎角40°抗力込み)
% %% 機体軸系に変換
% input.u0 = [                                                  cos(input.the0)*cos(input.psi0),                                                   cos(input.the0)*sin(input.psi0),                -sin(input.the0)]*[0;0;input.Vc0];
% input.v0 = [sin(input.the0)*cos(input.psi0)*sin(input.phi0) - cos(input.phi0)*sin(input.psi0), sin(input.the0)*sin(input.phi0)*sin(input.psi0) - cos(input.phi0)*cos(input.psi0), cos(input.the0)*sin(input.phi0)]*[0;0;input.Vc0];
% input.w0 = [sin(input.phi0)*sin(input.psi0) + sin(input.the0)*cos(input.phi0)*cos(input.psi0), sin(input.the0)*cos(input.phi0)*sin(input.psi0) - cos(input.psi0)*sin(input.phi0), cos(input.the0)*cos(input.phi0)]*[0;0;input.Vc0];
% 
% input.DR0  = 0;	% 初期ダウンレンジ [m]
% input.CR0  = 0;	% 初期クロスレンジ [m]

% Number of Vent Hole
input.an0 = 2; input.bn0 = 2; input.cn0 = 2; input.dn0 = 2; input.en0 = 2; input.fn0 = 2;
input.n0_ind = 4;

% input.an0 = 4; input.bn0 = 4; input.cn0 = 4; input.dn0 = 4; input.en0 = 4; input.fn0 = 4;

input.vent_s = 0;

input.break_pressure = input.p / 10^3 + 10;    % ベントホール開口圧力[kPa]           

% オリフィス流量係数(参考：https://yuruyuru-plantengineer.com/orifice-flow-coefficient/)
input.vent_Ca = 0.96; % 縮流係数
input.vent_Cv = 0.64; % 速度係数
input.vent_Cd = 0.62; % 流量係数

% バック高さ更新フラグ [-]
input.Vflag = 0;% + input.gn0 + input.hn0 + input.in0 + input.jn0 + input.kn0 + input.ln0;

% エアバッグがぺったんこになったら強制終了（(；゜Д゜)
input.bag_dh = 0.001;

% 実験値からの断面積の値
input.ex_bag_S = 0.512190218;

% 減衰係数
% 減衰項　d_c(1)*速度 + d_c(2)
input.d_c = [-72.946,33.811];
% input.d_c = [0,0];

% mode
input.calc_mode = ini.calc_mode;
