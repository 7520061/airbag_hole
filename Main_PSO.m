close all
clear all

%% PATH
% addpath(genpath('util'))
dir0 = pwd;
addpath(fullfile(dir0,'/util'));

global DIRECTORY
rng(1); % ��������
flag = 0; % 0:�����������_��, 1:��������1�Ɍ���̍œK�����̗p(opt.mat���Ȃ��ƃG���[)
% -----------------------------------
% ��������
ini.pos_select = 1;     % �O��I���F0 = front, 1 = rear
ini.area_select = 1;    % �ꏊ�I���F0 = �G�X�����W, 1 = ���b�y��, 2 = ���ԕc��
ini.T0_cel = 15;       % �o�b�O���������x�@�Z���V�E�X[c.]
ini.margin = 0;         % �K�X�}�[�W��[g]
ini.p_bag = 0 * 1000;   % �o�b�O��������[PaG]
ini.close_H = 0.3;       % �x���g�z�[�����鍂��[m]
% ini.the = -10;       % �����s�b�`�p [deg.]
% ini.phi = 0;         % �������[���p [deg.]
% ini.AoA = 90;        % �����}�p [deg.]
% ini.beta = 0;        % ����������p [deg.]
% ini.UE  = 0;         % �n���Œ���WX,���� [m/s]
% ini.VE  = 0;         % �n���Œ���WY,���� [m/s]
% ini.WE  = 0;         % �n���Œ���WZ,���� [m/s]
% -----------------------------------
% Constraint
ini.P_limit  = 40 * 10^3; % Maximum airbag inner pressure [Pa]
ini.lf_limit = 7.0; % Maximum load factor [G]
% -----------------------------------

input = Loadmodel(ini);

% input.h0 = 2.7;
% input.m0 = 252.025;
% input.Vc0 = 0;

% ���݂̓��t�Ǝ��Ԃ��擾
datetime_now = datetime('now');
date_str = datestr(datetime_now, 'yyyymmdd_HHMMSS');


if ~ini.pos_select
mkdir('PSO_result/front/',date_str)
    DIRECTORY = append('PSO_result/front/',date_str);
else
mkdir('PSO_result/rear/',date_str)
    DIRECTORY = append('PSO_result/rear/',date_str);
end

if ~flag
    status = mkdir(DIRECTORY);
    % delete('trajs/*')
    if ~isempty(dir([DIRECTORY '/iter*']))
        rmdir([DIRECTORY '/iter*'],'s');
    end
end
mkdir([DIRECTORY '/iter=0'])

input.method = 'PSO';
nGen = 40; % 100; % Number of iteration or generation
nPop = 200; % 100; % Number of population
input.nGen = nGen;
input.nPop = nPop;

tic

% Vent Hole Diamiter [m]
%                   a      b      c      d      e      f
l_bound = [         0.020; 0.020; 0.020; 0.020; 0.020; 0.020;]; % zeros(nVar,1); % Lower bound of design variable
u_bound = [         0.070; 0.070; 0.070; 0.070; 0.070; 0.070;]; % ones(nVar,1);  % Upper bound of design variable

% Vent Hole Location [m]
%                   a      b      c      d      e      f
l_bound = [l_bound; 0.000; 0.000; 0.000; 0.000; 0.000; 0.000;]; % zeros(nVar,1); % Lower bound of design variable
u_bound = [u_bound; 0.500; 0.500; 0.500; 0.500; 0.500; 0.500;]; % ones(nVar,1);  % Upper bound of design variable

% % ��ɊJ���Ă��錊
% l_bound =[l_bound; 0.000];
% u_bound =[u_bound; 0.050];

input.l_bound = l_bound;
input.u_bound = u_bound;

vent_num = 12;
input.vent_num = 12;
nVar = 12; % input.nVar; % Number of design variable
input.nVar = nVar;

input.DIRECTORY = DIRECTORY;

options = optimoptions('particleswarm');
options.Display           = 'iter';  % Information is displayed at each iteration
% options.FunctionTolerance = 0;  % Stopping tolerance
options.MaxIterations     = input.nGen;  % Number of generations
options.MaxStallIterations = 15;
options.OutputFcn         = @mypsofun;   %
options.SwarmSize         = input.nPop;  % Number of populations
options.PlotFcn           = 'pswplotbestf'; % Plots the best objective function value against iterations
options.UseParallel       = true;  % Parallel evaluation option

if flag
    load('_Rear Airbag_hd����_�璹_25mm\20230510_132455\opt.mat','x_opt');
    options.InitialSwarmMatrix = x_opt;
end

x_opt = particleswarm(@(TBD)ObjFunc_PSO(TBD, input),nVar,l_bound,u_bound,options);

% Best
for j = 1 : 1 : vent_num/2
    x_opt(j) = x_opt(j) * 1000;
    x_opt(j) = round(x_opt(j),TieBreaker="even");
    if rem(x_opt(j), 2) == 1
        x_opt(j) = x_opt(j) + 1;
    end
    x_opt(j) = x_opt(j) / 1000;
    if x_opt(j) < 0.006
        x_opt(j) = 0.000;
    end
end

for j = vent_num/2+1 : 1 : vent_num
    x_opt(j) = x_opt(j) * 1000;
    x_opt(j) = round(x_opt(j),TieBreaker="even");
    if rem(x_opt(j), 2) == 1
        x_opt(j) = x_opt(j) + 1;
    end
    x_opt(j) = x_opt(j) / 1000;
end

vent_D = (x_opt(1:vent_num/2))';
vent_H = (x_opt(vent_num/2+1:vent_num))';

[vent_H_sorted,I] = sort(vent_H,'descend');
vent_D_sorted = vent_D(I);

for j = 1 : 1 : vent_num/2
    x_opt(j) = vent_D_sorted(j);
end

for j = vent_num/2+1 : 1 : vent_num
    x_opt(j) = vent_H_sorted(j-6);
end

[obj_opt,landing_opt,pnlt_opt] = ObjFunc_PSO(x_opt, input);

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
