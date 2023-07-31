function vent_hole_plot_Straight(input, x_opt)

% close all
global DIRECTORY

if ~input.pos_select
    POSITION = 'Front Airbag';
else
    POSITION = 'Rear Airbag';
end

% %% 上半分で対称になるように配置
% des_num = length(x_opt); %設計変数の数
% 
% vent_H = x_opt(des_num/2+1:end);
% vent_H = [vent_H,flip(input.bag_Hs-vent_H)];
% 
% vent_D = x_opt(1:des_num/2);
% vent_D = [vent_D,flip(vent_D)];
% 
% x_opt = [vent_D,vent_H];

%% sort
x_opt = Para_sort(x_opt);  % 距離が近い順に設計変数をソート
para_calc_num = length(x_opt);
vent_D_sorted = x_opt(1:para_calc_num/2);
vent_H_sorted = x_opt(para_calc_num/2+1:end);


% 値整理
bag_D = input.bag_D * 1000;         % エアバッグ直径
bag_R = bag_D * 0.5;
bag_Hs = input.bag_Hs * 1000;       % エアバッグ側面高さ
bag_Hs_half = bag_Hs * 0.5;
bag_H = input.bag_H * 1000;         % エアバッグ高さ
bag_H_half = bag_H * 0.5;
bag_H_round = (bag_H-bag_Hs)/2;    % エアバッグ円弧部高さ(片側)

vent_D = [ x_opt(1);  x_opt(2);  x_opt(3);  x_opt(4);  x_opt(5);  x_opt(6)] * 1000;
vent_H = [ x_opt(7);  x_opt(8);  x_opt(9); x_opt(10); x_opt(11); x_opt(12)] * 1000;
vent_H = vent_H * 0.5;
[vent_H_sorted,I] = sort(vent_H,'descend');
vent_D_sorted = vent_D(I);

% ベントホールの有無を判定
vent_flag = [0 0 0 0 0 0];
for i = 1:1:6
    if vent_D_sorted(i) ~= 0
        vent_flag(i) = i;
    end
end

B = vent_flag == 0;
vent_H_sorted(B) = [];
vent_D_sorted(B) = [];
vent_r = vent_D_sorted * 0.5;
n = numel(vent_r);

% if n <= 3
%     toc
%     error('The number of vent holes is 3 or less.');
% end

%% グラフ設定
fName = 'Arial';
lWidth  = 1.0;  % 枠線
lWidth2 = 1.5;  % 枠線
fSize  = 20;
fSize2 = 16;

color = ["red" "#0072BD" "#77AC30" "#7E2F8E" "#EDB120" "#4DBEEE"];

nRow = 1;
nCol = 1;

r     = groot; % グラフィックスルートオブジェクトハンドルを格納
style = r.MonitorPositions(1,:);    % メインディスプレイのサイズ情報[x1 y1 width1 height1], (x1,y1)はモニタの原点座標, 左下隅原点

f(1) = figure('Name','Vent Hole') ; % clf(f(1));
f(1).OuterPosition = [style(3)/4 0.1*style(4) 0.4*style(3) 0.7*style(4)];    % [left bottom width height]

t = tiledlayout(nRow,nCol); % 複数のグラフをタイル表示, subplotと違い位置指定不要, グラフ間隔の変更が可能
nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
axis equal;
set(gca,'xlim',[-bag_R-170 bag_R+170],'ylim',[-bag_H_half-100 bag_H_half+100], ...
    'FontSize',fSize2,'Fontname',fName)
xlabel('x [mm]','FontSize',fSize2,'Fontname',fName);
ylabel('y [mm]','FontSize',fSize2,'Fontname',fName);
[m,s] = title(POSITION,'Unit: [mm]');
m.FontSize = 20;
s.FontSize = 14;

%% エアバッグ外形
% 側面
plot([-bag_R  bag_R],[-bag_Hs_half -bag_Hs_half],'linewidth',lWidth2,'color','k');
plot([-bag_R  bag_R],[ bag_Hs_half  bag_Hs_half],'linewidth',lWidth2,'color','k');
plot([-bag_R -bag_R],[-bag_Hs_half  bag_Hs_half],'linewidth',lWidth2,'color','k');
plot([ bag_R  bag_R],[-bag_Hs_half  bag_Hs_half],'linewidth',lWidth2,'color','k');

% 上下楕円部
a = 0:10:180;
x = bag_R * cosd(a);
y = bag_H_round * sind(a) + bag_Hs_half;
plot(x,y,'linewidth',lWidth2,'color','k')

a = 180:10:360;
y = bag_H_round * sind(a) - bag_Hs_half;
plot(x,y,'linewidth',lWidth2,'color','k')

% 寸法
% エアバッグ直径
plot([ bag_R  bag_R],[bag_Hs_half+20  bag_Hs_half+200],'linewidth',lWidth,'color','k');
plot([-bag_R -bag_R],[bag_Hs_half+20  bag_Hs_half+200],'linewidth',lWidth,'color','k');
q = quiver(0,bag_H_half+30,bag_R,0);
q.LineWidth = lWidth2;  q.Color = 'k';  q.AutoScale = "off";
q = quiver(0,bag_H_half+30,-bag_R,0);
q.LineWidth = lWidth2;  q.Color = 'k';  q.AutoScale = "off";
txt = text('String',num2str(bag_D),...
    'FontSize',12, 'FontName',fName, 'Position',[-25 bag_H_half+60 0],'Color','k');

% エアバッグ側面高さ
plot([ bag_R+10  bag_R+80],[ bag_Hs_half   bag_Hs_half],'linewidth',lWidth,'color','k');
plot([ bag_R+10  bag_R+80],[-bag_Hs_half  -bag_Hs_half],'linewidth',lWidth,'color','k');
q = quiver(bag_R+50,0,0,bag_Hs_half);
q.LineWidth = lWidth2;  q.Color = 'k';  q.AutoScale = "off";
q = quiver(bag_R+50,0,0,-bag_Hs_half);
q.LineWidth = lWidth2;  q.Color = 'k';  q.AutoScale = "off";
txt = text('String',num2str(bag_Hs),...
    'FontSize',12, 'FontName',fName,'Position',[bag_R+30 -25 0],'Color','k');
set(txt,'rotation',90)

% エアバッグ高さ
plot([ 150  bag_R+150],[ bag_H_half   bag_H_half],'linewidth',lWidth,'color','k');
plot([ 150  bag_R+150],[-bag_H_half  -bag_H_half],'linewidth',lWidth,'color','k');
q = quiver(bag_R+120,0,0,bag_H_half);
q.LineWidth = lWidth2;  q.Color = 'k';  q.AutoScale = "off";
q = quiver(bag_R+120,0,0,-bag_H_half);
q.LineWidth = lWidth2;  q.Color = 'k';  q.AutoScale = "off";
txt = text('String',num2str(bag_H),...
    'FontSize',12, 'FontName',fName,'Position',[bag_R+100 -25 0],'Color','k');
set(txt,'rotation',90)

%% ベントホール
% 2列間距離算出
vent_x = zeros(n-1,1);
for i = 1:1:n-1
    l_req = vent_r(i) + vent_r(i+1) + 25;
    h_req = vent_H_sorted(i) - vent_H_sorted(i+1);
    if l_req < h_req
        l_req = h_req;
    end
    vent_x(i) = sqrt(l_req^2 - h_req^2);
    vent_x(i) = ceil(vent_x(i));
end

x_req = max(vent_x);

% プロット
a = 0:10:360;
for i = 1:1:n
    x = vent_r(i) * cosd(a);
    y1 = vent_r(i) * sind(a) + vent_H_sorted(i);
    y2 = vent_r(i) * sind(a) - vent_H_sorted(i);
    plot(x,y1,'linewidth',lWidth2,'color',color(i))
    plot(x,y2,'linewidth',lWidth2,'color',color(i))

    txt = text('String',strcat('D=',num2str(vent_D_sorted(i))),...
        'FontSize',12, 'FontName',fName,'Position',[30 vent_H_sorted(i) 0],'Color',color(i));

    if rem(i, 2) == 0
    plot([-bag_R+30*i x_req],[ vent_H_sorted(i)  vent_H_sorted(i)],'linewidth',lWidth,'color',color(i));
    plot([-bag_R+30*i     0],[-vent_H_sorted(i) -vent_H_sorted(i)],'linewidth',lWidth,'color',color(i));
elseif rem(i, 2) == 1
    plot([-bag_R+30*i     0],[ vent_H_sorted(i)  vent_H_sorted(i)],'linewidth',lWidth,'color',color(i));
    plot([-bag_R+30*i x_req],[-vent_H_sorted(i) -vent_H_sorted(i)],'linewidth',lWidth,'color',color(i));
    end

q = quiver(-bag_R+50*i,0,0,vent_H_sorted(i));
q.LineWidth = lWidth2;  q.Color = color(i);  q.AutoScale = "off";   q.MaxHeadSize = 1;
q = quiver(-bag_R+50*i,0,0,-vent_H_sorted(i));
q.LineWidth = lWidth2;  q.Color = color(i);  q.AutoScale = "off";   q.MaxHeadSize = 1;
txt = text('String',num2str(vent_H_sorted(i)*2),...
    'FontSize',12, 'FontName',fName,'Position',[-25-bag_R+50*i -25 0],'Color',color(i));
set(txt,'rotation',90)

end

ymax = 300;
ymin = -300;
plot([0 0],[ymax ymin],'-.','linewidth',0.5,'color','k');


