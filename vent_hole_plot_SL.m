function vent_hole_plot_SL(input,para)

if ~input.pos_select
    POSITION = 'Front Airbag';
else
    POSITION = 'Rear Airbag';
end

% 値整理
bag_D = input.bag_D * 1000;         % エアバッグ直径
bag_R = bag_D * 0.5;
bag_Hs = input.bag_Hs * 1000;       % エアバッグ側面高さ
bag_Hs_half = bag_Hs * 0.5;
bag_H = input.bag_H * 1000;         % エアバッグ高さ
bag_H_half = bag_H * 0.5;
bag_H_round = (bag_H-bag_Hs)/2;    % エアバッグ円弧部高さ(片側)

%% グラフ設定
fName = 'Arial';
lWidth  = 1.0;  % エアバッグ枠線
lWidth2 = 1.5;  % ベントホール枠線
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
vent_num = length(para);
vent_D = 1000*para(1:vent_num/2);
vent_H = 1000*para(vent_num/2+1:end);
D = vent_D(1);
vent_r = vent_D/2;

% 最小穴ふち距離
dif_min = 25;

% 2列間距離算出
H_dif = abs(vent_H(1) - vent_H(2));
x_req = sqrt((D+dif_min)^2-H_dif^2);
x_req = ceil(x_req); 

% プロット
theta = 0:10:360;
for i=1:1:length(vent_H)
    x = vent_r(i)*cosd(theta);
    if rem(i,2) == 0
        x = x + x_req;
    end
    y = vent_r(i)*sind(theta) + vent_H(i);

    plot(x,y,'linewidth',lWidth2,'Color','b')
end

% 中心線
ymax = 300;
ymin = -300;
plot([0 0],[ymax ymin],'-.','linewidth',0.5,'color','k');



