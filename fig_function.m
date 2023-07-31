function  fig_function(ts,result,input,varargin)%(ts,result, input)

fName = 'Times New Roman'; % 'Arial';
lWidth = 1.0;  % 枠線
lWidth2 = 2.0;  % 枠線
fSize  = 20;
fSize2 = 16;

nRow = 4;
nCol = 3;

r     = groot; % グラフィックスルートオブジェクトハンドルを格納
style = r.MonitorPositions(1,:);    % メインディスプレイのサイズ情報[x1 y1 width1 height1], (x1,y1)はモニタの原点座標, 左下隅原点

f(1) = figure('Name','Result') ; clf(f(1));
f(1).OuterPosition = [style(3)/4 0.1*style(4) 0.7*style(3) 0.8*style(4)];    % [left bottom width height]

t = tiledlayout(nRow,nCol); % 複数のグラフをタイル表示, subplotと違い位置指定不要, グラフ間隔の変更が可能

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,1),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('x [m]','FontSize',fSize2,'Fontname',fName);
%         legend('3 vent holes','4 vent holes','pettern 3','pettern 4','pettern 5','pettern 6','pettern 7','pettern 8','pettern 9')

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,2),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('v [m/s]','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,7),'linewidth',lWidth2);
yline(7,'-r','7G');
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('a [G]','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
plot(ts,result(:,3) - input.p/1000,'linewidth',lWidth2);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('P [kPa]','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,5)-273.15,'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('Gas Temperature[oC]','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,13),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('mdot [kg/s]','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,9),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('Mech Number','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,6),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('Volume[m3]','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,4),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('Gas Mass [kg]','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,12),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('number of vent hole','FontSize',fSize2,'Fontname',fName);

nexttile;
hold on; box on; grid on;
set(gca,'fontname',fName,'fontsize',fSize,'linewidth',lWidth,'XColor','k');
set(gcf,'Color',[1 1 1]);
set(gca,'xtick',0:input.dt_plot:input.te_plot,'xlim',[0 inf],'FontSize',fSize2,'Fontname',fName)
plot(ts,result(:,10),'linewidth',lWidth2);
xlabel('Time [s]','FontSize',fSize2,'Fontname',fName);
ylabel('vent area [m2]','FontSize',fSize2,'Fontname',fName);


t.TileSpacing = 'compact';
t.Padding = 'compact';
drawnow;
clear ts result
end