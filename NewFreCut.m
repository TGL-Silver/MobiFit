

clc;
clear;
close all;

%% 参数设置
load E:\代码\QIAlgorithm\数据\步行\Data\20191211\corrider11_9428.txt;
Filename='tian7';
FilePath1='D:\6月24日代码数据汇总\代码数据汇总\数据\步行\特征\重新整理\20191104\Original';
FilePath2='D:\6月24日代码数据汇总\代码数据汇总\数据\步行\特征\重新整理\20191104\28';
flag=7;
% % load squat1;
%
x=1;
s=corrider11_9428(x:end);
s= hampel(s,256);

fs = 100;
T = 1/fs;
L = length(s);
t = (0:L-1)*T;
%% 原始图像
plot(t,s);title('','FontSize',18,'fontname','Times New Roman');
axis tight;set(gca,'FontSize',18);set(gca, 'LineWidth',1,'fontname','Times New Roman');
xlabel('Time(seconds)','FontSize',18,'Fontname', 'Times New Roman');ylabel('Amplitude','FontSize',18,'Fontname', 'Times New Roman');
set(gca,'FontSize',18)
set(gca, 'LineWidth',1.25)
%% 所有动作朝上切割 ，画图使用原始数据 stemp
 stemp=s;
 stempMean=mean(s(1:300));
 s=stempMean-s;
window = 512;                % number of section
step = 8;
overlap = window - step;         % number of overlap
nfft = window;                % number of fft

%% 预处理 归一化 STFT ,函数
sNorm = mapminmax(s',0,1);
stempNrom=mapminmax(stemp',0,1);
[S,W,T,tax]=NewFreCutSTFT(sNorm,window,overlap,nfft,fs);

%% 求一阶近似导数

[m, n] = size(S);

sPart1 = S(:,1:n - 4);
sPart2 = S(:,5:n);

searchMutations = sPart2 - sPart1;
radio = searchMutations;
radioNorm = mapminmax(radio,0.25,0.75);

%% 选取0.8Hz频率线
% radioNormfinal=(radioNorm(2,:)+radioNorm(3,:))/2;
e =smoothdata(radioNorm(2,:),'gaussian',80);
% e =radioNorm(3,:);
eSort =sort(e);

figure;
plot(tax,S(2,:),'LineWidth',1.25);axis tight;
hold on;plot(tax,S(3,:),'LineWidth',1.25);axis tight;
hold on;plot(tax,S(4,:),'LineWidth',1.25);axis tight;
legend('0.4Hz','0.8Hz','1.2Hz');
set(gca, 'LineWidth',1.25)
set(gca,'FontSize',18)
xlabel('Time(s)','FontSize',18,'Fontname', 'Times New Roman');
ylabel('Frequency','FontSize',18,'Fontname', 'Times New Roman');

%% 绘制图像 
figure;plot(t,sNorm);title('Walk Cutting Result','FontSize',18,'fontname','Times New Roman');
set(gca, 'LineWidth',1.25);

%% 绘制0.2Hz频率线
freLine = S(3,:);
freLineNorm = mapminmax(freLine,0,1);
hold on;plot(tax,freLineNorm,':m','LineWidth',2);

%% 绘制近似一阶导数
beta=0.88;
alpha=0.12;
hold on;plot(tax(1:n - 4),smoothdata(radioNorm(2,:),'gaussian',80),'black','LineWidth',2);axis tight;
% figure;plot(t,sNorm,'Color',[0.7734375 0.3046875 0.16796875]);
num = ceil(beta* length(eSort));
num2 = ceil(alpha* length(eSort));
mostLine = eSort(num);
lowLine = eSort(num2);
% hold on;line([0 t(end)],[mostLine mostLine],'Color',[1 0.5 0.5],'LineStyle',':','LineWidth',2);
% hold on;line([0 t(end)],[lowLine lowLine],'Color',[1 0.5 0.5],'LineStyle',':','LineWidth',2);
eIndexTime = [];
eIndexNum = [];
eNum = [];
eIndexFre=[];
for i =2:(length(e)-1)
    % % %     %过滤既是峰值又大于阈值的点
    if((e(i) > e(i-1)) && (e(i) > e(i+1)) && e(i) > mostLine)||((e(i) < e(i-1)) && (e(i) < e(i+1)) && e(i) < lowLine)
        eIndexTime = [eIndexTime tax(i)];
        eIndexNum = [eIndexNum fix(tax(i)*fs)];
        eNum = [eNum e(i)];
        temp=fix(tax(i)*fs);
        tempnow=(temp-512-255)/8;
        eIndexFre=[eIndexFre i];
    end
end


hold on;plot(eIndexTime,eNum,'*r');

hold on;hh = axis;
for i = eIndexTime
    plot([i,i], [hh(3),hh(4)],'Color',	'black','LineWidth',2);
end

set(gca,'FontSize',18)
title('  Cutting Result','FontSize',18,'fontname','Times New Roman');
xlabel('Time(s)','FontSize',18,'Fontname', 'Times New Roman');
ylabel('Amplitude','FontSize',18,'Fontname', 'Times New Roman');
legend('Normalized signal','Main frequency line ','Derivative','Peak and Valley');
set(legend, 'fontsize',10);
set(gca, 'LineWidth',1.25);

%% 若切割准确使用奇偶分类
for i = 1:length(eIndexNum) - 1
    cutSignal(i,:) = eIndexNum(i:i + 1);
end
cutSignal = [[1 eIndexNum(1)];cutSignal;[eIndexNum(end) length(sNorm)]];
cutNum=cutSignal(:,2)-cutSignal(:,1);
cutSignalFre_2n_Num = cutSignal(2:2:end,:);
cutSignalFre_2n_1_Num = cutSignal(1:2:end,:);
%%  计算各种能量值
for i = 1:length(eIndexFre) - 1
    cutSignalFre(i,:) = eIndexFre(i:i + 1);
end
cutSignalFre = [[1 eIndexFre(1)];cutSignalFre;[eIndexFre(end) length(S)]];
cut=cutSignalFre(:,2)-cutSignalFre(:,1);
cutSignalFre_2n = cutSignalFre(2:2:end,:);
cutSignalFre_2n_1 = cutSignalFre(1:2:end,:);
% cut_2n = cut(2:2:end,:);
% cut_2n_min=min(cut_2n);
% cut_2n_max=max(cut_2n);
% cut_2n_mean=mean(cut_2n);
% cut_2n_1 = cut(1:2:end,:);
% cut_2n_1_min=min(cut_2n_1);
% cut_2n_1_max=max(cut_2n_1);
% cut_2n_1_mean=mean(cut_2n_1);
% CutSignalFreSum=[];
% CutSignalFreSumRadio=[];
% Frenow=S(2,:);
% for i=1:size(cutSignalFre,1)
%     CutSignalFreSum(i)=(sum(Frenow(cutSignalFre(i,1):cutSignalFre(i,2))).^2)/(cutSignalFre(i,2)-cutSignalFre(i,1));
% end
% CutSignalFreSum=CutSignalFreSum';
% CutSignalFreSum_2n = CutSignalFreSum(2:2:end,:);
% CutSignalFreSum_2n_1 = CutSignalFreSum(1:2:end,:);
% CutSignalFreSum_2n_1=CutSignalFreSum_2n_1(2:end,1);
% CutSignalFreSumRadio=CutSignalFreSum_2n./CutSignalFreSum_2n_1;
% CutSignalFreSumRadio_min=min(CutSignalFreSumRadio);
% CutSignalFreSumRadio_max=max(CutSignalFreSumRadio);
% CutSignalFreSumRadio_mean=mean(CutSignalFreSumRadio);
%% 原始数据切割图像
titlenow='Counting Result';
xlablenow='Sampling points';
ylabelnow='Amplitude';
QiPlot(stempNrom,titlenow,xlablenow,ylabelnow);
hh = axis;hold on;
for i = eIndexNum
    plot([i,i], [hh(3),hh(4)],'Color','black','LineWidth',2);
end
%% 动态切割 以多少次STFT为标志位 频率位置
ActionALL=WalkenergyCalucate(20,80,cutSignalFre_2n,sNorm);
% 计算动作从持续时间
all1=mean(ActionALL(2,:)-ActionALL(1,:));
ActionALL=ActionALL';
%% 特征提取  以采样点位标志位  时域位置
% FeatureExtraction(s,sNorm,ActionALL,flag,fs,Filename,FilePath1);
% addpath 'D:\6月24日代码数据汇总\代码数据汇总\QIAlgorithm'
% FeatureExtraction(s,sNorm,cutSignalFre_2n_Num,flag,fs,Filename,FilePath2);


