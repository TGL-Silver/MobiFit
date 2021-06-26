

clc;
clear;
close all;

%% ��������
% load D:\6��24�մ������ݻ���\�������ݻ���\����ʶ��\20191111\Data\shuju3_teng.txt
% load shuju3_teng.mat;
% x=1;

% s=shuju3_teng(3000:11400);
% sfinal=s(8000:end);
% s=[s;sfinal];
% load mix(1).mat;
% load multiple.mat;
% s=shuju3_teng;
% s= hampel(s,256);
% load D:\6��24�մ������ݻ���\�������ݻ���\����ʶ��\����\data\fwc\qjb\2019-1-10\fwc1.txt
% load fwc.mat;
load C:\Users\92854\Desktop\QIAlgorithm\����\����\Data\20191106\officeteng.txt
x=1;
s=officeteng(x:end);
s= hampel(s,256);
fs = 200;
T = 1/fs;
L = length(s);
t = (0:L-1)*T;
%% ����ԭʼͼ��
stemp=s;
%% ����Ҫ�����и�
% stempMean=mean(s(1:300));
% s=stempMean-s;

%% ԭʼͼ��

window = 512;                % number of section
step = 8;
nov = window - step;         % number of overlap
nff = window;                % number of fft

%% Ԥ���� ��һ��
sNorm = mapminmax(s',0,1);
TStart=[];
TEnd=[];
Ftime=[];
CellularS=zeros(1,500);
CellularE=zeros(1,500);
count=2;
tempcount=2;
CellularE(1)=1;
CellularS(1)=1;

d=4;
Derivative=[];
alpha1=0.8;
alpha2=0.75;
beta=1.5;
% Restmin=57;
% Actionmin=30;
Restmin=156;
Actionmin=41;
tempRestmin=Restmin;
tempActionmin=Actionmin;
tempalpha=alpha;

% varargout = untitled1();
global gloablcount;
global sGui;
global CellularALL;
CellularALL=[];
gloablcount=count;
sGui=[];
varargout = MobiFit();
for qiposition=512:512:L
%         pause(1);
    %
    
    %     plot(stemp(1:qiposition));
    %     title(' Original signal  ','FontSize',18,'fontname','Times New Roman');
    %     axis tight;set(gca,'FontSize',18);set(gca, 'LineWidth',1,'fontname','Times New Roman');
    %     xlabel('Time(seconds)','FontSize',18,'Fontname', 'Times New Roman');
    %     ylabel('Amplitude','FontSize',18,'Fontname', 'Times New Roman');
    %     set(gca,'FontSize',18);
    %     set(gca, 'LineWidth',1.25);
    
    sGui=stemp(1:qiposition);
    varargout = MobiFit();
    
    for a=qiposition-511:qiposition
        if a-511>0
            if(mod((a-window),step)==0)
                CellularPs=(a-window)/step+1; %%����ڼ���STFT
                F(:,CellularPs)=CellularSTFT(200,sNorm(a-511:a)');
                
                midwindow=.5*window;
                if(CellularPs==1)
                    Ftime(CellularPs)=midwindow;
                else
                    Ftime(CellularPs)= Ftime(CellularPs-1)+8;
                end
                
                if(CellularPs>=5) %���4λ��λ���
                    CellularPd=CellularPs-d;
                    Derivative(:,CellularPd)=F(:,CellularPs)-F(:,CellularPd);%��������
                    if(size(Derivative,2)>=3)
                        
                        %               [pks,locs] =findpeaks(Derivative(2,:)); %%���ҳ���ֵ�͹�ֵ ����
                        %               [pksvalley,locsvalley] =findpeaks(-Derivative(2,:));%%���ҳ���ֵ�͹�ֵ ����
                        
                        
                        if (CellularPd>=174) && (CellularPd<=176)
                            CellularPd=CellularPd;
                        end
                        
                        Derivativetemp=Derivative(2,:); %�ڶ���������
                        %�����ֵ
                        %                           if (Derivativetemp(CellularPd-1)>0.02)
                        if (Derivativetemp(CellularPd-1)>Derivativetemp(CellularPd-2)&&Derivativetemp(CellularPd-1)>Derivativetemp(CellularPd))
                            %                   &&Derivativetemp(CellularPd-1)>0.02)
                            
                            peaknow=Derivativetemp(CellularPd-1);
                            % ��ǰ�� ��ȥ��һ�������Ľ����� Ҫ����һ����Ϣ����ķ���
                            if((CellularPd-1)-CellularE(count-1)>alpha1*Restmin)
                                
                                if CellularS(count)==0
                                    CellularS(count)=CellularPd-1;
                                else
                                    %���·�ֵ
                                    if(Derivativetemp(CellularPd-1)>Derivativetemp( CellularS(count)))
                                        CellularS(count)=CellularPd-1;
                                    end
                                end
                                TStart(count)=CellularS(count)*step+window/2-1;
                            end
                        else
                            if CellularS(count)~=0
                                if CellularPd==87
                                    CellularPd=87;
                                end
                                
                                %�����ֵ
                                %                                 if (Derivativetemp(CellularPd-1)<-0.02)
                                if(Derivativetemp(CellularPd-1)<=Derivativetemp(CellularPd-2)&&Derivativetemp(CellularPd-1)<Derivativetemp(CellularPd))
                                    %                       &&Derivativetemp(CellularPd-1)<-0.02
                                    valleynow=Derivativetemp(CellularPd-1);
                                    % ��ֵ��λ�ü�ȥ��һ�������Ŀ�ʼʱ��Ҫ����һ�������ĳ���ʱ��ķ�������ǰ�����������ͽ���ǰһ����Ϣ�ı�ֵҪ����alpha
                                    
                                    % �����һ�����������ֵ
                                    if count==2
                                        sumNow=((sum(F(2,CellularS(count):(CellularPd-1))).^2)/((CellularPd-1)-CellularS(count)))>alpha2*((sum(F(2,CellularS(count)-fix(Restmin*beta1):CellularS(count))).^2)/(CellularS(count)-(CellularS(count)-fix(Restmin*beta1))));
                                    else
                                        sumNow=((sum(F(2,CellularS(count):(CellularPd-1))).^2)/((CellularPd-1)-CellularS(count)))>alpha2*((sum(F(2,CellularE(count-1):CellularS(count))).^2)/(CellularS(count)-CellularE(count-1)))     ;
                                    end
                                    
                                    if ((CellularPd-1)-CellularS(count)>beta2*Actionmin) && sumNow
                                        %                                     if ((CellularPd-1)-CellularS(count)>beta2*Actionmin)  && ((sum(F(2,CellularS(count):(CellularPd-1))).^2)/((CellularPd-1)-CellularS(count)))>alpha*((sum(F(2,CellularE(count-1):CellularS(count))).^2)/(CellularS(count)-CellularE(count-1)))
                                        if CellularE(count)==0
                                            CellularE(count)=CellularPd-1;
                                            count=count+1;
                                            
                                        else
                                            %�ҵ��˸��͵Ĺ�ֵ����ǰ�ȵ�λ�ü�ȥ��һ����ֵ��λ��ҪС��һ����Ϣ��϶
                                            if Derivativetemp(CellularPd-1)<Derivativetemp( CellularE(count-1)) && (CellularPd-1)-CellularE(count-1)<beta3*Restmin
                                                CellularE(count-1)=CellularPd-1;
                                                
                                            end
                                            
                                        end
                                        TEnd(count)=CellularE(count)*step+window/2-1;
                                        
                                        
                                        CellularALL=[CellularS;CellularE];
%                                         if(count>tempcount)
%                                             %������С����ʱ��
%                                             tempActionmin=CellularALL(2,tempcount)-CellularALL(1,tempcount);
%                                             if tempActionmin<Actionmin
%                                                 Actionmin=tempActionmin;
%                                             end
%                                             %������С��Ϣʱ��
%                                             if tempcount>2
%                                                 tempRestmin=CellularALL(1,tempcount)-CellularALL(2,tempcount-1);
%                                                 if tempRestmin<Restmin
%                                                     Restmin=tempRestmin;
%                                                 end
%                                                 
%                                                 tempalpha=(sum(F(2,CellularALL(1,tempcount):CellularALL(2,tempcount)).^2))/((CellularALL(2,tempcount)-CellularALL(1,tempcount))) / (sum(F(2,CellularALL(2,tempcount-1):CellularALL(1,tempcount)).^2)) / ((CellularALL(1,tempcount)-CellularALL(2,tempcount-1)));
%                                                 
%                                                 if tempalpha<alpha
%                                                     alpha=tempalpha;
%                                                 end
%                                                 
%                                             end
%                                             %�������� alpha
%                                             tempcount=tempcount+1;
%                                             
%                                         end
%                                         
                                        
                                        gloablcount=count;
                                        varargout = MobiFit();
                                        %                                         hold on;
                                        %                                         value1 = 1;
                                        %                                         value2 = 1;
                                        %                                         i = 2;
                                        %                                         while value1~=0 &&value2~=0
                                        %                                             plot([CellularALL(1,i)*8+255,CellularALL(1,i)*8+255],[min(s),max(s)],'-g');
                                        %                                             plot([CellularALL(2,i)*8+255,CellularALL(2,i)*8+255],[min(s),max(s)],'-r');
                                        %                                             value1 = CellularALL(1,i+1);
                                        %                                             value2 = CellularALL(2,i+1);
                                        %                                             i=i+1;
                                        %                                         end
                                        %
                                        
                                        
                                        
                                        
                                    end
                                    
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
    end
    
end
% %%   ���Է�ֵ ��Ҫ +4 -4 ����׼ȷ
% %  plot(F(2,:));
% %  temp=F(2,:);
% % for i=1:length(locs);
% %     peaksnow(i)=temp(locs(i)-4);
% % end
% % hold on; plot(locs,peaksnow,'*r');
%

%% ����ż����
% cut=cutSignalFre(:,2)-cutSignalFre(:,1);
% cutSignalFre_2n = cutSignalFre(2:2:end,:);
% cutSignalFre_2n_1 = cutSignalFre(1:2:end,:);
% %���㶯��ʱ��
% cut_2n = cut(2:2:end,:);
% cut_2n_min=min(cut_2n);
% cut_2n_max=max(cut_2n);
% cut_2n_mean=mean(cut_2n);
% %������Ϣʱ��
% cut_2n_1 = cut(1:2:end,:);
% cut_2n_1_min=min(cut_2n_1);
% cut_2n_1_max=max(cut_2n_1);
% cut_2n_1_mean=mean(cut_2n_1);
% %��������
% CutSignalFreSum=[];
% CutSignalFreSumRadio=[];
% Frenow=F(2,:);
% for i=1:size(cutSignalFre,1)
%     CutSignalFreSum(i)=(sum(Frenow(cutSignalFre(i,1):cutSignalFre(i,2))).^2)/(cutSignalFre(i,2)-cutSignalFre(i,1));
% end
% CutSignalFreSum=CutSignalFreSum';
% CutSignalFreSum_2n = CutSignalFreSum(2:2:end,:);
% CutSignalFreSum_2n_1 = CutSignalFreSum(1:2:end,:);
% CutSignalFreSum_2n_1=CutSignalFreSum_2n_1(1:end,:);
% CutSignalFreSumRadio=CutSignalFreSum_2n./CutSignalFreSum_2n_1;
% CutSignalFreSumRadio_min=min(CutSignalFreSumRadio);
% CutSignalFreSumRadio_max=max(CutSignalFreSumRadio);
% CutSignalFreSumRadio_mean=mean(CutSignalFreSumRadio);

