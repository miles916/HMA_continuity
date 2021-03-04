%% go to folder, load files, paths
clear
clc
close all


hdir = 'C:\Users\miles\Documents\GitHub\Huss2010_flowparam\Committed_loss_calcs'; 
cd(hdir)

addpath(genpath('C:\Users\miles\Documents\GitHub\Huss2010_flowparam'))
addpath('\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\debris_fluxes')
% Rongbuk
% odir='Rongbuk_test_fixedSMB';
% cd C:\Users\miles\Documents\GitHub\Huss2010_flowparam\Committed_loss_calcs\test_Rongbuk_07-Aug-2020\
% DH=double(imread('Rongbuk_dH.tif'));
% DEM=double(imread('Rongbuk_DEM.tif'));
% THX=double(imread('Rongbuk_THX.tif'));
% SMB=double(imread('Rongbuk_zSMB.tif'));
% info=geotiffinfo('Rongbuk_DEM.tif');
% dx=info.PixelScale(1);
% SMBtype=2; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2

% % %Aletsch
% odir='Aletsch_hybridSMB';
% cd C:\Users\miles\Documents\GitHub\Huss2010_flowparam\Committed_loss_calcs\test_Aletsch_07-Aug-2020
% DH=double(imread('Aletsch_dH.tif'));
% DEM=double(imread('Aletsch_DEM.tif'));
% THX=double(imread('Aletsch_THX.tif'));
% SMB=double(imread('Aletsch_zSMB.tif'));
% info=geotiffinfo('Aletsch_DEM.tif');
% dx=info.PixelScale(1);
% SMBtype=4; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2

%24K
% odir='24k_cLoss_mixedSMB';
% cd C:\Users\miles\Dropbox\collaborations\24K_emergence\15.11758
% DH=double(imread('15.11758_dH.tif'));
% DEM=double(imread('15.11758_AW3D.tif'));
% THX=double(imread('15.11758_THX.tif'));
% SMB=double(imread('15.11758_zSMB.tif'));
% info=geotiffinfo('15.11758_AW3D.tif');
% dx=info.PixelScale(1);
% SMBtype=4; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2

% % Langtang
% odir='Langtang_cLoss_fixedSMB';
% cd \\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\ValidationSites_12-May-2020
% DH=double(imread('15.04121_dH.tif'));
% DEM=double(imread('15.04121_AW3D.tif'));
% THX=double(imread('15.04121_THX.tif'));
% SMB=double(imread('15.04121_zSMB.tif'));
% info=geotiffinfo('15.04121_AW3D.tif');
% dx=info.PixelScale(1);
% SMBtype=2; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2
% oname = '15.04121';

% %Satopanth
% odir='Satopanth_mixedSMB';
% cd \\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\ValidationSites_12-May-2020
% DH=double(imread('15.07122_dH.tif'));
% DEM=double(imread('15.07122_AW3D.tif'));
% THX=double(imread('15.07122_THX.tif'));
% SMB=double(imread('15.07122_zSMB.tif'));
% info=geotiffinfo('15.07122_AW3D.tif');
% dx=info.PixelScale(1);
% SMBtype=4; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2

% %7084
% odir='7084_mixedSMB';
% cd \\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\RGI131415_ALL2km_13-May-2020
% DH=double(imread('15.07084_dH.tif'));
% DEM=double(imread('15.07084_AW3D.tif'));
% THX=double(imread('15.07084_THX.tif'));
% SMB=double(imread('15.07084_zSMB.tif'));
% info=geotiffinfo('15.07084_AW3D.tif');
% dx=info.PixelScale(1);
% SMBtype=4; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2

% %17894
% odir='13.17894_mixedSMB';
% cd \\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\RGI131415_ALL2km_13-May-2020
% DH=double(imread('13.17894_dH.tif'));
% DEM=double(imread('13.17894_AW3D.tif'));
% THX=double(imread('13.17894_THX.tif'));
% SMB=double(imread('13.17894_zSMB.tif'));
% info=geotiffinfo('13.17894_AW3D.tif');
% dx=info.PixelScale(1);
% SMBtype=4; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2

% Kangjiaruo Glacier - opposite Langtang
odir='15.10299_fixedSMB';
cd C:\Users\miles\Documents\GitHub\Huss2010_flowparam\15.10299
DH=double(imread('15.10299_dH.tif'));
DEM=double(imread('15.10299_AW3D.tif'));
THX=double(imread('15.10299_THX.tif'));
SMB=double(imread('15.10299_zSMB.tif'));
info=geotiffinfo('15.10299_AW3D.tif');
dx=info.PixelScale(1);
SMBtype=2; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the 2
oname = '15.10299';

% go to output folder
cd(hdir)
mkdir(odir)
cd(odir)

gifout=1;
P2 = 0.5;% P2 is partition of mass gain - best guess is ~0.483 (48%) but equifinality with P4
P4 = 11.8; % P4 is the terminus longitudinal gradient allowing advance - best guess is 11.8deg with P2=50%
H2010=2; %configuration of flow-parameterization: 1 is with area-based curves (precisely in H2010) and 2 is with each glacier's DH curve

%measure glacier width of terminus
mask=THX>0;
widths=measureTerminusWidths(DEM,mask,dx);close
P1 = mean(widths)./dx;% P1 is N terminus pixels to average over - determine as glacier width
P3 = P1;% P3 is the number of pixels to advance into - best guess is equalt to P1
%%
% [PctComLoss,etime,a]=comLossFxn(DH,DEM,THX,SMB,dx,SMBtype,oname,gifout)
[PctComLoss,etime,a]=comLossFxn_complex_lgrad(DH,DEM,THX,SMB,dx,SMBtype,oname,gifout,P1,P2,P3,P4,H2010)

%% plots
% final plots
figure
yyaxis left
plot(a.Yr,a.pVol.*100,'-');hold on 
plot(a.Yr,a.pArea.*100,':');
% plot(a.Yr,a.dVolcl.*10000,':');
ylabel('Percent of present-day')
ylim([0,max([100,100.*a.pVol(end),100.*a.pArea(end)])])
yyaxis right
plot(a.Yr,a.gMB,'-');hold on 
legend('Volume','Area','MB','Location','east')
ylabel('Mass balance (m w.e. a^{-1})')
xlim([min(a.Yr),max(a.Yr)])
saveas(gcf,[oname '_glacier_evolution.png'])


%% plot ablbal

% final plots
figure
yyaxis left
area(a.Yr,-1.*[a.Balabl;a.Imbalabl]','FaceColor','flat','FaceAlpha',0.5,'LineStyle','none');
line(a.Yr,-1.*a.Totabl,'Linewidth',1.5);
ylabel('Ablation (m^3 w.e. a^{-1})')
% ylim([0,100])
yyaxis right
plot(a.Yr,a.AAR,'-','Linewidth',1.5);hold on 
legend('Balance','Imbalance','Total','AAR','Location','north')
ylabel('AAR (-)')
xlim([min(a.Yr),max(a.Yr)])
saveas(gcf,[oname 'glacier_evolution2.png'])
