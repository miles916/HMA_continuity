%% go to folder, load files, paths
clear
clc
close all

% addpath(genpath('C:\Users\miles\Documents\GitHub\Huss2010_flowparam'))
addpath(genpath('\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Huss2010_flowparam'))
addpath(('\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Committed_loss_calcs'))

% %set path to analyze, load list of files
% hdir = '\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\RGI131415_ALL2km_13-May-2020'; 

%load input data
glacier_postproc_results = '\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\POST_2km2_15-May-2020\'
X=2;
cd(glacier_postproc_results)
res=dir('ELA*.mat');
load(res(end).name);
VALID2 = ~g.issurging & g.ELAacc>0.5 & g.ELAs>g.Zmin & (g.DHstd)<3 & (g.MBstd_C)<3 &(g.MBmeansig<5) & (g.DHslope'>-0.5 );
VALID=VALID2 &g.Area>X;

%set output directory
% odir='\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Committed_loss_calcs\ComLoss_SMB2_MC_25-Aug-2020';
% odir='\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Committed_loss_calcs\ComLoss_SMB2_tgrad_width_MC_30-Oct-2020';
odir='\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Committed_loss_calcs\ComLoss_SMB2_tgrad_width_MC_17-Nov-2020';
% odir='\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Committed_loss_calcs\ComLoss_SMB2_tgrad_width_MC_23-Nov-2020';
% odir = 'C:\Users\miles\Desktop\CommLossTemp';
% % go to output folder
cd(odir)

%%

aCom.PctComLoss=NaN*g.AAR;
aCom.etime=NaN*g.AAR;
ers=0.*g.AAR;

tic
for ig=1:length(g.AAR)
    if VALID(ig)==1
        disp([ig length(g.AAR) sum(VALID(1:ig))])
        close all
        oname=g.RGIid{ig}(end-7:end);

        %load current output
        try
        gCom=load(fullfile(odir,[oname '.mat']));

        %% index
%         aCom.a(ig)=a;
        aCom.PctComLoss(ig)=gCom.M.PctComloss;
        aCom.PctVol2100(ig)=gCom.M.pVol(100);
        aCom.PctComLoss_u(ig)=gCom.U.PctComloss;
        aCom.PctVol2100_u(ig)=gCom.U.pVol(100);
        
        aCom.PctComLossA(ig)=gCom.M.pArea(end);
        aCom.PctA2100(ig)=gCom.M.pArea(100);
        aCom.PctComLossA_u(ig)=gCom.U.pArea(end);
        aCom.PctA2100_u(ig)=gCom.U.pArea(100);
        
        aCom.etime(ig)=gCom.M.etime;
        aCom.A0(ig)=gCom.M.Area(1);
        aCom.V0(ig)=gCom.M.Vol(1);
        aCom.gMB(ig,1:1000)=gCom.M.gMB(:);
        aCom.gMBU(ig,1:1000)=gCom.U.gMB(:);
        aCom.Vol(ig,1:1000)=gCom.M.Vol(:);
        aCom.VolU(ig,1:1000)=gCom.U.Vol(:);
        aCom.Area(ig,1:1000)=gCom.M.Area(:);
        aCom.AreaU(ig,1:1000)=gCom.U.Area(:);
        aCom.Balabl(ig,1:1000)=gCom.M.Balabl(:);
        aCom.BalablU(ig,1:1000)=gCom.U.Balabl(:);
        aCom.Imbalabl(ig,1:1000)=gCom.M.Imbalabl(:);
        aCom.ImbalablU(ig,1:1000)=gCom.U.Imbalabl(:);
        aCom.Totabl(ig,1:1000)=gCom.M.Totabl(:);
        aCom.TotablU(ig,1:1000)=gCom.U.Totabl(:);
        aCom.AAR(ig,1:1000)=gCom.M.AAR(:);
        aCom.AARu(ig,1:1000)=gCom.U.AAR(:);
        catch
             ers(ig)=1;
        end
    end
end
toc
disp(nansum(aCom.PctComLoss(VALID)>0))
% disp(sum(VALID(1:ig)))
%%
save('ers.mat','ers')
save(fullfile(odir,['summary_' date '.mat']))
% cd(odir)

% %% etime plots
% aCom.etime(aCom.etime==0)=NaN;
% figure;
% subaxis(1,2,1)
% scatter(g.X(VALID),g.Y(VALID),[],aCom.etime(VALID),'filled');c=colorbar
% ylabel(c,'e-folding time (yrs)')
% subaxis(1,2,2)
% histogram(aCom.etime(VALID),0:25:500);
% xlabel('e-folding time (yrs)')
% ylabel('N glaciers')
% 
% saveas(gcf,fullfile(odir,'efolding_summary.png'))
% 
%% comloss plots
% aCom.PctComLoss(isnan(aCom.etime))=NaN;
figure;
subaxis(1,2,1)
scatter(g.X(VALID),g.Y(VALID),[],-aCom.PctComLoss(VALID),'filled');c=colorbar;caxis([-100,100])
ylabel(c,'Committed mass change (%)')
subaxis(1,2,2)
histogram(-aCom.PctComLoss(VALID),-100:10:200);xlim([-100,200])
xlabel('Committed mass change (%)')
ylabel('N glaciers')

saveas(gcf,fullfile(odir,'pctComLoss_summary.png'))
% 
% %% total committed loss.
% ival = find(VALID); %index of glaciers actually processed
% % ival = find(VALID(1:ig)); %index of glaciers actually processed
% VolN=zeros(1,length(ival));
% VolF=zeros(1,length(ival));
% Vol2100=zeros(1,length(ival));
% Vol2200=zeros(1,length(ival));
% for iv=1:length(ival)
%     try
%     VolN(iv)=aCom.V0(iv);
%     VolF(iv)=aCom.V0(iv).*aCom.PctComLoss(iv);
%     Vol2100(iv)=aCom.V0(iv).*aCom.PctVol2100(iv);
%     Vol2100u(iv)=aCom.V0(iv).*aCom.PctVol2100_u(iv);
%     Vol2200(iv)=aCom.Vol(iv,200);
%     Vol2200u(iv)=aCom.VolU(iv,200);
%     catch
%     end
% end
% lig=zeros(1,max(ival));
% for ig2=1:length(aCom.A0) %will be shorter if only subset processed
%     try
%         lig(ig2)=length([aCom.gMB(ig2,:)]);
%     catch
%     end
% end
% ival=ival(lig(ival(ival<length(lig)))>0); %index relative to ig
% [ival2]=(find(lig>0)); %index relative to valid glaciers
% 
% % VolN = (g.VOL(ival));
% % VolF = VolN.*(1-[aCom.PctComLoss(ival)]/100);
% Vol2100((VolN>1e13)|isinf(VolN)|isinf(VolF))=NaN; %something wrong with 1 glacier 10^30 m3?!
% Vol2100u((VolN>1e13)|isinf(VolN)|isinf(VolF))=NaN; %something wrong with 1 glacier 10^30 m3?!
% % Vol2200((VolN>1e13)|isinf(VolN)|isinf(VolF))=NaN; %something wrong with 1 glacier 10^30 m3?!
% VolF((VolN>1e13)|isinf(VolN)|isinf(VolF))=NaN; %something wrong with 1 glacier 10^30 m3?!
% VolN((VolN>1e13)|isinf(VolN)|isinf(VolF))=NaN; %something wrong with 1 glacier 10^30 m3?!
% 
% tVolN=nansum(VolN);
% tVolF=nansum(VolF);
% tVol2100=nansum(Vol2100);
% tVol2100u=sqrt(nansum(Vol2100u.^2));
% tVol2200=nansum(Vol2200);
% tVol2200u=sqrt(nansum(Vol2200u.^2));
% % tVol2200=nansum(Vol2200);
% 
% aPctComLoss=(tVolN-tVolF)/tVolN
% aPctVol2100=(tVol2100)/tVolN
% aPctVol2100u=(tVol2100u)/tVolN
% aPctVol2200=(tVol2200)/tVolN
% aPctVol2200u=(tVol2200u)/tVolN
% 
% %% loss by 2100 and 2200
% 
% % aCom.PctComLoss(isnan(aCom.etime))=NaN;
% figure;
% subplot(1,2,1)
% scatter(g.X(VALID),g.Y(VALID),[],[100.*aCom.PctVol2100(VALID)],'filled');c=colorbar;caxis([0,200])
% ylabel(c,'Volume at 2100 relative to present (%)')
% subplot(1,2,2)
% histogram(100.*aCom.PctVol2100(VALID),0:5:200);
% xlabel('Volume at 2100 relative to present (%)')
% ylabel('N glaciers')
% 
% saveas(gcf,fullfile(odir,'pctVol2100_summary.png'))
% 
% % aCom.PctComLoss(isnan(aCom.etime))=NaN;
% figure;
% subplot(1,2,1)
% scatter(g.X(ival2),g.Y(ival2),[],100.*Vol2200(1:length(ival2))./VolN(1:length(ival2)),'filled');c=colorbar;caxis([0,200])
% ylabel(c,'Volume at 2200 relative to present (%)')
% subplot(1,2,2)
% histogram(100.*Vol2200(:)./VolN(:),0:5:200);
% xlabel('Volume at 2200 relative to present (%)')
% ylabel('N glaciers')
% 
% saveas(gcf,fullfile(odir,'pctVol2200_summary.png'))
% 
% %% scatter with AAR
% figure;
% scatter(g.AAR(ival),-aCom.PctComLoss(ival),[],g.meanMB(ival),'filled');c=colorbar;caxis([-1,1])
% xlabel('AAR')
% ylabel('Change in ice volume (%)')
% ylabel(c,'Mass balance (m/a)')
% ylim([-100,100])
% 
% saveas(gcf,fullfile(odir,'AAR_pctComLoss_MB.png'))
% 
% % figure;
% % scatter(g.meanMB(ival),-aCom.PctComLoss(ival),'filled');
% 
% figure; scatter(g.meanMB(ival),-aCom.PctComLoss(ival),[],log10(g.Area(ival)));xlim([-2,1.5]);ylim([-100,300])
% figure; scatter(g.Area(ival),aCom.etime(ival))
% figure; scatter(g.meanMB(ival),aCom.etime(ival))
% % binandplot()
% %% regional volume evolution
% ival=find(lig>0);
% aVol=[aCom.Vol;];
% aVolU=[aCom.VolU;];
% % aVol=reshape(aVol,1000,length(ival))';
% pVol=[aCom.Vol./aCom.Vol(:,1);];
% % pVol=reshape(pVol,1000,length(ival))';
% 
% aVol(aVol(:,1)>1e13,:)=NaN; %something wrong with 1 glacier 10^30 m3?!
% 
% 
% fvals=(sum(isfinite(aVol)==0,2))==0;
% rVol=nansum(aVol(fvals,:),1);
% rVolU=sqrt(nansum(aVolU(fvals,:).^2,1));
% 
% figure;
% plot(1:1000,rVol./rVol(1))
% xlabel('Model year')
% ylabel('Regional volume relative to present-day (-)')
% saveas(gcf,fullfile(odir,'regional_volume_3000.png'))
% 
% figure;
% % imagesc(pVol(:,1:200));c=colorbar;caxis([0,2])
% imagesc(pVol(VALID,:));c=colorbar;caxis([0,2])
% xlabel('Model year')
% ylabel('Glacier')
% ylabel(c,'Volume relative to present-day (-)')
% saveas(gcf,fullfile(odir,'glacier_volume_3000.png'))
% 
% 
% %% to 2100
% rpatchx=[1:100];rpatchx=[rpatchx fliplr(rpatchx)];
% rpatchy=repmat(rVol(1:100)./rVol(1),1,2) + [rVolU(1:100)./rVol(1) -rVolU(1:100)./rVol(1)];
% figure;
% a1=subaxis(1,2,2)
% pa=patch(rpatchx,rpatchy);hold on
% pl=plot(1:100,rVol(1:100)./rVol(1))
% xlabel('Model year')
% ylabel('Regional volume relative to present-day (-)')
% grid on
% a1.YAxisLocation='right'
% % saveas(gcf,fullfile(odir,'regional_volume_2100.png'))
% 
% 
% a1=subaxis(1,2,1)
% % figure;
% % imagesc(pVol(:,1:200));c=colorbar;caxis([0,2])
% imagesc(pVol(VALID,1:100));c=colorbar;caxis([0,2])
% xlabel('Model year')
% ylabel('Glacier')
% ylabel(c,'Volume relative to present-day (-)')
% saveas(gcf,fullfile(odir,'changes_2100.png'))
% 
% %% regional MB evolution
% ival=find(lig>0);
% aArea=[aCom.Area];
% % aArea=reshape(aArea,1000,length(ival))';
% aGMB=[aCom.gMB;];
% % aGMB=reshape(aGMB,1000,length(ival))';
% 
% adelV=aArea.*aGMB;
% 
% rArea=nansum(aArea,1);
% rdelV=nansum(adelV,1);
% 
% rMB=rdelV./rArea;
% 
% 
% figure;
% % imagesc(pVol(:,1:200));c=colorbar;caxis([0,2])
% imagesc(aGMB(VALID,:));c=colorbar;caxis([-1,1])
% xlabel('Model year')
% ylabel('Glacier')
% ylabel(c,'glacier-wide mass balance (m w.e. a^{-1})')
% saveas(gcf,fullfile(odir,'glacier_mb_3000.png'))
% 
% figure;
% % imagesc(pVol(:,1:200));c=colorbar;caxis([0,2])
% plot(2001:3000,rMB);
% xlabel('Model year')
% ylabel('glacier-wide mass balance (m w.e. a^{-1})')
% saveas(gcf,fullfile(odir,'region_mb_3000.png'))
% 
% %% regional ablation balance evolution
% ival=find(lig>0);
% aBal=[aCom.Balabl;];
% % aBal=reshape(aBal,1000,length(ival))';
% aImbal=[aCom.Imbalabl;];
% % aImbal=reshape(aImbal,1000,length(ival))';
% aTot=[aCom.Totabl;];
% % aTot=reshape(aTot,1000,length(ival))';
% 
% aPctBal=aBal./aTot;
% aPctImbal=aImbal./aTot;
% 
% rBal=nansum(aBal,1);
% rImbal=nansum(aImbal,1);
% rTot=nansum(aTot,1);
% 
% rPctBal=rBal./rTot;
% rPctImbal=rImbal./rTot;
% 
% 
% figure;
% % imagesc(pVol(:,1:200));c=colorbar;caxis([0,2])
% imagesc(100.*aPctImbal);c=colorbar;caxis([-50,100])
% xlabel('Model year')
% ylabel('Glacier')
% ylabel(c,'Imbalance ablation (% Total ablation)')
% saveas(gcf,fullfile(odir,'glacier_imbalabl_3000.png'))
% 
% figure;
% % imagesc(pVol(:,1:200));c=colorbar;caxis([0,2])
% plot(2001:3000,rPctImbal);
% xlabel('Model year')
% ylabel('Imbalance ablation (% Total ablation)')
% saveas(gcf,fullfile(odir,'region_imbalabl_3000.png'))

%%
save(fullfile(odir,['summary_' date '.mat']))
