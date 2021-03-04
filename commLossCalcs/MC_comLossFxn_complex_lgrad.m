function MC_comLossFxn_complex_lgrad(SMB,SMBe,adj,nit,DH,DEM,THX,dx,SMBtype,oname,gifout,P1,P2,P3,P4,H2010)
        try
        for it=1:nit
            SMBc=SMB+nanmean(SMBe(:)).*adj(it); %perturb SMB by error for the run
            evalc('[PctComLoss(it),e_time(it),a(it)]=comLossFxn_complex_lgrad(DH,DEM,THX,SMBc,dx,SMBtype,oname,gifout,P1,P2,P3,P4,H2010)'); %suppress outputs
        end
        M.PctComloss=nanmean(PctComLoss);
        U.PctComloss=nanstd(PctComLoss);
        M.etime=nanmean(e_time);
        U.etime=nanstd(e_time);
        M.gMB=nanmean(reshape([a.gMB],1000,nit)',1);
        U.gMB=nanstd(reshape([a.gMB],1000,nit)',1);
        M.Area=nanmean(reshape([a.Area],1000,nit)',1);
        U.Area=nanstd(reshape([a.Area],1000,nit)',1);
        M.Vol=nanmean(reshape([a.Vol],1000,nit)',1);
        U.Vol=nanstd(reshape([a.Vol],1000,nit)',1);
        M.Totabl=nanmean(reshape([a.Totabl],1000,nit)',1);
        U.Totabl=nanstd(reshape([a.Totabl],1000,nit)',1);
        M.Balabl=nanmean(reshape([a.Balabl],1000,nit)',1);
        U.Balabl=nanstd(reshape([a.Balabl],1000,nit)',1);
        M.Imbalabl=nanmean(reshape([a.Imbalabl],1000,nit)',1);
        U.Imbalabl=nanstd(reshape([a.Imbalabl],1000,nit)',1);
        M.AAR=nanmean(reshape([a.AAR],1000,nit)',1);
        U.AAR=nanstd(reshape([a.AAR],1000,nit)',1);
        M.pVol=nanmean(reshape([a.pVol],1000,nit)',1);
        U.pVol=nanstd(reshape([a.pVol],1000,nit)',1);
        M.pArea=nanmean(reshape([a.pArea],1000,nit)',1);
        U.pArea=nanstd(reshape([a.pArea],1000,nit)',1);
        save([oname '.mat'],'M','U')
        
        %% plots
%     % final plots
%     Yr=2001:3000;
%         figure
%         yyaxis left
%         p=plot(Yr,M.pVol.*100,'-','LineWidth',1.5);hold on 
%         h=area(Yr,[M.pVol-U.pVol;2.*U.pVol]'.*100,'FaceColor','flat','FaceAlpha',0.3,'LineStyle','none');
%         h(1).FaceColor=[1,1,1];h(1).FaceAlpha=0;
%         h(2).FaceColor=p.Color;
% % plot(Yr,(M.pVol-U.pVol).*100,':');hold on 
% %         plot(Yr,(M.pVol+U.pVol).*100,':');hold on 
% %         plot(Yr,M.pArea.*100,'-r');hold on 
% %         plot(Yr,(M.pArea-U.pArea).*100,':r');hold on 
% %         plot(Yr,(M.pArea+U.pArea).*100,':r');hold on 
%         ylabel('Volume (% of present-day)')
% %         ylim([0,max([100,100.*a.pVol(end),100.*a.pArea(end)])])
%         yyaxis right
%         p2=plot(Yr,M.gMB,'-','LineWidth',1.5);hold on 
%         h2=area(Yr,[M.gMB-U.gMB;2*U.gMB]','FaceColor','flat','FaceAlpha',0.3,'LineStyle','none');
%         h2(1).FaceColor=[1,1,1];h2(1).FaceAlpha=0;
%         h2(2).FaceColor=p2.Color;
% %         plot(Yr,0.*Yr,'-k');
% %         plot(Yr,(M.gMB-U.gMB),':');hold on 
% %         plot(Yr,(M.gMB+U.gMB),':');hold on 
% %         legend([p,p2],'V','MB','Location','east')
%         ylabel('Mass balance (m w.e. a^{-1})')
%         xlim([min(Yr)-0.5,max(Yr)+0.5])
%         saveas(gcf,[oname '_glacier_evolution.png'])
% 
% 
%         %% plot ablbal
% 
%         % final plots
%         figure
%         yyaxis left
%         a1=area(Yr,-1.*[M.Balabl;M.Imbalabl]','FaceColor','flat','FaceAlpha',0.5,'LineStyle','none');hold on
%         p=line(Yr,-1.*M.Totabl,'Linewidth',1.5);
%         h=area(Yr,-1.*[M.Totabl-U.Totabl;2*U.Totabl]','FaceColor','flat','FaceAlpha',0.3,'LineStyle','none');
%         h(1).FaceColor=[1,1,1];h(1).FaceAlpha=0;
%         h(2).FaceColor=p.Color;        
%         h2=area(Yr,-1.*[M.Balabl-U.Balabl;2*U.Balabl]','FaceColor',[0,0,0],'FaceAlpha',0.3,'LineStyle','none');
%         h2(1).FaceAlpha=0;
% %         line(Yr,-1.*(M.Balabl-U.Balabl),'Linestyle',':');
% %         line(Yr,-1.*(M.Balabl+U.Balabl),'Linestyle',':');
% %         line(Yr,-1.*(M.Totabl-U.Totabl),'Linestyle',':');
% %         line(Yr,-1.*(M.Totabl+U.Totabl),'Linestyle',':');
%         ylabel('Ablation (m^3 w.e. a^{-1})')
%         % ylim([0,100])
%         yyaxis right
%         p2=plot(Yr,M.AAR,'-','Linewidth',1.5);hold on 
%         h3=area(Yr,[M.AAR-U.AAR;2*U.AAR]','FaceColor','flat','FaceAlpha',0.3,'LineStyle','none');
%         h3(1).FaceColor=[1,1,1];h3(1).FaceAlpha=0;
%         h3(2).FaceColor=p2.Color;        
%         legend([a1(1),a1(2),p,p2],'Balance','Imbalance','Total','AAR','Location','north')
%         ylabel('AAR (-)')
%         xlim([min(Yr)-0.5,max(Yr)+0.5])
%         saveas(gcf,[oname '_glacier_evolution2.png'])
        end