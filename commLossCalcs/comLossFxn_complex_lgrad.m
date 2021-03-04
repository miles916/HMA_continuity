function [PctComLoss,etime,a]=comLossFxn_complex_lgrad(DH,DEM,THX,SMB,dx,SMBtype,oname,gifout,P1,P2,P3,P4,H2010)
%P1,P2,P3 are parameters needed for advancing glaciers in the complex version of the Huss parameterization
% P1 is N terminus pixels to average over - best guess is 40
% P2 is partition of mass gain - best guess is 0.3 (30%)
% P3 is the number of pixels to advance into - best guess is 20
% P4 is the terminus longitudinal gradient allowing advance - guess is 20deg

MASK=THX>0;

%% determine Huss param and SMB dist for this glacier - don't fit, use real nDH
gMB0=nanmean(SMB(:)); %m we a-1
if gMB0<-0.1
    [nELs,nDH,~]=invert_Huss2010(DEM,THX,DH);
else
    dEL=10; %10m interval
    ELs = min(DEM(MASK)):dEL:max(DEM(MASK));
    nELs=1-(ELs-min(DEM(MASK)))/(max(DEM(MASK))-min(DEM(MASK)));
end
%determine SMB vs elevation
ELs=min(DEM(:)):10:max(DEM(:));
eSMB=NaN.*ELs;
for iEL=1:length(ELs)-1
    cDEM=(DEM>=ELs(iEL))&(DEM<ELs(iEL+1));
    eSMB(iEL)=nanmean(SMB(cDEM));
%     eArea(iEL)=nansum(cDEM(:));
end

eSMB2=fillmissing(eSMB,'movmean',5);

Tthx=atan(P4*pi()/180)*dx; %calculate terminus thickness threshold based on dx
b(2)=0.007; %mean altitudinal rate across the region - mb gradient varies between 0.007-0.015 - see new results. 0.007 is mean
iEL0=find(isnan(eSMB2)==0,1); %index of lowest elevation with an SMB
b(1)= eSMB2(iEL0)-b(2).*ELs(iEL0); %find y intercept
eSMBint=b(1)+b(2).*ELs; %use mean MB gradient (can be negative slope - leads to indefinite advance)
eSMBint(eSMBint>1)=1; %cap the accumulation
eSMB2(isnan(eSMB2))=eSMBint(isnan(eSMB2));
eSMB2=movmean(eSMB2,10); %smooth further

SMBv = interp1(ELs,eSMB2,DEM(THX>0)); 
SMBel=0.*(THX);SMBel(THX>0)=SMBv; %insert into grid

SMBoffset=SMB-SMBel;SMBoffset(isnan(SMBoffset))=0;

%% set up iteration
% nYr=1000;
nYr=200;
V0=nansum(THX(:)).*dx.^2; %m3 
A0=nansum(THX(:)>0).*dx.^2;%m2
gMB0=nanmean(SMB(:)); %m we a-1
Totabl0=nansum(SMB(SMB<0)).*dx.^2; %m3 a-1
Imbalabl0=A0.*gMB0; %m3 a-1
AAR0=nansum(SMB(:)>0)./nansum(THX(:)>0);
Balabl0=Totabl0-Imbalabl0; %m3 a-1 %balance ablation is total ablation -imbal ablation


iYr=1;
gMB=gMB0;
cV=V0;
cA=A0;
cDEM=DEM;
cDH=DH;
cSMB=SMB;
cTHX=THX;
cMASK=cTHX>0;
dV=((nansum(DH(cMASK)).*dx.^2./cV));
Totabl=Totabl0;
Imbalabl=Imbalabl0;
Balabl=Balabl0;
AAR=AAR0;

    a.gMB=gMB;
    a.Area=cA;
    a.Vol=cV;
    a.dVol=dV;
    a.Yr=2000;
    a.Totabl=Totabl;
    a.Imbalabl=Imbalabl;
    a.Balabl=Balabl;
    a.AAR=AAR;

    cRB=cbrewer('div','RdBu',20);
        if gifout==1
            h=figure('units','normalized','position',[0.15,.2,.8,.4],'Color','w');
            a1=subaxis(1,3,1);
            imagesc(a1,cTHX,'AlphaData',cMASK);c1=colorbar;caxis([0,max(THX(:))]);ylabel(c1,'Ice thickness (m)')
            colormap(a1,'parula')

            axis equal tight
            a2=subaxis(1,3,2);
            imagesc(a2,SMB,'AlphaData',cMASK);c2=colorbar;caxis([-10,10]);ylabel(c2,'SMB (m we/a)')
            axis equal tight
            colormap(a2,cRB)
            sgtitle(num2str(iYr+2000))
            a3=subaxis(1,3,3);
            DHmax=ceil(max(abs(DH(:))));
            imagesc(a3,DH,'AlphaData',cMASK);c2=colorbar;caxis([-DHmax,DHmax]);ylabel(c2,'DH (m/a)')
            axis equal tight
            colormap(a3,(cRB))
            sgtitle(num2str(iYr+2000))
        %     saveas(gcf,[oname '_glacierstate_' num2str(iYr+2000) '.png'])
        end
    %% repeat for nYrs
while iYr<nYr 
    if gifout==1
        imagesc(a1,cTHX,'AlphaData',cMASK);c1=colorbar(a1);caxis(a1,[0,max(THX(:))]);axis(a1,'equal','tight');colormap(a1,'parula');ylabel(c1,'Ice thickness (m)')
        imagesc(a2,cSMB,'AlphaData',cMASK);c2=colorbar(a2);caxis(a2,[-10,10]);axis(a2,'equal','tight');colormap(a2,cRB);ylabel(c2,'SMB (m we/a)')
        DHmax=ceil(max(abs(cDH(:))));
        imagesc(a3,cDH,'AlphaData',cMASK);c2=colorbar(a3);ylabel(c2,'DH (m we/a)');axis(a3,'equal','tight');colormap(a3,(cRB));caxis(a3,[-DHmax,DHmax]);
        sgtitle(num2str(iYr+2000))
        drawnow 
        
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      filename = [oname '_glacierstate.gif'];
      if iYr==1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
      else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
      end
    end
     
      %% determine SMB across the glacier for current geometry
    %option1: determine SMB values based on elevation
    cSMBv = interp1(ELs,eSMB2,cDEM(cMASK),'linear','extrap'); 
    cSMB1=0.*(cMASK);cSMB1(cMASK)=cSMBv; %insert into grid
    
    %option2: use spatial distribution as before
    cSMB2=SMB.*cMASK; 
    
    switch SMBtype
        case 1 %altitudinal SMB
            cSMB=cSMB1.*cMASK;
        case 2 %spatial SMB from before
            cSMB=cSMB1.*cMASK;
            cSMB(MASK)=cSMB2(MASK);
        case 3 %mixture
            cSMB=(cSMB1+cSMB2).*cMASK/2;
        case 4 %offset approach
            cSMB=(cSMB1+SMBoffset.*cMASK);
    end
    gMB=nanmean(cSMB(cMASK));
    if(gMB0<-0.1) && H2010==2
        [cDH,~]=apply_Huss2010_complex(nELs,nDH,gMB,cDEM,cTHX,P1,P2,P3,Tthx);
    else
        [cDH,~]=apply_Huss2010_complex(nELs,cA,gMB,cDEM,cTHX,P1,P2,P3,Tthx);
    end
    cTHX(cMASK==0)=0; 
    cTHXt=cTHX+cDH;%temporary variable
    t1=cTHXt<0;
    cDH(t1)=-cTHX(t1); %identify areas with no thickness left

    t2=cSMB>0; %areas with positive SMB cannot become non-glacier through thinning (would cut off mass supply unrealistically)
    cDH(t2&(cTHXt<=0))=0;%nominally ensure that positive SMB areas remain glacier, with thickness of 1m...
    
    %rescale cDH to conserve mass gain (loss) after above corrections
    cDH(~(t1|t2))=cDH(~(t1|t2)).*gMB./nanmean(cDH(cMASK&~(t1|t2)));
    
    %then correct the results given the new DH
    cTHX=cTHX+cDH;
    cDEM=cDEM+cDH;
    cMASK=cTHX>0;

    iYr=iYr+1;
    
    cA=nansum(cMASK(:)).*dx.^2; 
    cV=nansum(cTHX(:)).*dx.^2; %current thickness
    dV=(gMB.*cA)./cV;
    Totabl=nansum(cSMB(cSMB<0)).*dx.^2; %m3 a-1
    Imbalabl=nanmin(cA.*gMB,0); %m3 a-1
    AAR=nansum(cSMB(:)>0)./nansum(cTHX(:)>0);
    Balabl=Totabl-Imbalabl; %m3 a-1 %balance ablation is total ablation -imbal ablation

    
    %display the run
    disp([num2str(iYr) ' ' num2str(gMB,3) ' ' num2str(Imbalabl./Totabl,3) ' ' num2str(cV,3) ' ' num2str(dV,4)])

    %save outputs
    a.gMB(iYr)=gMB;
    a.Area(iYr)=cA;
    a.Vol(iYr)=cV;
    a.dVol(iYr)=dV;
    a.Yr(iYr)=iYr+2000;
    a.Totabl(iYr)=Totabl;
    a.Imbalabl(iYr)=Imbalabl;
    a.Balabl(iYr)=Balabl;
    a.AAR(iYr)=AAR;

end
    
%% normalize relative to start
a.pArea=a.Area./A0;
a.pVol=a.Vol./V0;

%% e-folding time
a.delV=V0-a.Vol;
a.pdelV=a.delV./a.delV(end);
[~,etime]=min(abs(a.pdelV-(1-exp(-1))))
PctComLoss=100*(1-a.Vol(end)./a.Vol(1));

%%
save([oname '.mat'],'a','PctComLoss','etime')