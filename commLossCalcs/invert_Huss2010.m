function [nELs,nDH,Nparams] = invert_Huss2010(DEM,THX,DH)
% function [nELs,nDH,Nparams] = invert_Huss2010(DEM,THX,DH)

%% 
MASK=THX>0;

%% determine hypsometry
gArea=sum(MASK(:)); %just npixels, to normalize area!

dEL=10; %10m interval
ELs = min(DEM(MASK)):dEL:max(DEM(MASK));
pArea = NaN.*ELs;
pDH=pArea;
pTHX = pArea;

for iel=1:numel(ELs)
    cur=(DEM<ELs(iel)+dEL/2)&(DEM>=ELs(iel)-dEL/2); %current section of DEM
    pArea(iel)=sum(cur(:))./gArea;
    pTHX(iel)=nanmean(THX(cur));
    pDH(iel)=nanmean(DH(cur));
end

% figure
% yyaxis left
% plot(ELs,pArea);
% ylabel('Portion of area')
% yyaxis right
% plot(ELs,pTHX);
% ylabel('Glacier thickness')
% saveas(gcf,'glacier_geometry.png')

%% normalize and filter
nELs=1-(ELs-min(DEM(MASK)))/(max(DEM(MASK))-min(DEM(MASK)));

pDH2=fillmissing(pDH,'movmean',10,'EndValues','nearest');
pDH2(pDH2>0)=0; %eliminates thickening
pDH2=movmean(pDH2,30);
pDH2=fillmissing(pDH2,'linear','EndValues','extrap');
nDH=pDH2./nanmin(pDH2);
% nDH(nDH<0)=0; %eliminates thickening
nDH(end-5:end)=0; %limits headwall retreat

% figure
% subaxis(2,1,1)
% plot(ELs,pDH);hold on
% plot(ELs,pDH2);hold on
% ax1=gca;
% ylabel('Elevation change (m a~{-1})')
% xlabel('Elevation (m)')
% set(ax1,'Xdir','reverse')
% grid on
% legend('observed','smoothed')
% 
% ax2=subaxis(2,1,2)
% % ax2=axes('Position',ax1.Position,'Color','none');
% plot(ax2,nELs,nDH);
% % ax2.Color='none';
% ylabel('Normalized elevation change (m a~{-1})')
% xlabel('Normalized elevation')
% % % set(ax2,'XAxisLocation','top','YAxisLocation','right')
% % set(ax2,'YAxisLocation','right')
% set(ax2,'Ydir','reverse')
% % ylim([0.05,1])
% grid on
% saveas(gcf,'normalized_DH.png')

%% determine parameterization
%H2010 parameter sets
params = [-.3,.6,.09,2];
paramm = [-.05,.19,.01,4];
paraml = [-.02,.12,0,6];

% delHs = @(a,b,c,d) (nELs+a).^d+b.*(nELs+a)+c; %according to Huss 2010
% ft=fittype( @(a,b,c,d,x) (x+a).^round(d)+b.*(x+a)+c); %rounding to ensure coefficient is an integer!!
ft=fittype( @(a,b,c,d,x) (max(x+a,0)).^(d)+b.*(x+a)+c); %ensure X+a>=0!!
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-.5,0.1,0,0],...
               'Upper',[-0.01,1,0.15,7]);%,...
%                'StartPoint', paraml;
f1 = fit(nELs', nDH', ft, fo); %start with 'med' param set

%now the same but with d rounded to an integer
ft2=fittype( @(a,b,c,x) (max(x+a,0)).^round(f1.d)+b.*(x+a)+c); %ensure X+a>=0!!
f2 = fit(nELs', nDH', ft2, fo); %start with 'med' param set

%% compare with H2010 coeffs
           
delHs = @(hr) (hr+params(1)).^params(4)+params(2).*(hr+params(1))+params(3); %according to Huss 2010
delHm = @(hr) (hr+paramm(1)).^paramm(4)+paramm(2).*(hr+paramm(1))+paramm(3); %according to Huss 2010
delHl = @(hr) (hr+paraml(1)).^(paraml(4))+paraml(2).*(hr+paraml(1))+paraml(3); %according to Huss 2010
% delHNEW = @(hr) (hr+f1.a).^round(f1.d)+f1.b.*(hr+f1.a)+f1.c; %according to Huss 2010
% delHNEW = @(hr) (max(hr+f1.a,0)).^round(f1.d)+f1.b.*(hr+f1.a)+f1.c; %according to Huss 2010
delHNEW = @(hr) (max(hr+f2.a,0)).^round(f1.d)+f2.b.*(hr+f2.a)+f2.c; %according to Huss 2010


% figure
% plot(nELs,nDH);hold on;
% plot(nELs,delHs(nELs));hold on;
% plot(nELs,delHm(nELs));hold on;
% plot(nELs,delHl(nELs));hold on;
% plot(nELs,delHNEW(nELs));hold on;
% legend('local','small','med','large','fitted','Location','west')
% ylabel('Normalized elevation change (m a~{-1})')
% xlabel('Normalized elevation')
% set(gca,'Ydir','reverse')
% grid on
% ylim([-0.020,1])
% saveas(gcf,'normalized_DH2.png')

Nparams= [f2.a,f2.b,round(f2.c,2),round(f1.d)];
Neq=['delH = (hr+' num2str(f2.a,2) ')^' num2str(round(f1.d)) ' + ' num2str(f2.b,2) '*(hr+' num2str(f2.a,2) ') + ' num2str(round(f2.c,2))]

