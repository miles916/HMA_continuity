function zones = segment_Gmask_EL_debris(DEM,mask,debrismask)
% DEM is DEM for whole area
% mask denotes glacier area

DEM0=DEM;
    if sum(mask(:))>500
    DEM1 = DEM0;
    DEM1(mask==0)=NaN;

    DEM2 = inpaint_nans(DEM1,4);

    DEM3 = imgaussfilt(DEM,2);
    DEM3 = imgaussfilt(DEM3,2);

    % mask = isnan(DEM);

    DEM = DEM3;
end
DEM(mask==0)=NaN;

% [SLO,DIR] = imgradient(DEM./dx);
% [Fx,Fy] = gradient(DEM3);
% SLO = sqrt(Fx.^2+Fy.^2)/dx;
% SLO(mask==0)=NaN;

% %% elevation range 
% 
ELmax = nanmax(DEM(mask));
ELmin = nanmin(DEM(mask));

Cint = 50;
ELvs = Cint.*[floor(ELmin/Cint):1:ceil(ELmax/Cint)];

zones = zeros(size(DEM));
%% iterate through zones
iZ = 0;

for iEL=2:length(ELvs)
    cur = (DEM>ELvs(iEL-1))&(DEM<=ELvs(iEL));
    cur = imfill(cur,'holes');
    cz = bwlabel(cur); %label each segment individually
    cn = max(cz(:)); %number of new zones
    zones = zones+cz+cur.*iZ;
    iZ = iZ+cn;
    imagesc(zones.*double(mask)); colorbar;
end

% zones=remove_small_zones(uint16(zones),mask,30);

zones = (zones+(iZ+1).*double(debrismask)).*double(mask);

figure;
imagesc(zones.*double(mask)); colorbar;
axis square
% set(gca,'ydir','normal')
colormap([0,0,0;rand(256,3)])
% % caxis([0,length(zonelist)])
% % linkaxes([a1,a2],'xy')

