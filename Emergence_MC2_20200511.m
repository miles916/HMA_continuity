clear 
close all

RAVENpath = '\\wsl.ch\fe\gebhyd\8_Him';
[codedir,~,~]=fileparts(mfilename('fullpath'));
addpath(genpath(codedir))
% addpath(genpath([RAVENpath '\Personal_folders\Evan\EmergenceSMB\code']))

RGIR='15'; %RGI region to process
X3 = 2; %minimum glacier size (km2) to consider
plotouts=0; %output plots or not
exports=1; %save geotiffs or not
val=1; %validation sites or all

THXopt='c'; %[1:5,'c']

% outdir = [RAVENpath '\Personal_folders\Evan\EmergenceSMB\RGI131415_ALL2km_' date]
outdir = [RAVENpath '\Personal_folders\Evan\EmergenceSMB\ValidationSites_' date]
mkdir(outdir)

%% prepare velocity loading info data (Dehecq2019)
V.path = [RAVENpath '\Remote_Sensing_Data\Regional_Global\Dehecq2019_ITS_LIVE\HMA_G0240_0000.nc'];
V.info = ncinfo(V.path);
V.Rdata = (str2num(V.info.Variables(13).Attributes(7).Value))

V.gtinfo = geotiffinfo([RAVENpath '\Personal_folders\Evan\EmergenceSMB\NP_HMA_G0240_vx.tif']); %easiest way to get correct projection details
V.R = V.Rdata([5,6;2,3;1,4]);% V.EPSG = V.info.Variables(13).Attributes(10).Value;
V.xp = [1:V.info.Dimensions(1).Length]-0.5;
V.yp = [1:V.info.Dimensions(2).Length]-0.5; %pixel coordinates are cell centers
[V.xm,~] = pix2map(V.R,ones(size(V.xp)),V.xp);
[~,V.ym] = pix2map(V.R,V.yp,ones(size(V.yp)));

%% load RGI for region
RGIdir = [RAVENpath '\Remote_Sensing_Data\Regional_Global\00_rgi60'];
cd(RGIdir);
RGIRdir = dir([RGIR '*']);
RGIRdir=RGIRdir([RGIRdir(:).isdir]==1);
cd(RGIRdir.name)
RGIRfile = dir([RGIR '*.shp']);
RGI=shaperead(RGIRfile.name);

if val==1 %only run for validation sites
if str2double(RGIR)==15
% igs = 4770;
    igs = [...
        2291, ... %Ganju La % NOT IN BRUN DATA (or Dehecq!)!!
        3416, ... %pokalde
        3448, ... %trakarding
        3507, ... %??
        3585, 3584, 3586, 3581, 3583, 3592, ... %Mera/Naulek % 
        3733,3473,3417,3742,3743,9991,10055,... %Everest region
        3734, ... %Changri
        3954, ... %Yala
        4045,4121,4119,4176,9474, ... %Langtang valley
        4770, ... %Milarepa
        4847, ... %Rikha Samba
        7122,6881,6777, ... %Satopanth, Gangotri, Dunagiri 
        7143, ... %Chorabari
        7605, ... %Dokriani
        7645, ... %Baishui % not in Velocity data
        7886, ... %Hailuogou
        9026, ... %Naimon'nyi
        10263, ... %Kangwure
        11758,11752,11909,11963,11693,11962,12707, ... %Parlung area
        11973, ... %Parlung N4
        ];
elseif str2double(RGIR)==14
    igs = [...
        04291,... %Virjerab     
        06794, ... %Baltoro
        11566,... %Naradu
        15536,... %Hamtah
        15990,... %Chhota Sigri % VALIDATION SITES
        16041, ... %Sutri Dhaka
        16042, ... %Batal Dhaka
        19394,... %Shaigiri
        20156,... %Rakhiot
        ];
elseif str2double(RGIR)==13
    igs = [...
%             05000, ... %Inylchek
            33853,...
            19758, ...%unknown problem glacier in pamirs!
            06361, ... %Karak-batkak
            06974, ... %Batysh Sook
            07064, ... %No354 - Ak-Shiyrak
            08054, ... %Bordu N
            08055, ... %Sary-Tor
            08624, ... %TsTuyuksuyskiy
            10093, ... %No599 
            11609, ... %Golubin
            18096, ... %Abramov
            28839, ... %Gurenhekou
            24602, ... %Xiaodongkemadi
            32330, ... %Qiyi
            41896, ... %Mustagh Ata
            43165, ... %72
            43174, ... %74
            43207, ... %Tuomer
            43232, ... %Koxcar
            45335, 45336, ... %Urumqi E and W
            49645, ... %Zhongxi
            49754, ... %Zhadang
        ];
else
    error('wrong domain')
end
    %to index to correct RGI outline (needed for zone 13 especially)
    RGIid = zeros(length(RGI),1); 
    for ir = 1:length(RGI)
        RGIid(ir) = str2double(RGI(ir).RGIId(end-4:end));
    end

else
    %trim to glaciers over X km2??
    RGI=RGI([RGI(:).Area]>X3);

    %to index to correct RGI outline (needed for zone 13 especially)
    RGIid = zeros(length(RGI),1); 
    for ir = 1:length(RGI)
        RGIid(ir) = str2double(RGI(ir).RGIId(end-4:end));
    end

    igs = RGIid;
end

%% load Scherler debris for region
DEBdir = [RAVENpath '\Remote_Sensing_Data\Regional_Global\Scherler2018_global_debris\S2_2015-2017_NDSI'];
cd(DEBdir);
DEBfile = dir([RGIR '*.shp']);
schDEB=shaperead(DEBfile.name);
dRGIid = zeros(length(schDEB),1);
for ir = 1:length(schDEB)
    dRGIid(ir) = str2double(schDEB(ir).RGIId(end-4:end));
end
%% working directory for region
Fdir0 = [RAVENpath '\Remote_Sensing_Data\Regional_Global\Farinotti2019_composite_thickness_RGI60-all_regions'];

Hdir = fullfile(Fdir0,['RGI60-' RGIR]);

%% create parallel pool if desired

% delete(gcp('nocreate'))
% % 
% poolobj = parpool(10)

%% iterate through RGI

profile on
% parfor ig = 1:length(igs) %for parallel
for ig = 1:1%length(igs) 
    delete a m %intermediate variables from MC analysis
     disp([ig length(igs)])
    
    cd(outdir)

%  basic glacier properties
    cur=RGI(RGIid==igs(ig));
    cur.RGIID = cur.RGIId(7:end);
        
    cur.LL = [cur.Y; cur.X]';
    cur.LL(isnan(cur.LL(:,1)),:)=[];
    cur.LL = unique(floor(cur.LL),'rows');
    cur.utmZone= utmzone(cur.CenLat,cur.CenLon);
    
    if (cur.Area)<15
        dx=50;
    elseif (cur.Area)>80
        dx=200;
    else
        dx=100;
    end
        outtitle1 = [num2str(dx) 'mgrid'];

    %% load DEM
    buffdist = 1000; %area around glacier to extract
    [cur.DEMall,cur.DEMallLon,cur.DEMallLat] = load_GDEM3_v2(RAVENpath,cur.LL,cur.BoundingBox,buffdist);
    
   %% read, subset,preprocess velocity
    cur.V = load_Dehecq2019(V,cur.BoundingBox);
        try
                cur.V.VmagE(cur.V.VmagE(:)==0)=NaN;
                cur.paramsVEL=lognfit(cur.V.VmagE(~isnan(cur.V.VmagE(:))&~cur.V.ERs(:)));
                cur.V.SEraw=logninv(0.68,cur.paramsVEL(1),cur.paramsVEL(2));
                vGd=(cur.V.ERs==0);
        catch
            cur.V.VD=0; %don't break the loop. if no data, skip the rest 
        end
        
        if cur.V.VD==1
  
    %% load dH (Brun2017)
    [cur.dH,cur.dHE,cur.dHlat,cur.dHlon,cur.DHE] = load_Brun2017(RAVENpath,cur.LL,cur.BoundingBox); 
%         paramsDH=lognfit(cur.dHE(cur.dHE>0))
%         ecdf(cur.dHE0(:));hold on
        
%         plot(0.01:.01:1,logncdf(0.01:.01:1,paramsDH(1),paramsDH(2)),'k')
%         cur.dHE0=cur.dHE;
% %         cur.dHE=logninv(0.68,paramsDH(1),paramsDH(2)); %this is the systematic uncertainty
% %         cur.dHE=prctile(cur.dHE0(:),68); %this is the systematic uncertainty
%         cur.dHE=rms(cur.dHE); %this is the systematic uncertainty
        
%%
    if cur.DHE>0
        %% load Farinotti 2019 data
        
        switch THXopt
            case 'c' %composite case
                try
                    cur.Ffile = fullfile(Hdir,[cur.RGIId '_thickness.tif']);
                    [cur.F.thxR,cur.R3] = geotiffread(cur.Ffile);
                    [cur.R3b] = geotiffinfo(cur.Ffile);
                catch %drive connection sometimes fails - wait and try again
                    pause(10)
                    cur.Ffile = fullfile(Hdir,[cur.RGIId '_thickness.tif']);
                    [cur.F.thxR,cur.R3] = geotiffread(cur.Ffile);
                    [cur.R3b] = geotiffinfo(cur.Ffile);
                end
            otherwise %any option model
%                 if str2num(THXopt)>=1 & str2num(THXopt)<=5
                    h= '\\wsl.ch\fe\gebhyd\1_Alle\3_Projekte\HIMAL\Remote_Sensing_Data\Regional_Global\Farinotti2019_allModels';

                    cur.Ffile = fullfile(h,['results_model_' THXopt '/RGI60-' RGIR '/thickness_' cur.RGIId '.tif']);
                    if (exist(cur.Ffile,'file')==2)
                        [cur.F.thxR,cur.R3] = geotiffread(cur.Ffile);
                        [cur.R3b] = geotiffinfo(cur.Ffile);
                    end
                
        end
        
       if (exist(cur.Ffile,'file')==2)
     
        cur.F.thxerrR = Farin2019_uncertainty(cur.RGIId,cur.F.thxR); %loads individual model results and calculates error as std
        cur.F.thxerrR(cur.F.thxerrR==0)=NaN;
        %thx error model
        paramsTHX=lognfit(cur.F.thxerrR(~isnan(cur.F.thxerrR(:))))
        cur.F.thxerrR0=cur.F.thxerrR;
        cur.F.thxerrR=logninv(0.68,paramsTHX(1),paramsTHX(2)); %this is the systematic uncertainty
        
        %% REPROJECT ALL (to Farinotti CS)

        %coordinates from Farinotti
        cur.F.x3=cur.R3.XWorldLimits(1)+cur.R3.CellExtentInWorldX.*[cur.R3.XIntrinsicLimits(1):cur.R3.XIntrinsicLimits(2)-1];
        cur.F.y3=cur.R3.YWorldLimits(2)-cur.R3.CellExtentInWorldY.*[cur.R3.YIntrinsicLimits(1):cur.R3.YIntrinsicLimits(2)-1];
        [cur.F.x3g,cur.F.y3g]=meshgrid(cur.F.x3,cur.F.y3);

        %new coordinates based on Farinotti, but at specified interval
        x3 = dx.*[(floor((cur.F.x3(1)-buffdist)/dx)):((ceil(cur.F.x3(end)+buffdist)/dx))];
        y3 = dx.*[(ceil((cur.F.y3(1)+buffdist)/dx)):-1:(floor((cur.F.y3(end)-buffdist)/dx))];
        [x3g,y3g] = meshgrid(x3,y3);

        Rout = [0,-dx;dx,0;x3(1),y3(1)];

        %RESAMPLE THICKNESS
        cur.F.thx = griddata(cur.F.x3g(:),cur.F.y3g(:),double(cur.F.thxR(:)),x3g(:),y3g(:),'cubic');
        cur.F.thx = reshape(cur.F.thx,size(x3g));
        cur.F.thxerr=repmat(cur.F.thxerrR,size(x3g));

        %resample velocity data
        [cur.V.x0,cur.V.y0] = projfwd(cur.R3b,cur.V.Vlat,cur.V.Vlon); %compute projected coordinates
        [cur.V.e2.x0,cur.V.e2.y0] = projfwd(cur.R3b,cur.V.e2.lat,cur.V.e2.lon); %compute projected coordinates
        cur.e2.Vu=cur.V.e2.x0-cur.V.x0; %project vector ends
        cur.e2.Vv=cur.V.e2.y0-cur.V.y0;%project vector ends
        cur.V.U = griddata(cur.V.x0(:),cur.V.y0(:),double(cur.e2.Vu(:)),x3g(:),y3g(:),'cubic');
        cur.V.U = reshape(cur.V.U,size(x3g));
        cur.V.UE = griddata(cur.V.x0(:),cur.V.y0(:),double(cur.V.VuE(:)),x3g(:),y3g(:),'cubic');
        cur.V.UE = reshape(cur.V.UE,size(x3g));
        cur.V.V = griddata(cur.V.x0(:),cur.V.y0(:),double(cur.e2.Vv(:)),x3g(:),y3g(:),'cubic');
        cur.V.V = reshape(cur.V.V,size(x3g));
        cur.V.VE = griddata(cur.V.x0(:),cur.V.y0(:),double(cur.V.VvE(:)),x3g(:),y3g(:),'cubic');
        cur.V.VE = reshape(cur.V.VE,size(x3g));
        cur.V.UE = repmat(cur.V.SEraw,size(x3g));
        cur.V.VE = repmat(cur.V.SEraw,size(x3g));
        cur.V.S = sqrt(cur.V.V.^2+cur.V.U.^2);
        
        %---------no longer used --------------
        cur.V.ERRORs = griddata(cur.V.x0(:),cur.V.y0(:),double(cur.V.ERs(:)),x3g(:),y3g(:),'nearest'); %get all affected pixels
        cur.V.ERRORs = reshape(cur.V.ERRORs,size(x3g));
        cur.V.ERRORs = imdilate(cur.V.ERRORs,strel('square',3));

        %resample AW3D
        [cur.DEMallX,cur.DEMallY] = projfwd(cur.R3b,cur.DEMallLat,cur.DEMallLon); %compute projected coordinates
        cur.DEM = griddata(cur.DEMallX(:),cur.DEMallY(:),double(cur.DEMall(:)),x3g(:),y3g(:),'cubic');
        cur.DEM = reshape(cur.DEM,size(x3g));

        %resample dH
        [cur.dHx0,cur.dHy0] = projfwd(cur.R3b,cur.dHlat,cur.dHlon); %compute projected coordinates
        cur.dH2 = griddata(cur.dHx0(:),cur.dHy0(:),double(cur.dH(:)),x3g(:),y3g(:),'cubic');
        cur.dH2 = reshape(cur.dH2,size(x3g));
        cur.dHE2 = griddata(cur.dHx0(:),cur.dHy0(:),double(cur.dHE(:)),x3g(:),y3g(:),'cubic');
        cur.dHE2 = reshape(cur.dHE2,size(x3g));

        %% load Scherler data, reproject, set mask, debris mask
        cur.mask = cur.F.thx>0;
        cur.DEB = schDEB(dRGIid==igs(ig));
        if sum(dRGIid==igs(ig))>0
            cur.debrismask = load_proc_Scherler2018(cur.DEB,cur.R3b,cur.mask,Rout);
        else
            cur.debrismask= 0.*cur.mask;
        end


        cur.dH3 = cur.dH2;
%         final inputs ... set masked pixels to 0. NaN will break the
%         boundary calculations with regards to mass conservation
        cur.dH3((cur.mask==0))=0;
        cur.dHE2((cur.mask==0))=0;
        cur.dH3(abs(cur.dH3)>50)=0;
        cur.V.U((cur.mask==0))=0;
        cur.V.V((cur.mask==0))=0;
        cur.V.UE((cur.mask==0))=0;
        cur.V.VE((cur.mask==0))=0;
        cur.F.thx((cur.mask==0))=0;

%     % column-averaged velocity values
    [umult,sig_umult,f] = f_probability2(cur.F.thx(cur.mask)); % uniform distribution for u_b/u_s but gives normal dist out
    paramsUMult=fitdist(f','normal');

    cur.V.Smean=umult.*cur.V.S;
    cur.V.Umean=umult.*cur.V.U;
    cur.V.Vmean=umult.*cur.V.V;
    

%%  1000 random draws for systematic (THXscale,umult,Vscale) and random (10m thx, sig m/a vel, sig density & dH) uncertainty
    nr= 1000;
    rns=randn(nr,3); %values
    prns=cdf('Normal',rns,1,1);
        
    % initialize vars
    tDIVFx=zeros(size(cur.F.thx));
    tDIVFy=zeros(size(cur.F.thx));
    tDIVF=zeros([size(cur.F.thx),nr]);
    tDH=zeros([size(cur.F.thx),nr]);
    
    dx = mode(diff(x3));
    dy = mode(diff(y3));
    
%     loop through runs
    tic
    for irun=1:nr
%         determine perturbations to data
        u_mult=icdf(paramsUMult,prns(irun,1));
        THXmult=1+icdf('lognormal',(prns(irun,2)),paramsTHX(1),paramsTHX(2)).*(rns(irun,2)./abs(rns(irun,2))); %extra to determine scaling up or down
        Vmult=1+icdf('lognormal',(prns(irun,3)),cur.paramsVEL(1),cur.paramsVEL(2)).*(rns(irun,3)./abs(rns(irun,3))); %extra to determine scaling up or down

        tUmean=cur.V.U.*u_mult.*Vmult;
        tVmean=cur.V.V.*u_mult.*Vmult;
        tTHX=cur.F.thx.*THXmult+10.*randn(size(cur.F.thx)).*cur.mask;

        tDIVFx(:,2:end-1) = (tUmean(:,1:end-2).*tTHX(:,1:end-2)-tUmean(:,3:end).*tTHX(:,3:end)).*abs(dy)/2;
        tDIVFy(2:end-1,:) = (-tVmean(1:end-2,:).*tTHX(1:end-2,:)+tVmean(3:end,:).*tTHX(3:end,:)).*abs(dx)/2; %note y coords are reversed
        tDIVF(:,:,irun) = (tDIVFx+tDIVFy)/abs(dy*dx);%total and normalize to area, m/yr
        tDH(:,:,irun) = cur.dH3+cur.dHE2.*randn(size(cur.F.thx)).*cur.mask;%total and normalize to area, m/yr
        
    end
    
    Qdensity = 0.9+0.05.*randn(size(cur.mask,1),size(cur.mask,2),nr).*cur.mask; %900 kg m3 everywhere plus random perturbations of 50 kg/m3
    Hdensity = repmat(0.9.*cur.mask,1,1,nr); %initial value everywhere of 900kg m3
    ind1=(tDIVF>0);
    ind2=(tDH>0);
    ind3=(abs(tDH)>abs(tDIVF));
    
    Hdensity(ind1&~ind2)=0.9;
    Hdensity(~ind1&ind2)=0.6;
    Hdensity(ind1&ind2&ind3)=0.6;
    Hdensity(ind1&ind2&~ind3)=0.85;
    Hdensity(~ind1&~ind2&ind3)=0.9;
    Hdensity(~ind1&~ind2&~ind3)=0.85;
    Hdensity = Hdensity+0.05.*randn(size(cur.mask,1),size(cur.mask,2),nr).*cur.mask; %no longer 900 kg m3 everywhere, add random perturbations of 50 kg/m3

    SMBstack1 = cur.mask.*(Hdensity.*tDH-Qdensity.*tDIVF);
    %%temporary structures to hold MC outputs
    a=struct;
    m=struct;
    a.SMBmean = nanmean(SMBstack1,3);
    a.SMBsig = nanstd(SMBstack1,[],3);
    a.accstack1=SMBstack1>0;
    a.accfreq1=nanmean(a.accstack1,3);
    a.accsig1=nanstd(a.accstack1,[],3);
    
    %% secondary outputs from MC
    
    cur.DIVF=nanmean(tDIVF,3); cur.DIVF(cur.mask==0)=NaN;
    cur.DIVFe=nanstd(tDIVF,[],3); cur.DIVFeS(cur.mask==0)=NaN; %absolute magnitude of error
       
    % density correction
    cur.Qdensity = cur.mask.*nanmean(Qdensity,3); %900 kg m3 everywhere
    cur.Hdensity = cur.mask.*nanmean(Hdensity,3); %initial value everywhere of 900kg m3
    sig_HDens= nanstd(Hdensity,[],3);
    sig_QDens= nanstd(Qdensity,[],3);
    
    % SMB and uncertainty
    cur.V.Smean=cur.V.S.*umult;
    cov_H_S=cov(cur.F.thx(cur.mask),cur.V.Smean(cur.mask)); %covariance matrix
    
    
    cur.SMB = a.SMBmean;
    cur.SMB(cur.mask==0)=NaN;
    cur.SMBe=a.SMBsig;

    %% metrics from MC
    [m.ELA,m.ELAacc,EstABL]= calc_ELA_errmin(cur.DEM,a.SMBmean,cur.mask); %ELA from mean data
    m.AAR=nansum(cur.DEM(cur.mask)>m.ELA)./nansum(cur.mask(:)); %AAR from mean data
    a.AARestMC=squeeze(nansum(a.accstack1,[1,2])./sum(cur.mask(:)));
    a.AARsigMC=std(a.AARestMC);
    
    %ELA sigma from AAR sigma
    els=min(cur.DEM(cur.mask)):10:max(cur.DEM(cur.mask));
    hyps=0.*els;
    for iEL=1:length(els)
        hyps(iEL)=nansum(cur.DEM(cur.mask)>els(iEL))./nansum(cur.mask(:));
    end
    [hyps, ih] = unique(hyps);
    els=els(ih);
    a.ELA_l=interp1(hyps,els,m.AAR-a.AARsigMC);
    a.ELA_u=interp1(hyps,els,m.AAR+a.AARsigMC);
    m.sig_ELA = (a.ELA_l-a.ELA_u)/2;
    
    %ablation balance
    a.totalabl=squeeze(nansum(SMBstack1.*(SMBstack1<0),[1,2]));
    a.imbalabl=squeeze(nansum(SMBstack1,[1,2]));
    a.balabl=a.totalabl-a.imbalabl;
    a.bal2tot=a.balabl./a.totalabl;
    
    m.totalabl = nansum(a.SMBmean(a.SMBmean<0)); %PIXELS
    m.imbalabl = nansum(a.SMBmean(:)); %PIXELS
    m.balabl=m.totalabl-m.imbalabl;
    m.bal2tot=m.balabl./m.totalabl;
    
%%     mask before plotting
    cd(outdir)
    cur.dH3((cur.mask==0))=NaN;
    cur.V.U((cur.mask==0))=NaN;
    cur.V.V((cur.mask==0))=NaN;
    cur.F.thx((cur.mask==0))=NaN;
    cur.DIVF((cur.mask==0))=NaN;
    cur.DIVFx((cur.mask==0))=NaN;
    cur.DIVFy((cur.mask==0))=NaN;
    cur.SMB((cur.mask==0))=NaN;
    cur.SMBh2((cur.mask==0))=NaN;

% aggregate to zones
    cur.zones2 = (segment_Gmask_EL_debris(cur.DEM,cur.mask,cur.debrismask));
    
%     dH
    cur.z2DH = zonal_aggregate_v2(cur.zones2,cur.dH3,'mean'); %simply aggregates values in the zone - FAST
    cur.z2DHe = zonal_aggregate_v2(cur.zones2,abs(cur.dHE2.*cur.dH3),'rms'); %simply aggregates values in the zone to mean value

%     fdiv
    cur.z2fdiv = zonal_aggregate_v2(cur.zones2,cur.DIVF,'mean'); %simply aggregates values in the zone - FAST
    cur.z2fdivE = zonal_aggregate_v2(cur.zones2,cur.DIVFe,'rms'); %simply aggregates values in the zone to mean value

%     smb
    cur.z2smb = zonal_aggregate_v2(cur.zones2,cur.SMB,'mean'); %simply aggregates values in the zone - FAST
    cur.z2smbE = zonal_aggregate_v2(cur.zones2,cur.SMBe,'rms'); %simply aggregates values in the zone to mean value

    %
if plotouts==1
cmap = [0,0,0;cbrewer('div','RdBu',21);0,0,0];

    curout = 'grid-emergence'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(cur.DIVF)
    colorbar;colormap(flipud(cmap))
    title(curtitle)
    caxis([-5,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    curout = 'EL-zones'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(cur.zones2)
    colorbar;colormap([0,0,0;lines(256)])
    title(curtitle)
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    curout = 'grid-emergence-unc'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc((cur.DIVFe))
    colorbar;colormap(flipud(gray))
    title(curtitle)
    caxis([0,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    curout = 'grid-Hdensity'
    curtitle = [curout ' (1000 km m^{-3}), ' outtitle1];
    figure
    imagesc((cur.Hdensity))
    colorbar;
    title(curtitle)
    caxis([.5,1])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    curout = 'grid-SMB-SNR'
    curtitle = [curout ' (m w.e. a^{-1}), ' outtitle1];
    figure
    imagesc((abs(cur.SMB./cur.SMBe)));c=colorbar;
    title(curtitle)
    caxis([0.5,20])
    set(gca,'colorscale','log')
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    
    curout = 'EL-zone-avg-emergence'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(cur.z2fdiv)
    colorbar;colormap(flipud(cmap))
    title(curtitle)
    caxis([-5,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])

    
    curout = 'EL-zone-avg-SMB-SNR'
    curtitle = [curout ' (m w.e. a^{-1}), ' outtitle1];
    figure
    imagesc((abs(cur.z2smb./cur.z2smbE)));c=colorbar;
    title(curtitle)
    caxis([0.5,20])
    set(gca,'colorscale','log')
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])

    curout = 'EL-zone-SMB-unc'
    curtitle = [curout ' (m w.e. a^{-1}), ' outtitle1];
    figure
    imagesc((cur.z2fdivE))
    colorbar;colormap(flipud(gray))
    title(curtitle)
    caxis([0,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    curout = 'grid-dh'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(cur.dH3)
    colorbar;colormap(cmap)
    title(curtitle)
    caxis([-5,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])

    curout = 'grid-SMB'
    curtitle = [curout ' (m w.e. a^{-1}), ' outtitle1];
    figure
    imagesc(cur.SMB)
    colorbar;colormap(cmap)
    title(curtitle)
    caxis([-5,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    curout = 'grid-SMB-unc'
    curtitle = [curout ' (m w.e. a^{-1}), ' outtitle1];
    figure
    imagesc(cur.SMBe)
    colorbar;%colormap(cmap)
    title(curtitle)
    caxis([0,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
    
    curout = 'thx'
    curtitle = [curout ' (m), ' outtitle1];
    figure
    imagesc(cur.F.thx)
    colorbar;
    title(curtitle)
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])

    
    curout = 'EL-zonal-SMB'
    curtitle = [curout ' (m w.e. a^{-1}), ' outtitle1];
    figure
    imagesc(cur.z2smb)
    colorbar;colormap(cmap)
    title(curtitle)
    caxis([-5,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])

    
    curout = 'EL-zonal-SMB-unc'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(cur.z2smbE)
    colorbar;colormap(flipud(gray))
    title(curtitle)
    caxis([0,5])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])

    
    curout = 'surface-speed'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(cur.V.Smean)
    colorbar;
    title(curtitle)
    caxis([0,20])
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])
        
    
    curout = 'SMB-acczones'
    curtitle = [curout ' (m a^{-1}), ' outtitle1];
    figure;
    imagesc(a.accfreq1)
    c= colorbar;
    caxis([0,1])
    ylabel(c,'SMB>0 freq')
    saveas(gcf,[cur.RGIID '_' curout '_' outtitle1 '.png'])

end
    %% export 
if exports==1
%save output of MC data
    savetofile({a,m},[cur.RGIID '_data.mat'])

%     Debris
    geotiffwrite([cur.RGIID '_debris.tif'],uint8(cur.debrismask+cur.mask),Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     DEM
    geotiffwrite([cur.RGIID '_AW3D.tif'],uint16(cur.DEM),Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     THX
    geotiffwrite([cur.RGIID '_THX.tif'],cur.F.thx,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     VEL
    geotiffwrite([cur.RGIID '_Smean.tif'],cur.V.Smean,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     dH
    geotiffwrite([cur.RGIID '_dH.tif'],cur.dH3,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % DIVF
    geotiffwrite([cur.RGIID '_FDIV.tif'],cur.DIVF,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)])

    % SMB
    geotiffwrite([cur.RGIID '_SMB.tif'],cur.SMB,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)])

%     dHe
    geotiffwrite([cur.RGIID '_dHe.tif'],cur.dHE2,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % DIVFe
    geotiffwrite([cur.RGIID '_FDIVe.tif'],cur.DIVFe,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)])

    % SMBe
    geotiffwrite([cur.RGIID '_SMBe.tif'],cur.SMBe,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)])

%     Zones
    geotiffwrite([cur.RGIID '_zones.tif'],uint16(cur.zones2),Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))
    
%     Hdensity
    geotiffwrite([cur.RGIID '_Hdensity.tif'],uint16(1000.*cur.Hdensity),Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     zDIVF
    geotiffwrite([cur.RGIID '_zFDIV.tif'],cur.z2fdiv,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     zSMB
    geotiffwrite([cur.RGIID '_zSMB.tif'],cur.z2smb,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     zDIVFe
    geotiffwrite([cur.RGIID '_zFDIVe.tif'],cur.z2fdivE,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

%     zSMBe
    geotiffwrite([cur.RGIID '_zSMBe.tif'],cur.z2smbE,Rout,'CoordRefSysCode',['EPSG:' num2str(cur.R3b.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

end

     else
        disp(['No thx data for ' cur.RGIID])
       end
    else
        disp(['No dH data for ' cur.RGIID])
    end
    else
        disp(['No velocity data for ' cur.RGIID])
    end
    close all

end
%%
cd(outdir)

disp('finished')

% delete(poolobj);
