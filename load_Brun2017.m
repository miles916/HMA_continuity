function [dH,dHE,dHlat,dHlon,DHE] = load_Brun2017(RAVENpath,LL,BBox)
   
    %initialize
    DHE=0;%0 DEM tiles loaded
    dH=[];
    dHE=[];
    dHlon=[];
    dHlat=[];
        
    bboxCoords.lat = BBox([1,2,2,1,1], 2);
    bboxCoords.lon = BBox([1,1,2,2,1], 1);
    
    for Itile = 1:size(LL,1)
        
    % load dH tiles
    dH0dir = [RAVENpath '\Remote_Sensing_Data\Regional_Global\Brun2017_ASTER_shared_HMA'];
    dH1dir = ['n' num2str(LL(Itile,1),'%02d') '_e' num2str(LL(Itile,2),'%03d')];
    dH2file = ['dh_dt_2000-2016_ASTER_' dH1dir '.tif'];
    fdHpath = fullfile(dH0dir,dH1dir,dH2file);
    %err
    dH2Efile = ['dh_dt_2000-2016_ASTER_' dH1dir '_err.tif'];
    fdHEpath = fullfile(dH0dir,dH1dir,dH2Efile);
    
    if exist(fdHpath,'file')>0
%         [dHraw] = geotiffread(fdHpath);
%         [dHEraw] = geotiffread(fdHEpath);
        [R2] = geotiffinfo(fdHpath);
        dHx1 = ([1:R2.Width]-0.5).*R2.RefMatrix(2,1)+R2.RefMatrix(3,1);
        dHy1 = ([1:R2.Height]-0.5).*R2.RefMatrix(1,2)+R2.RefMatrix(3,2);

        %project bounding box
        [bboxCoords.xH,bboxCoords.yH] = projfwd(R2, bboxCoords.lat, bboxCoords.lon);

        idHx = (dHx1>=min(bboxCoords.xH))&(dHx1<=max(bboxCoords.xH));
        idHx = imdilate(idHx,strel('diamond',5)); %expand bounding box by 2 pixels
        idHx1 = find(idHx,1,'first');
        idHxn = sum(idHx);
        idHy = (dHy1>=min(bboxCoords.yH))&(dHy1<=max(bboxCoords.yH));
        idHy = imdilate(idHy,strel('diamond',5)); %expand bounding box by 2 pixels
        idHy1 = find(idHy,1,'first');
        idHyn = sum(idHy);
        dHx1 = dHx1(idHx);
        dHy1 = dHy1(idHy);
        [dHxg,dHyg] = meshgrid(dHx1,dHy1);
        [dHlat1,dHlon1] = projinv(R2,dHxg,dHyg);

        %load only relevant area
        ySt = find(idHy,1,'first');
        yEn = find(idHy,1,'last');
        xSt = find(idHx,1,'first');
        xEn = find(idHx,1,'last');
        dHraw2=imread(fdHpath,'PixelRegion',{[ySt,yEn],[xSt,xEn]});
        dHEraw2=imread(fdHEpath,'PixelRegion',{[ySt,yEn],[xSt,xEn]});

%         dHraw2 = dHraw(idHy,idHx); %subset
%         dHEraw2 = abs(dHEraw(idHy,idHx)); %subset errors

        %    imagesc(cur.DEM2x,cur.DEM2y,cur.DEM2);hold on; line(cur.X,cur.Y);set(gca,'Ydir','normal')
        % combine tiles
            dH = [dH; dHraw2(:)];
            dHE = [dHE; dHEraw2(:)];
            dHlat = [dHlat; dHlat1(:)];
            dHlon = [dHlon; dHlon1(:)];
        dH(abs(dH)>30)=NaN;
%         dH(dHE>1)=NaN; %no need?
        %remove nans
        dHlat(isnan(dH))=[];
        dHlon(isnan(dH))=[];
        dHE(isnan(dH))=[];
        dH(isnan(dH))=[];
        DHE = DHE+1;%a dH exists
    if numel(dHlat)<20 %possible that all pixels are NaN - no dH data
        DHE = 0;
    end
    end
    end