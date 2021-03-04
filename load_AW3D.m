function [DEMall,DEMallX,DEMallY] = load_AW3D(RAVENpath,LL,BBox,buffdist)

    DEMall=[];DEMallX=[];DEMallY=[];
%     if size(cur.LL,1)>1 %ONLY RERUNNING GLACIERS!!!!
    for Itile = 1:size(LL,1)
        
    % load AW3D tiles
    DEM0dir = [RAVENpath '\Remote_Sensing_Data\Regional_Global\AW3D'];
    DEM1dir = ['N' num2str(LL(Itile,1),'%03d') 'E' num2str(LL(Itile,2),'%03d')];
    DEM2file = [DEM1dir '_AVE_DSM.tif'];
    fDEMpath = fullfile(DEM0dir,DEM1dir,DEM2file);
    
    [DEM,R] = geotiffread(fDEMpath);
    DEMx = R.LongitudeLimits(1)+R.CellExtentInLongitude.*([R.XIntrinsicLimits(1):R.XIntrinsicLimits(2)-1]);
    DEMy = R.LatitudeLimits(2)-R.CellExtentInLatitude.*([R.YIntrinsicLimits(1):R.YIntrinsicLimits(2)-1]);


    %determine appropriate buffer
%         buffdist = 2000;
        %calc lat distance at latitude - approx
        dlat = (buffdist/111000); %111km per degree latitude
        %calc lon distance at latitude
        CenLat=mean(BBox(:,2));
        dlon = (buffdist/(111000*cos(CenLat*pi()/180))); %111km per degree longtitude at equator

    
    % extract area of glacier
    ix = (DEMx>=BBox(1,1)-dlon) & (DEMx<=BBox(2,1)+dlon);
    ix = imdilate(ix,strel('square',64)); %grow by a few pixels
    iy = (DEMy>=BBox(1,2)-dlat) & (DEMy<=BBox(2,2)+dlat);
    iy = imdilate(iy,strel('square',64)); %grow by a few pixels
    DEM2 = DEM(iy,ix);
    DEM2x = DEMx(ix);
    DEM2y = DEMy(iy);
    
%    imagesc(cur.DEM2x,cur.DEM2y,cur.DEM2);hold on; line(cur.X,cur.Y);set(gca,'Ydir','normal')
%     clear cur.DEM cur.DEMx cur.DEMy
% .DEM=[];cur.DEMx=[];cur.DEMy=[];
    
    % project DEM
    [DEM2lon,DEM2lat]=meshgrid(DEM2x,DEM2y);
    [DEM3x,DEM3y]=ll2utm(DEM2lat,DEM2lon);
    
    % combine DEMS
%     if Itile == 1
%         cur.DEMall = cur.DEM2(:);
%         cur.DEMallX = cur.DEM3x(:);
%         cur.DEMallY = cur.DEM3y(:);
%     else
        DEMall = [DEMall; DEM2(:)];
        DEMallX = [DEMallX; DEM3x(:)];
        DEMallY = [DEMallY; DEM3y(:)];
%     end
%         clear cur.DEM2 cur.DEM2x cur.DEM2y cur.DEM3x cur.DEM3y cur.DEM2lon cur.DEM2lat
%         DEM2=[];DEM2x=[];DEM2y=[];DEM3x=[];DEM3y=[];cur.DEM2lon=[];cur.DEM2lat=[];
    end
    