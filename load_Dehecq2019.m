function [V] = load_Dehecq2019(Vdata,BBox)

% corner coordinates for bounding box
    bboxCoords.lat = BBox([1,2,2,1,1], 2);
    bboxCoords.lon = BBox([1,1,2,2,1], 1);
    [bboxCoords.xV,bboxCoords.yV] = projfwd(Vdata.gtinfo, bboxCoords.lat, bboxCoords.lon);
    
    iVx = (Vdata.xm>=min(bboxCoords.xV))&(Vdata.xm<=max(bboxCoords.xV));
    iVx = imdilate(iVx,strel('diamond',5)); %expand bounding box by 2 pixels
    iVx1 = find(iVx,1,'first');
    iVxn = sum(iVx);
    iVy = (Vdata.ym>=min(bboxCoords.yV))&(Vdata.ym<=max(bboxCoords.yV));
    iVy = imdilate(iVy,strel('diamond',5)); %expand bounding box by 2 pixels
    iVy1 = find(iVy,1,'first');
    iVyn = sum(iVy);
    Vx = Vdata.xm(iVx);
    Vy = Vdata.ym(iVy);
    [Vxg,Vyg] = meshgrid(Vx,Vy);
    [V.Vlat,V.Vlon] = projinv(Vdata.gtinfo,Vxg,Vyg); %revert to lat/lon
    
    try
    V.Vu = ncread(Vdata.path,'vx',[iVx1 iVy1],[iVxn iVyn])';
    V.Vv = ncread(Vdata.path,'vy',[iVx1 iVy1],[iVxn iVyn])'; %must be transposed - don't know why!!!
    V.Vmag = ncread(Vdata.path,'v',[iVx1 iVy1],[iVxn iVyn])';
    
    V.VuE = ncread(Vdata.path,'vx_err',[iVx1 iVy1],[iVxn iVyn])';
    V.VvE = ncread(Vdata.path,'vy_err',[iVx1 iVy1],[iVxn iVyn])'; %must be transposed - don't know why!!!
    V.VmagE = ncread(Vdata.path,'v_err',[iVx1 iVy1],[iVxn iVyn])';
        
    % correct vector endpoints
    V.e2.x=Vxg+V.Vu; %adjusted vector end-coordinates
    V.e2.y=Vyg+V.Vv; %adjusted vector end-coordinates
    
    [V.e2.lat,V.e2.lon]  = projinv(Vdata.gtinfo,V.e2.x,V.e2.y); %revert to lat/lon
    
    %% clean velocity data - identify areas where err>thresh
    VEthresh = 5;
    
    ERs = (V.VuE>VEthresh)|(V.VvE>VEthresh)|(V.VmagE>VEthresh);
%     ERs = imdilate((V.VuE>VEthresh)|(V.VvE>VEthresh)|(V.VmagE>VEthresh),strel('square',3));
    
    Vu0=V.Vu;
%     V.Vu(ERs) = NaN;
% %     V.Vu = inpaint_nans(V.Vu);
    
    Vv0 = V.Vv;
%     V.Vv(ERs) = NaN;
% %     V.Vv = inpaint_nans(V.Vv);
    
    Vmag0 = V.Vmag;
%     V.Vmag(ERs) = NaN;
% %     V.Vmag = sqrt(V.Vv.^2+V.Vu.^2);
    V.VD=1;
    
    V.ERs=ERs;
    if sum((V.Vmag(:)==-32767)|(V.Vmag(:)==0))==numel(V.Vmag)
        V.VD=0; %empty array where the glacier is
    end
    catch
%         disp(['problem with velocity data'])
        V.VD=0;
    end
    
