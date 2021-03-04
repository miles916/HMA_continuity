function debrismask = load_proc_Scherler2018(DEB,R,mask,Rout)

%load Scherler and process to consistent Dmask

%remove small rings less than 40 nodes
          partends = find(isnan(DEB.X));
          partlengths = diff([0,partends]);
          ikeep = find(partlengths>40);

          DEB.X2=[];
          DEB.Y2=[];

          for irem=1:length(ikeep)
              DEB.X2=[DEB.X2,DEB.X(partends(ikeep(irem))-partlengths(ikeep(irem))+1:partends(ikeep(irem)))];
              DEB.Y2=[DEB.Y2,DEB.Y(partends(ikeep(irem))-partlengths(ikeep(irem))+1:partends(ikeep(irem)))];
          end

        %   scatter(cur.DEB.X,cur.DEB.Y);hold on
        %   scatter(cur.DEB.X2,cur.DEB.Y2,'r')

          [DEB.Xn,DEB.Yn] = projfwd(R,DEB.Y2,DEB.X2); %compute projected coordinates from LL
          DEB.Xn=DEB.Xn+50; %one-pixel shift due to center vs corner coordinates
          DEB.Yn=DEB.Yn-50; %one-pixel shift due to center vs corner coordinates
          [DEB.r,DEB.c] = map2pix(Rout,DEB.Xn,DEB.Yn); %convert to pixel coordinates
      

    partends2 = find(isnan(DEB.c));
    partlengths2 = diff([0,partends2]);

    debrismask0 = 0.*mask;

    for ip=1:length(partends2)
        debrismask = poly2mask(DEB.c(partends2(ip)-partlengths2(ip)+1:partends2(ip)-1),DEB.r(partends2(ip)-partlengths2(ip)+1:partends2(ip)-1),size(mask,1),size(mask,2));
        debrismask0 = debrismask0|debrismask;
    end

    %trim slivers and holes
    nodebris = (mask-debrismask0).*mask;
    nodebris2 = bwareaopen(nodebris,5);

    %smooth and again remove small holes
    nodebris3 =medfilt2(nodebris2,[5,5]).*mask;
    % debrismask2 = medfilt2(debrismask,[5,5]);

    debrismask=mask-nodebris3;
    debrismask = bwareaopen(debrismask,10);
    debrismask = imfill(debrismask,'holes');
