function [newDHg,newDEM] = apply_Huss2010_complex(nELs,input1,gMB,DEM,THX,P1,P2,P3,P4)
% P1 is N terminus pixels to average over - best guess is 40
% P2 is partition of mass gain - best guess is 0.3 (30%)
% P3 is the number of pixels to advance into - best guess is 20
% P4 is the terminus thickness allowing advance - guess is 10m
% [newDHg,newDEM] = apply_Huss2010(nELs,nDH,gMB,DEM,THX,dx)
% [newDHg,newDEM] = apply_Huss2010(nELs,Area,gMB,DEM,THX,dx)
% [newDHg,newDEM] = apply_Huss2010(nELs,[param],gMB,DEM,THX,dx)
if nansum(THX(:))>0
    MASK=THX>0;

    if numel(input1)>4
        nDH=input1;
        idealized=0;
    elseif numel(input1)>1
        param=input1;
        idealized=0;
    else
        Area=input1;
        idealized=1;
    end

    rho_g = 0.85; %glacier-wide specific gravity of ice at the glacier surface (see Huss et al 2013);

    %% normalize elevation and determine hypsometry
    nDEM=(DEM-nanmin(DEM(MASK)))./(nanmax(DEM(MASK))-nanmin(DEM(MASK)));
    nDEM=(1-nDEM); %Huss2010 uses inverted normalized elevation!
    nDEM(MASK==0)=NaN;

    gArea=sum(MASK(:)); %just npixels, to normalize area!

    nArea = NaN.*nELs;

    for iel=1:numel(nELs)-1
        cur=(nDEM>nELs(iel+1))&(nDEM<=nELs(iel)); %current section of DEM
        nArea(iel)=nansum(cur(:))./gArea;
    end

    % figure
    % yyaxis left
    % plot(nELs,nArea);
    % ylabel('Portion of area')
    % yyaxis right
    % plot(nELs,nTHX);
    % ylabel('Glacier thickness')
    % saveas(gcf,'normalized_geometry.png')


    %% determine fs and estimate new surface
    if idealized==1
        % determine parameters based on size
        if exist('param')==0
            if Area<5e6 %m2
                param = [-.3,.6,.09,2];
            elseif Area<20e6
                param = [-.05,.19,.01,4];
            else
                param = [-.02,.12,0,6];
            end
        end

        %apply equation
        delH = @(hr) (hr+param(1)).^param(4)+param(2).*(hr+param(1))+param(3); %according to Huss 2010

        relDH = delH(nELs);
        fs=gMB./(rho_g.*nansum(relDH.*nArea));
        newDH=fs.*relDH; %altitudinal
    else
        fs=gMB./(rho_g.*nansum(nDH.*nArea));
        newDH=fs.*nDH; %altitudinal
    end

    %%
    % figure
    % plot(nELs,newDH)
    % if exist('nDH')
    %     hold on
    %     plot(nELs,nDH)
    %     legend('Huss2010 DH','obs DH')
    % else
    %     legend('Huss2010 DH')
    % end
    % ylabel('elevation change (m a~{-1}')
    % xlabel('Normalized elevation (inverted)')    
    % saveas(gcf,'DH_profile_compare.png')

    %% distribute spatially
    newDHv=interp1(nELs,newDH,nDEM(MASK));
    newDHg=double(MASK);newDHg(MASK)=newDHv;

    newDHg(isnan(newDHg))=0;
    newDHg(MASK==0)=0;

    % figure; 
    % subaxis(1,2,1)
    % plot(nELs,newDH)
    % subaxis(1,2,2)
    % imagesc(newDHg);colorbar
    %% deal with positive mass balance; similar to Huss and Hock 2015 but distinct routine for when a glacier advances: terminus must be 10m thick

        %determine thickness of lowest 10 pixels - this is the thickness to be added
        [~,il10]=sort(DEM(MASK));THXv=THX(MASK);ie = min(length(THXv),P1);
        l10thx=nanmean(THXv(il10(1:ie)));

        if l10thx>P4 %if the terminus is thicker than 10m and MB is positive, then advance. otherwise just thicken
            advVol=min(P2.*gMB.*nansum(THX(:)>0),l10thx.*P3); %volume attirbutable to advance is 1/10 the total mass gain or terminus thickness for 20 adjacent pixels; note in pix-m
            gMBn=gMB-advVol./(nansum(THX(:)>0)); %remaining mass gain to distribute; note advVol is in pix-m
            newDHg=newDHg.*gMBn./gMB; %scale the DH pattern accordingly

            %identify terminus pixels as lowest-elevation P3 pixels bordering glacier
            PotentialPix = logical(imdilate(MASK,strel('diamond',1))-MASK);
            [~,iRank]=sort(DEM(PotentialPix)); %the index to the ranked values
            [~,Rank]=sort(iRank); %the rank of each value
            rankedPotentialPix=0.*DEM;rankedPotentialPix(PotentialPix)=Rank;
            Advance=(PotentialPix &(rankedPotentialPix<=P3));

            newDHg(Advance)=advVol./P3; %distribute volume attributable to advance 

        end

    newDEM=DEM+newDHg;
else
    newDEM=DEM;
    newDHg=0.*THX;
end