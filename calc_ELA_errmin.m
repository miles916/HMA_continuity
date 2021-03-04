function [ELA,acc,estABL,ELvals,Accuracy] = calc_ELA_errmin(DEM,SMB,MASK)

% a1=figure;
% scatter(DEM(MASK),SMB(MASK));
% xlabel('ELevation (m)')
% ylabel('SMB (m/a)')
% % a1=gca;

    % el-bin segmentation
%     ELbins= (50*floor(min(cur.F.AW3Del(iN))/50)):50:(50*ceil(max(cur.F.AW3Del(iN))/50));
    ELvals=unique(floor(DEM(MASK)));
    TAB = NaN(size(ELvals));TAC = TAB;
    FAC =TAB;FAB = TAB;
    
    ACC=SMB>0;
    ABL=SMB<0;
    
    estABL=0; %assume we can derive an ELA
% go through all values and determine error rates for each
% 
N =sum(~isnan(SMB(:)));

    for iEL=1:length(ELvals)
        cBelow = ((DEM<ELvals(iEL))&MASK);
        cAbove = ((DEM>=ELvals(iEL))&MASK);
        valBelow(iEL) = nanmean(SMB(cBelow));
        valAbove(iEL) = nanmean(SMB(cAbove));
        
        cTAC= cAbove&ACC;
        cFAC= cAbove&ABL;
        cTAB= cBelow&ABL;
        cFAB= cBelow&ACC;
        
        TAC(iEL)=sum(cTAC(:));
        FAC(iEL)=sum(cFAC(:));
        TAB(iEL)=sum(cTAB(:));
        FAB(iEL)=sum(cFAB(:));
    end
    
    Accuracy = (TAC+TAB)./(N);
    Precision = TAC./(TAC+FAC);
    Recall= TAB./(TAB+FAB);
    Dice = 2*Precision.*Recall./(Precision+Recall);
    
    Accuracy=Accuracy.*double((valAbove>-0.1)&(valBelow<0.10))';
    [~,imax]=max(Accuracy);
    
    [~,idx,w,p] = findpeaks(smooth(Accuracy,50)); %possible spots
%     [~,idx,w,p] = findpeaks(smooth(Dice,50)); %possible spots
%     Above = NaN(size(idx));Below=Above;
%     for ix=1:length(idx)
    Above=valAbove(idx);
    Below=valBelow(idx);
    ABdiff=Above-Below; %want to maximize this
    
    poss = (Above>-0.1)&(Below<0.10);
%     pks=pks(poss);
%     idx=idx(poss);
%     w=w(poss);
%     p=p(poss);
    if (sum(poss)>0)%&&(imax>1)
    
    %     plot(ELvals(idx(poss)),Accuracy(idx(poss)),'+b')
%         idxF=idx(p==max(p(poss)));
        idxF=idx(ABdiff==max(ABdiff(poss)));

        if numel(idxF)>1
            [~,ix2]=max(Accuracy(idxF));
            idxF=idxF(ix2);
        end

%         figure;
%         plotyy(ELvals,smooth(Accuracy,50),ELvals,smooth(Dice,50));hold on
%         plot(ELvals(idx),Accuracy(idx),'xr')
%         plot(ELvals(idxF),Accuracy(idxF),'sk')
%         xlabel('Eevation (m)')
%         ylabel('Accuracy (L) / Dice (R)')
%         a2=gca;

        ELA=ELvals(idxF);
        acc=Accuracy(idxF);
    else
        
%         if nanmean(SMB(MASK))<0
%             ELA=max(DEM(MASK));
            disp(['ELA above ' num2str(max(DEM(MASK)))])
            %limit data to above mean elevation in case of debris-cover
            mEL=nanmean(DEM(MASK));
            SEL=(DEM>mEL) & MASK;
            y=SMB(SEL);
            x=[ones(sum(SEL(:)),1),DEM(SEL)];
%             scatter(DEM(MASK),SMB(MASK));
%             lsline;
            [b,~,~,~,stats]=regress(y,x);
%             if (stats(1)>0.5)&(b(1)<0)&(stats(3)<0.01)
            if (b(2)>0)
                ELA=-b(1)./b(2);
                estABL=stats(3);
                [~,idxF] = min(abs(ELvals-ELA));
%                 acc=Accuracy(idxF);
                acc=stats(1);
            else
                estABL=b(1);
                ELA=NaN;acc=NaN;
            end
%         else
%             ELA=min(DEM(MASK));
%             ELA=NaN;acc=NaN;
%             disp(['ELA below ' num2str(min(DEM(MASK)))])
%             estABL=-1;
%         end
    end