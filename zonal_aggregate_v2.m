function [zVAR] = zonal_aggregate_v2(zones,VAR,method)

mask = isnan(zones);

zVAR = zeros(size(VAR));

zonelist = unique(zones(:));

switch method
    case 'mean'
        for iz = 1:length(zonelist)
            cz = zonelist(iz);
            curzone = (zones==cz);
            cval = nanmean(VAR(curzone(:)==1));
            zVAR(curzone(:)==1)=(cval);
        end
    case 'rms'
        for iz = 1:length(zonelist)
            cz = zonelist(iz);
            curzone = (zones==cz);
            cval = rms(VAR(curzone(:)==1));
            zVAR(curzone(:)==1)=(cval);
        end
    case 'rssn'
        for iz = 1:length(zonelist)
            cz = zonelist(iz);
            curzone = (zones==cz);
            select=(curzone)&(isnan(VAR)==0);
            cval = sqrt(nansum(VAR(select).^2))./sum(select(:)); %root-sum-squared/N
            zVAR(curzone(:)==1)=(cval);
        end
end

zVAR(mask) = NaN;