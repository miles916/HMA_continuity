function thxerr = Farin2019_uncertainty(RGIId,thx)

% h= '\\wsl.ch\fe\gebhyd\1_Alle\3_Projekte\HIMAL\Remote_Sensing_Data\Regional_Global\Farinotti2019_allModels';
h= '\\wsl.ch\fe\gebhyd\8_Him\Remote_Sensing_Data\Regional_Global\Farinotti2019_allModels';

RGIR=RGIId(7:8);

allthx=NaN([size(thx),5]);

T1n = fullfile(h,['results_model_1/RGI60-' RGIR '/thickness_' RGIId '.tif']);
if (exist(T1n))==2
    allthx(:,:,1)=imread(T1n);
end

T2n = fullfile(h,['results_model_2/RGI60-' RGIR '/thickness_' RGIId '.tif']);
if (exist(T2n))==2
    allthx(:,:,2)=imread(T2n);
end

T3n = fullfile(h,['results_model_3/RGI60-' RGIR '/thickness_-' RGIId '.tif']);
if (exist(T3n))==2
    allthx(:,:,3)=imread(T3n);
end

T4n = fullfile(h,['results_model_4/RGI60-' RGIR '/thickness_' RGIId '.tif']);
if (exist(T4n))==2
    allthx(:,:,4)=imread(T4n);
end

T5n = fullfile(h,['results_model_5/RGI60-' RGIR '/thickness_' RGIId '.tif']);
if (exist(T5n))==2
    allthx(:,:,5)=imread(T5n);
end

thxstd=nanstd(allthx,[],3);
thxrng=nanmax(allthx,[],3)-nanmin(allthx,[],3);
thxct=nansum(allthx>0,3);

thxstd(thxct==2)=thxrng(thxct==2)./sqrt(2);
thxstd(thxct<2)=nanmean(thxstd(:));

thxse=thxstd./sqrt(thxct);
thxerr =thxse./thx; % normalise to estimated value for percent error
thxerr(isinf(thxerr))=NaN;

