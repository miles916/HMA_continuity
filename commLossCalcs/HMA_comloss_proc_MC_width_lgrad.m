%% go to folder, load files, paths
clear
clc
close all

diary log

addpath(genpath('\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Huss2010_flowparam'))
addpath(('\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Committed_loss_calcs'))
addpath(genpath('\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\code'))
addpath(('\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\debris_fluxes'))

%set path to analyze, load list of files
hdir = '\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\RGI131415_ALL2km_13-May-2020'; 

%load input data
glacier_postproc_results = '\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\EmergenceSMB\POST_2km2_15-May-2020\'
X=2;
cd(glacier_postproc_results)
res=dir('ELA*.mat')
load(res(end).name);
VALID2 = ~g.issurging & g.ELAacc>0.5 & g.ELAs>g.Zmin & (g.DHstd)<3 & (g.MBstd_C)<3 &(g.MBmeansig<5) & (g.DHslope'>-0.5 );
VALID=VALID2 &g.Area>X;

%set output directory
odir = ['\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\Committed_loss_calcs\ComLoss_SMB2_tgrad_width_MC_' date]

%settings
gifout=0;
SMBtype=2; %1 uses fixed altitudinal profile; 2 uses fixed spatial distribution, 3 is mean of the first 2; 4 is a fixed pixel-based deviation to the altitudinal profile, plus the altitudinal profile
N=5; %domain expansion for positive MB glaciers (km) so that growth is not domain-limited
P2 = 0.5;% P2 is partition of mass gain - best guess is ~0.483 (48%) but equifinality with P4
P4 = 11.8; % P4 is the terminus longitudinal gradient allowing advance - best guess is 11.8deg with P2=50%
H2010=2; %configuration of flow-parameterization: 1 is with area-based curves (as in H2010) and 2 is with each glacier's DH curve

% go to output folder
mkdir(odir)
cd(odir)

set(0, 'DefaultFigureWindowStyle', 'docked');

%% set up MC
nit=30;
adj=randn(nit,1);


delete(gcp('nocreate'))
 
poolobj = parpool(10);

ival=find(VALID);
%%
tic
parfor iig=1:length(ival)
    ig=ival(iig);
        disp([iig length(ival) ig])
        close all
        oname=g.RGIid{ig}(end-7:end);
        if exist(fullfile(odir,[oname '.mat']),'file')==0 %skip if results already exist
            DH=imread(fullfile(hdir,[oname,'_dH.tif']));
            DHe=imread(fullfile(hdir,[oname,'_dHe.tif']));
            DEM=double(imread(fullfile(hdir,[oname,'_AW3D.tif'])));
            THX=imread(fullfile(hdir,[oname,'_THX.tif']));
            SMB=imread(fullfile(hdir,[oname,'_zSMB.tif']));
%             SMBe=imread(fullfile(hdir,[oname,'_zSMBe.tif']));
            info=geotiffinfo(fullfile(hdir,[oname,'_dH.tif']));
            dx=info.PixelScale(1);

%             try

            %measure glacier width of terminus
            mask=THX>0;
            widths=measureTerminusWidths(DEM,mask,dx);close
            P1 = mean(widths)./dx;% P1 is N terminus pixels to average over - determine as glacier width
            P3 = P1;% P3 is the number of pixels to advance into - best guess is equalt to P1

                %% extend domain for positive MB glaciers
                if (nanmean(SMB(:))>0)
                    DEM0=DEM;
                    DH0=DH;
                    THX0=THX;
                    SMB0=SMB;
                    SMBe0=DHe;%use DHe as it represents GMB uncertainty well

                    % extract DEM for expanded area, reproject
                    BBox0=info.BoundingBox + N.*1e3.*[-1,-1;1,1]; %extend bounding box
                    BBoxUTMx=BBox0([1,1,2,2,1]);
                    BBoxUTMy=BBox0([3,4,4,3,3]);
                    [BBoxLLy,BBoxLLx]=projinv(info,BBoxUTMx,BBoxUTMy);
                    BBox=[min(BBoxLLx) min(BBoxLLy);max(BBoxLLx) max(BBoxLLy)];
                    LL = [BBoxLLy; BBoxLLx]';
                    LL(isnan(LL(:,1)),:)=[];
                    LL = unique(floor(LL),'rows');

                    [DEMall,DEMallLon,DEMallLat] = load_GDEM3_v2('\\wsl.ch\fe\gebhyd\8_Him',LL,BBox,0);
                    [DEMallX,DEMallY] = projfwd(info,DEMallLat,DEMallLon); %compute projected coordinates

                    %setup new refmat, extending N km
                    Refmat=info.RefMatrix;
                    Refmat(3,1)=Refmat(3,1)-N*1e3;
                    Refmat(3,2)=Refmat(3,2)+N*1e3;
                    offxy=N.*1e3./Refmat(2,1); %Npixels to extend

                    %interpolate to new grid
                    [xg,yg]=pixcenters(Refmat,size(THX)+2*offxy,'makegrid');
                    DEM = griddata(DEMallX(:),DEMallY(:),double(DEMall(:)),xg(:),yg(:),'cubic');
                    DEM = reshape(DEM,size(xg));

                    %insert everything into new large domain
                    ins=zeros(size(DEM));
                    ins(offxy+1:offxy+size(THX,1),offxy+1:offxy+size(THX,2))=1;
                    ins=logical(ins);
                    DH=NaN.*DEM;DH(ins)=DH0;
                    THX=NaN.*DEM;THX(ins)=THX0;
                    SMB=NaN.*DEM;SMB(ins)=SMB0;
                    SMBe=NaN.*DEM;SMBe(ins)=SMBe0;

                end
                
                %the main calculation....
                
            MC_comLossFxn_complex_lgrad(SMB,DHe,adj,nit,DH,DEM,THX,dx,SMBtype,oname,gifout,P1,P2,P3,P4,H2010); 


        end
end
toc
% disp(sum(VALID(1:ig)))
