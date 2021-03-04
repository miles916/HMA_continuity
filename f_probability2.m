function [umult,SIG_umult,f] = f_probability2(THXdist)

% THXdist=cur.F.thx(cur.mask);


%% set up param distributions

N=10000;
rs1=rand(N,3);
% 
% u_b = rs1(:,1);
% thx=50+400.*rs1(:,2);
% n=2.8+0.4.*rs1(:,3)

rs2=randn(N,3);

rs3=lognrnd(-2,0.75,N,1);
% u_b = 0.5+(0.5/3).*rs2(:,1);
% u_b = abs((0.8/3).*rs2(:,1));
% u_b = abs((1/1.5).*(rs3(:,1)-0.75));
u_b = (rs1(:,1));
u_b(u_b>1)=1; u_b(u_b<0)=0;
% thx=250+(200/3).*rs2(:,2);
thx=THXdist(ceil(length(THXdist).*rs1(:,1))); %index values randomly from the distribution
n=3+(0.2/3).*rs2(:,3);

V.S=20;

%%
tic
for ip=1:N
        cTHX=thx(ip);
        if cTHX>1
            for z=1:(cTHX)
               cSt(z) = u_b(ip)*V.S+V.S*(1-u_b(ip))*(1-((cTHX-z)/cTHX)^(n(ip)+1));
            end
            V.Sm(ip)=nanmean(cSt);
        else
            V.Sm(ip)=NaN;
        end
%         ip
end
toc
f=V.Sm./V.S;
%%
% histogram(f)
umult=mean(f);
SIG_umult=std(f);



