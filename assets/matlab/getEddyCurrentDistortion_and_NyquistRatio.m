function [avevoxsh,errvoxsh,NyqRatio]=getEddyCurrentDistortion_and_NyquistRatio(DWI,nb0,PAR)

close all

[Nx, Ny, numimgs] = size(DWI);
numr=ceil(sqrt(numimgs))
numc=numr;

mask=makemask(DWI,nb0,PAR); % gives figs 1 & 2

SigM=mask;
SigM1=SigM(:,:,1);

for i=1:numimgs
    DY(:,:,i)=abs(SigM(:,:,i)-SigM1);
    clear tmp
    tmp(:,:)=DY(:,:,i);
    tmp2=tmp;
    tmp2(find(tmp==2))=0;
    DY(:,:,i)=tmp2(:,:);
end

clear tmp*

h3=figure(3)
set(h3, 'Visible', 'off');
set(h3, 'Units', 'inches');
set(h3, 'Position', [10 7 20 10]);


for i=2:numimgs
    subplot(numr,numc,i)
    imagesc(DY(:,:,i))
    set(gca, 'fontsize', 6)
    title(['Image #',num2str(i)])
    colormap(gray)
    axis image
    axis off
end


[rmask,cmask]=find(SigM1);

cmaskmax=max(cmask)
cmaskmin=min(cmask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% go thru each col and find #adj mask pixels for thickness along y up and down
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xctr=floor(Nx/2)

for i=2:numimgs

    tmp(:,:)=DY(:,:,i);

    %% add this to remove left&right edge effects (SCedit-June2016)
    tmp(:,cmaskmin-1:cmaskmin+1)=0;
    tmp(:,cmaskmax-1:cmaskmax-1)=0;

    %new method (June2016)
    tmpup=tmp;
    tmpup(floor(xctr)+1:Nx,:)=0;
    tmpdown=tmp;
    tmpdown(1:floor(xctr),:)=0;

%     figure(10+i)
%     subplot(2,2,1)
%     imagesc(tmpup)
%     axis image
%     subplot(2,2,2)
%     imagesc(tmpdown)
%     axis image

    numcol2(i)=(cmaskmax-2)-(cmaskmin+1)+1;
    pixshcolup{i}=sum(tmpup,1);
    pixshcoldown{i}=sum(tmpdown,1);
    avecolpixsh2{i}=mean([pixshcolup{i};pixshcoldown{i}],1);

    avepixsh2(i)=sum(avecolpixsh2{i},2)/numcol2(i);

end

h4=figure(4)
set(h4,'Visible', 'off');
set(h4, 'Units', 'inches');
set(h4, 'Position', [10 7 15 10]);

avevoxsh=mean(avepixsh2(nb0+1:numimgs))
errvoxsh=mean(avepixsh2(2:nb0))
plot(avepixsh2,'b*-')
title(['Pixel Shifts for Phantom: avevoxshift(err)=',num2str(avevoxsh),'(',num2str(errvoxsh),')'])
%
%
% fig3name=strcat(DIRfig,'DiffMasks-',PAR)
% print('-f3',fig3name,'-djpeg')
%
% fig4name=strcat(DIRfig,'Plot-EddyCurrentDist-',PAR)
% print('-f4',fig4name,'-djpeg')

%% get Nyquist Ratio

%get Nyquist Ratio for each b=0, then average
h5=figure(5)
set(h5, 'Visible', 'off');
set(h5, 'Units', 'inches');
set(h5, 'Position', [10 7 20 5]);

for i=1:nb0
    tmp=DWI(:,:,i);
    [ROImaskPE,ROImaskRO,PEaveALL(i),ROaveALL(i)]=getNyquist(tmp,mask(:,:,i));

    NyqRatioALL(i)=PEaveALL(i)/ROaveALL(i);

    tmpfig=tmp;
    tmpfig(find(ROImaskPE))=max(max(tmp));
    tmpfig(find(ROImaskRO))=max(max(tmp))*1.2;

    subplot(1,nb0,i)
    imagesc(tmpfig)
    axis image
    title(['nDWI #',num2str(i),'// NyquistRatio=',num2str(NyqRatioALL(i))])
    colormap(gray)
end

NyqRatio=mean(NyqRatioALL)

