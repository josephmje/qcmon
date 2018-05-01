function [RatioB0,diax,diay]=getB0Distortion(DWI,nb0,PAR)

close all

[Nx, Ny, numimgs] = size(DWI);

DWIb0=DWI(:,:,1:nb0);
DWIb01=squeeze(DWI(:,:,1));
DWIb0ave=mean(DWIb0,3);

figure
subplot(1,2,1)
imagesc(DWIb0ave)
axis image
colormap(gray)
subplot(1,2,2)
imagesc(DWIb01)
axis image

maskb0=makemask(DWIb0ave,1,PAR);

[rmask,cmask]=find(maskb0);

rmasksort=sort(rmask);
lensm=length(rmask)
rmaskminave=mean(rmasksort(1:10))
rmaskmaxave=mean(rmasksort(lensm-9:lensm))
diay=rmaskmaxave-rmaskminave

cmasksort=sort(cmask);
cmaskminave=mean(cmasksort(1:10))
cmaskmaxave=mean(cmasksort(lensm-9:lensm))
diax=cmaskmaxave-cmaskminave

RatioB0=diay/diax

e=edge(maskb0,'canny');
showb0dist=squeeze(DWIb01);
maxv=max(max(showb0dist))
for rr=1:Nx
    for cc=1:Ny
        if(e(rr,cc))
            showb0dist(rr,cc)=maxv*1.1;
        end
    end
end
showb0dist(round(rmaskminave),:)=maxv;
showb0dist(round(rmaskmaxave),:)=maxv;
showb0dist(:,round(cmaskminave))=maxv;
showb0dist(:,round(cmaskmaxave))=maxv;

h1=figure(1)
set(h1, 'Visible', 'off');
set(h1, 'Units', 'inches');
set(h1, 'Position', [10 7 10 10]);

imagesc(showb0dist)
axis image
strtitle=strcat('RatioB0=',num2str(RatioB0),'/diax=',num2str(diax),'/diay=',num2str(diay))
title(strtitle)

% fig1name=strcat(DIRfig,'B0Distortion-',PAR)
% print('-f1',fig1name,'-djpeg')
