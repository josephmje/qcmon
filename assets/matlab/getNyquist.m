function [roimsk,roimsk2,PEave,ROave]=getNyquist(Y,mask)

% added signal mask that gets passed in to remove any pixels from phantom
% so that can safely move away from edges and avoid phantom signal
% if phantom is postioned strangely some day


[Nx,Ny]=size(Y)
SE=strel('square',3);
Smask=imdilate(mask,SE); % need to dilate it here to make not to catch "halo" 
[rM,cM]=find(Smask);
Frame=Smask*0;
Frame(find(Y==0))=1;
Frame=imdilate(Frame,SE);

mincolSmask=min(cM);
maxcolSmask=max(cM);
minrowSmask=min(rM);
maxrowSmask=max(rM);

roimsk=Smask*0;
roimsk(1:minrowSmask-5,:)=1;
roimsk(maxrowSmask+5:Ny,:)=1;
%clean up for areas to left and right of phantom (where noise in FE-direction is calculated)
roimsk(:,1:mincolSmask-3)=0;
roimsk(:,maxcolSmask+3:Ny)=0;

roimsk2=Smask*0;
roimsk2(:,1:mincolSmask-3)=2;
roimsk2(:,maxcolSmask+3:Ny)=2;

%remove edge effects 
roimsk(find(Frame))=0;
roimsk2(find(Frame))=0;

nd=Y(find(roimsk));

tmp=Y;

tmp(find(roimsk))=max(max(Y));
tmp(find(roimsk2))=max(max(Y))*1.2;

nd1=Y(find(roimsk));
nd2=Y(find(roimsk2));
% 

PEave=mean(nd1)
ROave=mean(nd2)

NyqRatio=PEave/ROave


% h5=figure(5)
% h5.Units='inches'
% h5.Position=[10 7 10 10];
% imagesc(tmp)
% axis image
% title(['NyquistRatio=',num2str(NyqRatio)])
% colormap(gray)