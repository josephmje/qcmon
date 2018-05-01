function mask=makemask(DWI,nb0,PAR)

close all
disp('in makemask..')

[Nx, Ny, numimgs] = size(DWI);
N=Nx;
mask=DWI*0;
lenfillb0(1:nb0)=0;
lenfillchk=Nx*Ny/2;

numr=ceil(sqrt(numimgs))
numc=numr;

chknumimgs=numr*numc
numimgs

for i=1:numimgs

    DWItmp=DWI(:,:,i);
    DWItmpsm=medfilt2(DWItmp,[11 11]);
    for ll=1:2
    DWItmpsm2=medfilt2(DWItmpsm,[11 11]);
    DWItmpsm=DWItmpsm2;
    end

    e=edge(DWItmpsm,'canny'); %finds lots of edges
%
    eorig=e;

    e(1:15,:)=0;
    e(:,1:15)=0;
    e(N-15:N,:)=0;
    e(:,N-15:N)=0;
    radc=30;
    circ=makecirc(N,N/2,N/2,radc);
    e(find(circ))=0; % remove some central edges from 'canny'

    ef2=imfill(e,[N/2,N/2]); % flood from central points (to overcome closed edges within phantom)
    lenfill=length(find(ef2));
    if i<nb0+1
        lenfillb0(i)=lenfill;
        AVElenfill=mean(lenfillb0(find(lenfillb0)));

    end

    se=strel('disk',3,0);
    if lenfill>lenfillchk


        'dilating then filling...'
        e2=imdilate(e,se);
        locations=find(circ);
        ef2=imfill(e2,locations);
        ef2=imfill(ef2,'holes');
        ef2=imerode(ef2,se); % erode the first imdilate around phantom
        lenfill=length(find(ef2))
%

        lenfillb0(i)=lenfill;
        AVElenfill=mean(lenfillb0(find(lenfillb0)));
%
        if lenfill>lenfillchk
            'still not closing borders...'

            e2=imdilate(e,se);% try to dilate x2
            e2=imdilate(e2,se);
            ef2=imfill(e2,locations);
            ef2=imerode(ef2,se); % erode x2
            ef2=imerode(ef2,se);
            lenfill=length(find(ef2))
            if i==1
                lenfillb0(i)=lenfill;
                AVElenfill=lenfill;
            end

        end
    end

    minsize=AVElenfill*0.92;
    if lenfill<minsize
            'not full'

            clear locations circ
            radc=radc+2;
            circ=makecirc(N,N/2,N/2,radc);
            locations=find(circ);
            ef3=imfill(ef2,locations);
            lenfill=length(find(ef3))
            lenfill=length(find(ef3))

        if lenfill>lenfillchk
               'got a problem!'
               'borders still not closing ...'

                ef3d=imdilate(ef2,se);% try to dilate x2
                ef3d=imdilate(ef3d,se);
                ef3f=imfill(ef3d,locations);
                ef3f=imerode(ef3f,se); % erode x2
                ef3f=imerode(ef3f,se);
                lenfill=length(find(ef3f))
                ef3=ef3f;
        end

        while lenfill<minsize
            'still not full'



            clear locations circ
            radc=radc+2
            circ=makecirc(N,N/2,N/2,radc);
            locations=find(circ);
            ef3f=imfill(ef3,locations);
            lenfill=length(find(ef3f))
            if lenfill>lenfillchk
               'got a problem!'
                 'borders still not closing ...'

                ef3d=imdilate(ef3,se);% try to dilate x2
                ef3d=imdilate(ef3d,se);
                ef3f=imfill(ef3d,locations);
                ef3f=imerode(ef3f,se); % erode x2
                ef3f=imerode(ef3f,se);
                lenfill=length(find(ef3f))

            end

            ef3=ef3f;
        end

        ef2=ef3;

    end

    ef=ef2;
    ef1=imerode(ef,se); % erode outer edges
    ef2=imdilate(ef1,se); % get back the phantom edge
    ef2=imfill(ef2,'holes'); % fill any leftover holes in phantom
    lenfill=length(find(ef2));
    if lenfill>lenfillchk
        'edge NEVER closed!'
         pause
    end


    if i<nb0+1, lenfillb0(i)=lenfill; AVElenfill=mean(lenfill(find(lenfill))); end

    mask(:,:,i)=ef2(:,:);
    %pause(0.1)
end

if numr>2

    h1=figure(1)
    set(h1, 'Visible', 'off');
    set(h1, 'Units', 'inches');
    set(h1, 'Position', [10 7 20 10]);
    
    for i=1:numimgs
	i
        subplot(numr,numc,i)
        imagesc(DWI(:,:,i))
        if i<nb0+1
            title('nDWI')
        else
            title('DWI')
        end
        axis image
    end

    h2=figure(2)
    set(h2, 'Visible', 'off');
    set(h2, 'Units', 'inches')
    set(h2, 'Position', [10 7 20 10]);
    
    for i=1:numimgs
        subplot(numr,numc,i)
        imagesc(mask(:,:,i))
        if i<nb0+1
            title('nDWI')
        else
            title('DWI')
        end
        axis image
    end

%     fig1name=strcat(DIRfig,'CentralSlice-',PAR)
%     print('-f1',fig1name,'-djpeg')
%
%     fig2name=strcat(DIRfig,'MaskCentralSlice-',PAR)
%     print('-f2',fig2name,'-djpeg')

end

