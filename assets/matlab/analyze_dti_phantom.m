%  analyze_dti_phantom(dwi, fa, bval, output, nyqopt)
%
%  'dwi':    4D diffusion weighted image
%  'fa':     FA map from DTIfit
%  'bval':   B value files from dcm2nii
%  'output': full path to output prefix
%  'accel': ('y', 'n') 'n' to measure nyquist ghost on non-accelerated data.

function analyze_dti_phantomSC(dwi, fa, bval, output_prefix, accel)
try
    %% Part 1: load data
    disp('loading data');
    bval = dlmread(bval);
    dwi = load_nifti(dwi);
    fa = load_nifti(fa);

    if accel=='y'
        PAR='PAR';
    elseif accel=='n'
        PAR='NPAR';
    else
        disp('problem with accel value')
        exit(1)
    end

    % calculate number of directions / b0 volumes
    ndir = length(find(bval>0));
    nb0 = length(find(bval==0));
    bvalueall=bval(find(bval>0));
    bvalue=bvalueall(1);

    % initalize outputs
    disp('opening output files');
    outname01 = strcat(output_prefix, 'Section2.3.1_SNR_ADC.csv');
    outname02 = strcat(output_prefix, 'Section2.3.2_B0DistortionRatio.csv');
    outname03 = strcat(output_prefix, 'Section2.3.3_EddyCurrentDistortions.csv');
    outname04 = strcat(output_prefix, 'Section2.3.4_AveNyqRatio.csv');
    outname05 = strcat(output_prefix, 'Section2.3.5_FAvalues.csv');

    fid01 = fopen(outname01, 'w');
    fid02 = fopen(outname02, 'w');
    fid03 = fopen(outname03, 'w');
    fid04 = fopen(outname04, 'w');
    fid05 = fopen(outname05, 'w');

    % print headers
    fprintf(fid01, '%s,%s,%s,%s,%s\n', ...
                   'AVE(SNR)0', '%CV(SNR)0', 'AVE(SNR)DWI','%CV(SNR)DWI', 'ADC');
    fprintf(fid02, '%s\n', ...
                   'B0 Distortion Ratio');
    fprintf(fid03, '%s,%s\n', ...
                   'avevoxshift','%errorvoxshift');
    fprintf(fid04, '%s\n', ...
                   'Nyquist Ratio');
    fprintf(fid05, '%s,%s\n', ...
                   'AVE(FA)','STD(FA)');

    % load DWI, averaging over 3 central slices. dim1,2=AXIAL, dim3=all directions)
    dims = size(dwi.vol);
    central_slice = ceil(dims(3)/2);
    DWI = mean(dwi.vol(:,:, central_slice-1:central_slice+1, :), 3);
    
    DWIrot=imrotate(DWI,90); % get PE along vertical direction
    clear DWI
    DWI=DWIrot;
    
    [Nx, Ny, numimgs] = size(DWI);

    % load FA, taking central slice only
    FA = fa.vol(:,:,central_slice);

    clear fa dwi

    %% 2.3.1 SNR Measurements
    [AVEsnr0, CVsnr0, AVEsnrDWI, CVsnrDWI] = getSNR(DWI, nb0, PAR);

    R = AVEsnrDWI/AVEsnr0;
    ADC = (log(R)/bvalue)*1000; % in units of 10^-3 s/mm^2

    DIRfig=output_prefix;
    fig1name=strcat(DIRfig,'DiffImgs-',PAR);
    print('-f1',fig1name,'-djpeg');

    fig2name=strcat(DIRfig,'StdPlotsHist-',PAR);
    print('-f2',fig2name,'-djpeg');

    fig3name=strcat(DIRfig,'SNRImgs-',PAR);
    print('-f3',fig3name,'-djpeg');

    fig4name=strcat(DIRfig,'SNRplots-',PAR);
    print('-f4',fig4name,'-djpeg');

    disp('done 2.3.1')
    pause(0.5)

    %% 2.3.2 B0 inhomogeneity
    [RatioB0, diax, diay] = getB0Distortion(DWI, nb0, PAR);

    fig1name=strcat(DIRfig,'B0Distortion-',PAR)
    print('-f1',fig1name,'-djpeg')

    disp('done 2.3.2')
    pause(0.5)

    %% 2.3.3 Eddy Current Distortions & 2.3.4 Nyquist Ratio
    [avevoxsh, errvoxsh, NyqRatio]=getEddyCurrentDistortion_and_NyquistRatio(DWI, nb0, PAR);

    fig1name = strcat(DIRfig, 'CentralSlice-', PAR);
    print('-f1', fig1name, '-djpeg')

    fig2name = strcat(DIRfig, 'MaskCentralSlice-', PAR);
    print('-f2', fig2name, '-djpeg')

    fig3name = strcat(DIRfig, 'DiffMasks-', PAR);
    print('-f3', fig3name, '-djpeg')

    fig4name = strcat(DIRfig, 'Plot-EddyCurrentDist-', PAR);
    print('-f4', fig4name, '-djpeg')

    fig5name = strcat(DIRfig, 'NyquistRatio-', PAR);
    print('-f5', fig5name, '-djpeg')

    disp('done 2.3.3 and 2.3.4')
    pause(0.5)

    %% 2.3.5 FA
    close all

    fav = FA(find(FA));
    aveFA = mean(fav);
    stdFA = std(fav);

    h1 = figure(1)
    set(h1, 'Visible', 'off');
    set(h1, 'Units', 'inches');
    h1.Position=[10 7 15 10];

    subplot(2,2,1)
    imagesc(FA, [0 0.1]);
    set(gca,'DataAspectRatio',[1 1 1]); colorbar;
    title('FA map');

    subplot(2,2,2)
    plot(fav);
    title('FA values');

    subplot(2,2,3)
    [n,x] = hist(fav,30); plot(x, n, 'k-');
    title(['FA mean(std)=', num2str(aveFA, '%5.3f'), '(', num2str(stdFA, '%5.3f'), ')']);

    %% Print cumulative values to files
    fprintf(fid01, '%5.2f,%5.2f,%5.2f,%5.2f,%5.2f\n', ...
                   AVEsnr0, CVsnr0, AVEsnrDWI, CVsnrDWI, ADC)
    fprintf(fid02, '%5.3f\n', ...
                   RatioB0)
    fprintf(fid03, '%7.4f,%7.4f\n', ...
                   avevoxsh, errvoxsh)
    fprintf(fid04, '%5.3f\n', ...
                   NyqRatio);
    fprintf(fid05, '%7.4f,%7.4f\n', ...
                   aveFA, stdFA);
    close all
    exit
catch
    exit(1)
end
end
