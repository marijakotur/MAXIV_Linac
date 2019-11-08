% bunch length vs charge
close all
clear all

charge = [12 20:20:180];
% charge = [6 9 15:10:45 65:10:95 115:10:175];

no_images = 5;
L01phase = 5;

%dispersion
dispersion = 0.53; %in m
C = dispersion * sind(L01phase)*2*pi*300e6*cosd(L01phase);
    
pix_per_m = 331.5/5*1e3;
xaxis = (1:1280)/pix_per_m;
% xaxis = (1:1280);
xaxis = xaxis*(-1)*cotd(L01phase)/(dispersion*2*pi*3e9)*1e12;
xaxis = xaxis - xaxis(end);


figure
hold on
cm = jet(length(charge));
for i=1:length(charge)
    fnamestr = [num2str(charge(i)) 'pC_img*.tiff'];
    fnames = dir(fnamestr);
    
    fnamestr_bkgnd = [num2str(charge(i)) 'pC_bkgnd_*.tiff'];
    fnames_bkgnd = dir(fnamestr_bkgnd);
    
    lineout = zeros(1,1280);
    lineout_bkgnd = zeros(1,1280);
    
    str{i} = [num2str(charge(i)) ' pC'];
    
        subplot(2,5,i)
        hold on
        title([num2str(charge(i)) 'pC'])
    
    for j=1:length(fnames)
        fname = fnames(j).name;
        data = imread(fname);
        lineout = lineout/j + sum(data(455:480,:),1);
        
        fname_bkgnd = fnames_bkgnd(j).name;
        bkgnd = imread(fname_bkgnd);
        lineout_bkgnd = lineout_bkgnd/j + sum(bkgnd(455:600,:),1);
        %         plot(sum(data(455:485,:),1))
    end
    f = sum(lineout(end-10:end))/sum(lineout_bkgnd(end-10:end));
    lineout = lineout/f-lineout_bkgnd;
    lineout_smoothed = smoothbox(9,lineout/max(lineout));
    fwhm_(i) = abs(fwhm(xaxis,lineout_smoothed));
    
    plot(xaxis,lineout_smoothed,'Color',cm(i,:))
%         plot(lineout_bkgnd,'r')
end
axis tight
box on
legend(str,'location','best')
grid on
ylabel('Signal [arb.]')
xlabel('Time [ps]')

figure
plot(charge,fwhm_,'o-')
xlabel('Charge [pC]')
ylabel('Bunch duration [ps]')


