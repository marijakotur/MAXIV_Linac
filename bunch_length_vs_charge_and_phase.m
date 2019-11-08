% bunch length vs charge
close all
clear all

charge = [23];
L01phase = 5:5:25;

no_images = 3;

%dispersion
dispersion = 0.53; %in m
C = dispersion * sind(L01phase)*2*pi*300e6.*cosd(L01phase);

pix_per_m = 331.5/5*1e3;
% xaxis = (1:1280)/pix_per_m;
% xaxis = (1:1280);

figure
hold on
cm = jet(length(L01phase));

for i=1:length(L01phase)
    fnamestr = [num2str(charge) 'pC_' num2str(L01phase(i)) 'deg_img*.tiff'];
    fnames = dir(fnamestr);
    
    lineout = zeros(1,1280);
    xaxis = (1:1280)/pix_per_m*(-1)*cotd(L01phase(i))/(dispersion*2*pi*3e9)*1e12;
    xaxis = xaxis - xaxis(end);
    
    str{i} = [num2str(L01phase(i)) ' pC'];
    
    subplot(2,5,i)
    hold on
    
    for j=1:length(fnames)
        fname = fnames(j).name;
        data = imread(fname);
        lineout = lineout/j + sum(data(440:500,:),1);
        
        %         plot(sum(data(455:485,:),1))
    end
    f = mean(lineout(end-10:end));
    lineout = lineout-f;
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
plot(L01phase,fwhm_,'o-')
xlabel('Charge [pC]')
ylabel('Bunch duration [ps]')


