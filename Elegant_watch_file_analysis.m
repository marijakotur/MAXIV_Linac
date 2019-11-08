% read Elegant watch file
clear all
close all


m0 = 0.510998910; % MeV
E_min = 50; % MeV, minimum energy to be analyzed
q0 = 200e-12; % charge in C

% fname = 'max4-I_BC1_DIA_SCRN_1.w1.txt';
% fname = 'max4-I_MS1_DIA_SCRN_1.w1';
% fname = 'w-init.sdds';
fname1 = 'w-lh1.sdds';
fname2 = 'w-lh2.sdds';
fname3 = 'w-lh3.sdds';
fname4 = 'w-lh4.sdds';

% fname = 'w-achr1a.sdds';

e = 1.60217662e-19;
c = 299792458;

watch_files = dir('w-lh*.sdds');
for j = 1:length(watch_files)
    
    watch_files(j).name
    command = ['sddsconvert ' watch_files(j).name ' ' watch_files(j).name '.txt -ascii'];
    system(command);
    
    dataFile = [watch_files(j).name];
    
    M = dlmread([dataFile '.txt'], '', 47, 0); % Might need to open this file IN MATLAB to check where matrix starts (row 48 atm)
    M = M(M(:,6)>E_min/m0,:); % Remove any energies that are <E_min MeV
    M_sort_by_t = sortrows(M,5);     % 6-D particle data sorted by time
    M_sort_by_p = sortrows(M,6); % We also sort by energy for transverse phase space plots
    p_sort_by_p = M_sort_by_p(:,6)';  % Dimensionless momentum p = gamma*beta
    nbrOfParticles = length(p_sort_by_p);
    
    t_sort_by_t = M_sort_by_t(:,5)';  % time
    t_sort_by_t = t_sort_by_t - mean(t_sort_by_t);
    p_sort_by_t = M_sort_by_t(:,6)';  % Dimensionless momentum p = gamma*beta
    x_sort_by_t = M_sort_by_t(:,1)';  % x-coordinate
    y_sort_by_t = M_sort_by_t(:,3)';  % y-coordinate
    xp_sort_by_t = M_sort_by_t(:,2)'; % x-prime
    yp_sort_by_t = M_sort_by_t(:,4)'; % y-prime
    
    figure(1)
    hold on
    subplot(2,2,j)
    plot(x_sort_by_t,y_sort_by_t,'.')
    axis equal
    axis([-3.5, 3.5,-3.5,3.5]*1e-4)

    energy_sort_by_t = p_sort_by_t*m0;
    
    
    if j==1
        p = polyfit(t_sort_by_t,energy_sort_by_t,3);
        
        E_smooth = p(1)*t_sort_by_t.^3 + p(2)*t_sort_by_t.^2 + p(3)*t_sort_by_t + p(4);
        
    end
    
%     energy_sort_by_t = energy_sort_by_t - 0*E_smooth;
    
    figure(2)
    hold on
    subplot(2,2,j)
    plot(t_sort_by_t,energy_sort_by_t,'.')
    box on
    xlabel('time [s]')
    ylabel('Energy (MeV)')
%     axis([-1e-14,1e-14,285,285.5])
    axis([-2e-14,2e-14,274.5,275.5])
%     axis tight

    % % Calculate the bunch and slice parameters
    % nBins = 500;
    % beamParameters = bunchParameters(M_sort_by_t, q0, nBins);
    
    % convert t to phase of the laser field
    lambda_laser = 790*1e-9;
    T_laser = lambda_laser/c;
    phi = 2*pi*t_sort_by_t/T_laser;
    phiWrapped = wrapTo2Pi(phi);
    
%     figure(3)
%     subplot(length(watch_files),1,j)
%     hold on
%     plot(phiWrapped,energy_sort_by_t,'.')
%     box on
%     xlabel('Phase (rad)')
%     ylabel('Energy [MeV]')
%     axis tight
    
end


%analyze sig file
% sigFile = 'max4.sig';
sigFileName = 'max4.sig.txt';
fid = fopen(sigFileName);
dataSubLine = fgetl(fid);
k = 1;
while ~isequal(dataSubLine, -1) % Find beginning of data
    dataSubLine = fgetl(fid);
    if strfind(dataSubLine, '_BEG_') % Breaks loop if _BEG_ is found
        break;
    end
    k = k + 1;
end
% fclose(sigFileName);
k = k + 1;
SIG = importdata(sigFileName, ' ', k); % File has k headerlines
s_coord = str2num(cell2mat(SIG.textdata(k+1:end,1))); % First column has coordinates
elementType = SIG.textdata(k+1:end, 4); % 4th column has element type info
s_x = SIG.data(:,50)*1e6; % Sigma x, um
s_xp = SIG.data(:,51); % Sigma x', rad
s_y = SIG.data(:,52)*1e6; % Sigma y, um
s_yp = SIG.data(:,53); % Sigma y', rad
s_dp = SIG.data(:,55); % Sigma dp/p, rms
s_t = SIG.data(:,56)*10^15; % Sigma t, fs
s_e_x = SIG.data(:,57); % Geometric emittance in x, m rad
s_e_nx = SIG.data(:,58)*1e6; % Normalized emittance in x, mm mrad
% s_e_x_noDisp = SIG.data(:,59); % Geometric emittance in x, m rad, no dispersive contribution
% s_e_nx_noDisp = SIG.data(:,60)*1e6; % Normalized emittance in x, mm mrad, no dispersive contribution
s_e_y = SIG.data(:,61); % Geometric emittance in y, m rad
s_e_ny = SIG.data(:,62)*1e6; % Normalized emittance in y, mm mrad
% s_e_y_noDisp = SIG.data(:,63); % Geometric emittance in y, m rad, no dispersive contribution
% s_e_ny_noDisp = SIG.data(:,64)*1e6; % Normalized emittance in y, mm mrad, no dispersive contribution
s_beta_x = SIG.data(:,65); % Statistical Twiss parameters for the whole beam, no dispersive contribution
s_alpha_x = SIG.data(:,66);
s_beta_y = SIG.data(:,67);
s_alpha_y = SIG.data(:,68);


figure(4)
hold on
box on
plot(s_coord,s_beta_x,'b-')
plot(s_coord,s_beta_y,'r-')
legend('\beta_x','\beta_y')
xlabel('s [m]')
