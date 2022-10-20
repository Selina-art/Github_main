%% ASSIGNMENT2 - Fourier Transformation and Filtering for AESM1511, WS22, TU Delft
% Two-dimensional Fourier transformation from the time-space domain to the
% frequency-wavenumber domain.
%
% Other m-files required: none
% MAT-files required: none
%
% Author: Selina Jäckle, Louise
% Date: 07.10.22

%% SETTINGS
clear all;  % clear all previous variables
close all;  % close all figures
clc;        % clear commandline

%% PARAMETERS
filename = 'refl_3layers_fp50_dx0p5_500_rvz.bin'; % file containing the reflection data.
traces = 401; % number of traces
time_samples = 1001; % number of time samples
time_sampling = 1; % time sampling in ms
frequency_range = 1:200; % new frequency range after Question 2.
filter1_frequencies = 1:60; % frequencies set to zero by filter 1 (Question 4)
filter_offset = 1:175; % uppper limit visually determined from reflection data both in time and frequency domain.
lower_frequencies = 1:30; % frequencies set to zero by filter 2 and 3
filter_2_and_3_frequencies = 30:60; % frequencies filtered by filter 2 and 3
frequency_sampling = 1; % frequency sampling in Hz
frequency_number = length(frequency_range); % number of frequencies
space_sampling = 2.5; % space sampling in m

%% Question 1 - Read the file
fid = fopen(filename, 'r'); % open the file specified in parameters
temp = fread(fid, traces*time_samples, 'float32', 0, 'ieee-le'); % read the file
fclose(fid); % close the file
data_refl(:,:) = reshape(temp,time_samples,traces); % assign the data from the file to a new variable
% Plot the reflection data
figure
imagesc(data_refl) % image of the reflection data
colormap("gray") % makes colormap gray
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title(gca, 'Common-source gather')
xlabel('Distance (m)')
ylabel('Two-way travel time (s)')

%% Question 2
frequency_data = fft(data_refl,[],1)*time_sampling; % Fast Fourier Transform of the reflection data from time to frequency domain
disp('Q2: We cannot plot complex numbers on a 2D plot, so we take the absolute value.')
frequency_data_positive = frequency_data(frequency_range,:); % select the positive frequencies from 0 to 200 only.
% Plot positive frequencies only
figure
imagesc(abs(frequency_data_positive)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Reflection data in frequency-space domain')
xlabel('Distance (m)')
ylabel('Frequency (Hz)')

%% Question 3
disp('Surface waves? But unsure about its velocity (180m/s !?!)')

%% Question 4
disp('The low-cut filter is not applied to all horizontal distances, only until the amplitude drops.')
% Filter 1) Set frequencies lower than 60 Hz to zero.
frequency_filter1 = frequency_data_positive; % Copy the frequency data.
frequency_filter1(filter1_frequencies,filter_offset) = 0; % Set all frequencies below 60 to zero, up to a certain offset.

% Filter 2) Set frequencies from 60 Hz to 30 Hz to go linearly to zero
frequency_filter2 = frequency_data_positive; % Copy the frequency data.
frequency_filter2(lower_frequencies,filter_offset) = 0; % Set frequencies below 30 Hz to zero.
linear_filter = linspace(0,1,31); % Create a linear filter from 0 to 1, with 60-30+1 values.
frequency_filter2(filter_2_and_3_frequencies,filter_offset).*transpose(linear_filter); % Filter data to go lineary to zero between 60 and 30 Hz.

% Filter 3) Set frequencies from 60Hz to 30Hz to go to zero following a
% sinusoidal curve (quarter of a period)
frequency_filter3 = frequency_data_positive; % Copy the frequency data.
frequency_filter3(lower_frequencies,filter_offset) = 0; % Set frequencies below 30 Hz to zero.
sin_curve = sin(linspace(pi/2,pi,31)); % Create a sinusoidal filter (quarter of a period), from 1 to 0.
frequency_filter3(filter_2_and_3_frequencies,filter_offset).*transpose(sin_curve); % Filter data to go to zero between 60 and 30 Hz following a sinusoidal curve.

% Transform data back to the time-space domain. Three times, for each
% filter.
time_domain_filter1 = 2*real((frequency_number)*ifft(frequency_filter1,[],1)*(frequency_sampling));
time_domain_filter2 = 2*real((frequency_number)*ifft(frequency_filter2,[],1)*(frequency_sampling));
time_domain_filter3 = 2*real((frequency_number)*ifft(frequency_filter3,[],1)*(frequency_sampling));
%% Plot filtered data.
figure
t = tiledlayout(2,2);
nexttile
imagesc(abs(frequency_data_positive)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Reflection data in frequency-space domain')
nexttile
imagesc(abs(frequency_filter1)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Filter 1')
nexttile
imagesc(abs(frequency_filter2)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Filter 2: Linear filter')
nexttile
imagesc(abs(frequency_filter3)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Filter 3: Sinusoidal filter')
t.Padding = 'compact';
t.TileSpacing = 'compact';

%% Plot the filtered data in time domain.
figure
t = tiledlayout(2,2);
nexttile
imagesc(abs(data_refl)) % Image of the reflection data in the frequency domain
colormap("gray") % makes colormap gray
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title(gca, 'Common-source gather')
xlabel('Distance (m)')
ylabel('Two-way travel time (s)')
nexttile
imagesc(abs(time_domain_filter1)) % Image of the reflection data in the frequency domain
colormap("gray") % makes colormap gray
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Filter 1')
xlabel('Distance (m)')
ylabel('Two-way travel time (s)')
nexttile
imagesc(abs(time_domain_filter2)) % Image of the reflection data in the frequency domain
colormap("gray") % makes colormap gray
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Filter 2: Linear filter')
xlabel('Distance (m)')
ylabel('Two-way travel time (s)')
nexttile
imagesc(abs(time_domain_filter3)) % Image of the reflection data in the frequency domain
colormap("gray") % makes colormap gray
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Filter 3: Sinusoidal filter')
xlabel('Distance (m)')
ylabel('Two-way travel time (s)')
t.Padding = 'compact';
t.TileSpacing = 'compact';

disp('Q4: Our filters do not give very different results, and most events seem to be preserved. Perhaps our filters our not working properly.')
%% Question 5
frequency_wavenumber_data = fftshift(fft(frequency_data_positive,[],1)*space_sampling,1);
%
figure
imagesc(abs(frequency_wavenumber_data)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Question 5')
ylabel('Frequency (Hz)')
xlabel('Wavenumber (1/m)')

disp('Q5: The energy is concentrated in a single trace similar to the one in time-sapce and frequency-space domains.')
%% Question 6

frequency_wavenumber_data_2ndtrace = fftshift(fft(frequency_data_positive(:,1:2:end),[],2)*space_sampling,2);
frequency_wavenumber_data_4thtrace = fftshift(fft(frequency_data_positive(:,1:4:end),[],2)*space_sampling,2);
frequency_wavenumber_data_8thtrace = fftshift(fft(frequency_data_positive(:,1:8:end),[],2)*space_sampling,2);

%%
figure
subplot(2,2,1)
imagesc(abs(frequency_wavenumber_data)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('Frequency wavenumber data')
ylabel('Frequency (Hz)')
xlabel('Wavenumber (1/m)')
subplot(2,2,2)
imagesc(abs(frequency_wavenumber_data_2ndtrace)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('2nd trace')
ylabel('Frequency (Hz)')
xlabel('Wavenumber (1/m)')
subplot(2,2,3)
imagesc(abs(frequency_wavenumber_data_4thtrace)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('4th trace')
ylabel('Frequency (Hz)')
xlabel('Wavenumber (1/m)')
subplot(2,2,4)
imagesc(abs(frequency_wavenumber_data_8thtrace)) % Image of the reflection data in the frequency domain
colormap("default") % makes colormap full color
colorbar % colorbar visible
set(gca,'linewidth',2, 'fontsize', 14) % changes linewidth to 2, fontsize to 14 for the entire figure.
title('8th trace')
ylabel('Frequency (Hz)')
xlabel('Wavenumber (1/m)')

disp('Q6: The main energy is split as the trace number increases.')
disp('Q7: See attachement.')