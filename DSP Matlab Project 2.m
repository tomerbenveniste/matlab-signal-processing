clc;
clear all;
close all;
 
% -------------------------- Question 1A ----------------------------------
% defining parameters
theta1 = pi/11.25;
theta2 = pi/3;
N = 30; % length of the discrete signals
n = 0:N-1; % discrete time vector
% defining the signals s[n],v[n],x[n]
s_n = 3*cos(theta1*n);
v_n = 5*sin(theta2*n);
x_n = s_n + v_n;

% showing the signals s[n],v[n],x[n] on one plot - Discrete Time Domain
figure; % new figure window
plot(n,s_n,'b-o','LineWidth',2);
hold on;
plot(n,v_n,'k-o','LineWidth',2);  
hold on;
plot(n,x_n,'r-o','LineWidth',2);   
title('The Signals s[n],v[n],x[n] in Time Domain','FontSize',18); 
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
legend([{'s[n]'};{'v[n]'};{'x[n]'}],'FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% Calculating the DFT's S[k],V[k],X[k]
Sd_k = fft(s_n); 
Vd_k = fft(v_n); 
Xd_k = fft(x_n); 

k = 0:N-1; % defining the frequncy vector

% showing the signals S[k],V[k],X[k] on one plot - Frequency Domain
figure; % new figure window
stem(k,abs(Sd_k),'b-o','LineWidth',2);
hold on;
stem(k, abs(Vd_k),'k-o','LineWidth',2);
hold on;
stem(k, abs(Xd_k),'r--o','LineWidth',2);
xlabel('k [Frequency Index]','FontSize',16);
ylabel('Magnitude','FontSize',16);
title('DFT Magnitude','FontSize',18);
legend([{'S^d[k]'};{'V^d[k]'};{'X^d[k]'}],'FontSize',16);
% defining ticks for w axis
xticks([pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi]); 
labels = {'pi/3', '2*pi/3', 'pi', '4*pi/3', '5*pi/3', '2*pi'};
xticklabels(labels);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% -------------------------- Question 1B ----------------------------------
% Zero Padding
num_zeros = 15; % number of zeros to pad
xz_n = [x_n, zeros(1,num_zeros)]; % creatind the zero padded signal
Nz = length(xz_n); % new length of the signal: Nz = 45
Xdz_k = fft(xz_n); % DFT of the zero padded signal

% creating new normalized axises for comparison
omega_N = (0 : N-1) * (2*pi / N); % for the original signal
omega_Nz = (0 : Nz-1) * (2*pi / Nz); % for the zero padded signal

% showing the X[k],Xz[k] on one plot - Frequency Domain
figure; % new figure window
stem(omega_N,abs(Xd_k),'b-o','LineWidth',2);
hold on;
stem(omega_Nz,abs(Xdz_k),'k-o','LineWidth',2);    
title('Comparison of Original DFT and Zero-Padded DFT','FontSize',18); 
xlabel('Normalized Frequency - \omega [rad/sec]','FontSize',16); 
ylabel('Magnitude [V]','FontSize',16);
legend([{'|X^d[k]| (N=30)'};{'|X_z^d[k]| (N=45)'}],'FontSize',16);
% defining ticks for w axis
xticks([pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi]); 
labels = {'$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$','$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20); % defining scale of axis size
xlim([0, 2*pi]); % limit x-axis to [0,2*pi]
grid on;

% -------------------------- Question 1C ----------------------------------
% additional samplings effect
add_samps = 15;
N2 = N + add_samps; % new length of the discrete signals: N=45
n2 = 0:N2-1; % discrete time vector
% defining the signals s2[n],v2[n],x2[n]
s_n2 = 3*cos(theta1*n2);
v_n2 = 5*sin(theta2*n2);
x_n2 = s_n2 + v_n2;

% Calculating the DFT's S2[k],V2[k],X2[k]
S2d_k = fft(s_n2); 
V2d_k = fft(v_n2); 
X2d_k = fft(x_n2); 

k2 = 0:N2-1; % defining the frequncy vector

% creating new normalized axis for comparison
omega_N2 = (0 : N2-1) * (2*pi / N2); % for the zero padded signal

% showing the X[k],X2[k] on one stem - Frequency Domain
figure; % new figure window
stem(omega_N,abs(Xd_k),'k-o','LineWidth',2);
hold on;
stem(omega_N2,abs(X2d_k),'r-o','LineWidth',2);    
title('Comparison of Original DFT and 15 added samplings DFT','FontSize',18); 
xlabel('\omega [rad/sec]','FontSize',16); 
ylabel('Magnitude [V]','FontSize',16);
legend([{'|X^d[k]| (N=30)'};{'|X_2^d[k]| (N=45)'}],'FontSize',16);
% defining ticks for w axis
xticks([pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi]); 
labels = {'$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$','$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20); % defining scale of axis size
xlim([0, 2*pi]); % limit x-axis to [0,2*pi]
grid on;

% -------------------------- Question 1D ----------------------------------
% calculating Parseval equation for both original and zero padded signals

% calculating Parseval equation for the original signal (N=30)
energy_time_orig = x_n * x_n'; % left side - time domain energy
energy_freq_orig = (Xd_k * Xd_k') / N; % right side - frequency domain energy

% calculating Parseval equation for the zero padded signal (Nz=45)
energy_time_zpad = xz_n * xz_n'; % left side - time domain energy
energy_freq_zpad = (Xdz_k * Xdz_k') / Nz; % right side - frequency domain energy

% -------------------------- Question 1E ----------------------------------
% defining the continous theta grid
theta_grid = linspace(0, 2*pi, 1000); 
% analytic DTFT of h1[n]
H1_theta = (1/3) * (1 + exp(-1j * theta_grid) + exp(-1j * 2 * theta_grid));
% defining magnitudes
H1_theta_mag = abs(H1_theta); % DTFT of H1 magnitude
Xd_k_mag = abs(Xd_k); % Xd[k] magnitude
norm_const = 1 / max(Xd_k_mag); % const for normalizing Xd[k]
Xd_k_mag_norm = norm_const * Xd_k_mag; % normalized Xd[k]

% showing the Xd[k],DTFT{H1} on one plot - Frequency Domain
figure; % new figure window
plot(theta_grid,abs(H1_theta),'b','LineWidth',2);
hold on;
stem(omega_N,abs(Xd_k_mag_norm),'r-o','LineWidth',2);    
title('Normalized Spectrum vs Filter Response','FontSize',18); 
xlabel('Normalized Frequency - \omega [rad/sec]','FontSize',16); 
ylabel('Normalized Magnitude (Max=1) [V]','FontSize',16);
legend([{'|H_1^f(\theta)| (Filter Response)'};{'|X_z^d[k]| (Normalized)'}],'FontSize',16);
% defining ticks for w axis
xticks([pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi]); 
labels = {'$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$','$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20); % defining scale of axis size
xlim([0, 2*pi]); % limit x-axis to [0,2*pi]
grid on;

% defining DFT{y1[n]}
h1_n = [1/3, 1/3, 1/3]; % h1 filter
Nh = length(h1_n); % filter length
N_lin_conv = N + Nh - 1; % DFT length for linear convolution

% defining DFT for zero padded x[n] and for zero padded h1[n] filter 
Xd_k_pad = fft(x_n, N_lin_conv); % DFT of zero padded x[n]
Hd_k_pad = fft(h1_n, N_lin_conv); % DFT of zero padded h1[n] filter
Y1d_k = Xd_k_pad .* Hd_k_pad; % DFT{y1[n]} after linear convolution

% calculating y1[n] with IFFT
y1_n = ifft(Y1d_k); 

% Normalized omega axis for Y_1^d[k] (Ny=32)
omega_Ny = (0 : N_lin_conv-1) * (2*pi / N_lin_conv);

% showing the X[k],Y1[k] on one stem - Frequency Domain
figure; % new figure window
stem(omega_N,abs(Xd_k),'k-o','LineWidth',2);
hold on;
stem(omega_Ny,abs(Y1d_k),'r-o','LineWidth',2);    
title('Final Comparison: Input vs Output DFT','FontSize',18); 
xlabel('\omega [rad/sec]','FontSize',16); 
ylabel('Magnitude [V]','FontSize',16);
legend([{'|X^d[k]| (N=30)'};{'|Y_1^d[k]| (N=32)'}],'FontSize',16);
% defining ticks for w axis
xticks([pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi]); 
labels = {'$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$','$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20); % defining scale of axis size
xlim([0, 2*pi]); % limit x-axis to [0,2*pi]
grid on;

% defining first 30 samples of the 32 samples of y1[n] 
y1_n_display = y1_n(1:N);

% showing the signals s[n],v[n],x[n],y1[n] on 4 subplots - Discrete Time Domain
figure('Name','Time Domain Decomposition'); % new figure window

% s[n] - plot - slow cosine
subplot(2,2,1); 
plot(n,s_n,'b-o','LineWidth',2);
title('1. Desired Signal s[n] (Low Freq)','FontSize',18);
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% v[n] - plot - fast sine
subplot(2,2,2);
plot(n,v_n,'k-o','LineWidth',2);  
title('2. Interference v[n] (High Freq)','FontSize',18);
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% x[n] - plot - input signal to filter
subplot(2,2,3);
plot(n,x_n,'r-o','LineWidth',2);   
title('3. Total Input x[n] = s[n] + v[n]','FontSize',18); 
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% y1[n] - plot - output signal from filter
subplot(2,2,4);
plot(n,y1_n_display,'m-o','LineWidth',2);   
title('4. Filtered Output y_1[n]','FontSize',18); 
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% -------------------------- Question 1F ----------------------------------
% analytic DTFT of h2[n]
H2_theta = (1 + exp(-1j * theta_grid));
norm_H2_theta = 0.5 * H2_theta;

% showing the Xd[k],DTFT{H2} on one plot - Frequency Domain
figure; % new figure window
plot(theta_grid,abs(norm_H2_theta),'b','LineWidth',2);
hold on;
stem(omega_N,abs(Xd_k_mag_norm),'r-o','LineWidth',2);    
title('Normalized Spectrum vs Filter Response','FontSize',18); 
xlabel('Normalized Frequency - \omega [rad/sec]','FontSize',16); 
ylabel('Normalized Magnitude (Max=1) [V]','FontSize',16);
legend([{'|H_2^f(\theta)| (Filter Response - Normalized)'};{'|X_z^d[k]| (Normalized)'}],'FontSize',16);
% defining ticks for w axis
xticks([pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi]); 
labels = {'$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$','$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20); % defining scale of axis size
xlim([0, 2*pi]); % limit x-axis to [0,2*pi]
grid on;

% defining DFT{y2[n]}
h2_n = [1, 1]; % h2 filter
Nh2 = length(h2_n); % filter length
N_lin_conv_2 = N + Nh2 - 1; % DFT length for linear convolution

% defining DFT for zero padded x[n] and for zero padded h2[n] filter 
Xd_k_pad_2 = fft(x_n, N_lin_conv_2); % DFT of zero padded x[n]
Hd_k_pad_2 = fft(h2_n, N_lin_conv_2); % DFT of zero padded h1[n] filter
Y2d_k = Xd_k_pad_2 .* Hd_k_pad_2; % DFT{y2[n]} after linear convolution

% calculating y1[n] with IFFT
y2_n = ifft(Y2d_k); 

% Normalized omega axis for Y_2^d[k] (Ny=31)
omega_Ny_2 = (0 : N_lin_conv_2 - 1) * (2*pi / N_lin_conv_2);

% showing the X[k],Y2[k] on one stem - Frequency Domain
figure; % new figure window
stem(omega_N,abs(Xd_k),'k-o','LineWidth',2);
hold on;
stem(omega_Ny_2,abs(Y2d_k),'r-o','LineWidth',2);    
title('Final Comparison: Input vs Output DFT','FontSize',18); 
xlabel('\omega [rad/sec]','FontSize',16); 
ylabel('Magnitude [V]','FontSize',16);
legend([{'|X^d[k]| (N=30)'};{'|Y_2^d[k]| (N=31)'}],'FontSize',16);
% defining ticks for w axis
xticks([pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi]); 
labels = {'$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$','$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20); % defining scale of axis size
xlim([0, 2*pi]); % limit x-axis to [0,2*pi]
grid on;

% defining first 30 samples of the 31 samples of y2[n] 
y2_n_display = y2_n(1:N);

% showing the signals s[n],v[n],x[n],y2[n] on 4 subplots - Discrete Time Domain
figure('Name','Time Domain Decomposition'); % new figure window

% s[n] - plot - slow cosine
subplot(2,2,1); 
plot(n,s_n,'b-o','LineWidth',2);
title('1. Desired Signal s[n] (Low Freq)','FontSize',18);
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% v[n] - plot - fast sine
subplot(2,2,2);
plot(n,v_n,'k-o','LineWidth',2);  
title('2. Interference v[n] (High Freq)','FontSize',18);
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% x[n] - plot - input signal to filter
subplot(2,2,3);
plot(n,x_n,'r-o','LineWidth',2);   
title('3. Total Input x[n] = s[n] + v[n]','FontSize',18); 
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% y2[n] - plot - output signal from filter
subplot(2,2,4);
plot(n,y2_n_display,'m-o','LineWidth',2);   
title('4. Filtered Output y_2[n]','FontSize',18); 
xlabel('n [sec]','FontSize',16); 
ylabel('Amplitude [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;


% -------------------------- Question 2A-C --------------------------------

% these sections did not require MATLAB coding.
% in the meantime added 4 variables from data2025.mat file:
% x_test, size: 1x123482
% y,      size: 1x167580
% y_test, size: 1x167580
% y_z,    size: 1x167580
load('data2025.mat');

% -------------------------- Question 2D ----------------------------------
% preparations for DFT - zero padding the signals for linear convolutions
L_min_required = length(y_test); % assuming the output is already...->
% ...-> includes all limitations for linear convolution
N_FFT = 2^nextpow2(L_min_required); % next power of 2 rounding

% creating DFT for all the signals
Y_z_k = fft(y_z, N_FFT); % Y_z[k] - noise
Y_test_k = fft(y_test, N_FFT); % Y_test[k]
X_test_k = fft(x_test, N_FFT); % X_test[k]

% calculating H_total in Frequency domain
Y_test_clean_k = Y_test_k - Y_z_k; % noise "clean-up"
H_total_k = Y_test_clean_k ./ (X_test_k + eps); % eps = tiny epsilon is..->
% ..-> making sure we are not dividing by zero

% building the reconstruction filter H_rec[k]
H_rec_k = 1 ./ H_total_k;

% IFFT{H_total} to get h_total in time domain
h_rec_n = ifft(H_rec_k, 'symmetric');

% limiting time to 1 sec (0.45 + 0.45 = 0.9 < 1)
time_limit = 1.0;
Fs = 44100; % sampling frequency
samples_to_display = round(time_limit * Fs); % number of samplings to show  
time_axis_sec = (0:length(h_rec_n)-1) / Fs; % defining the time axis
n1 = time_axis_sec(1:samples_to_display);
h_total_n_1 = h_rec_n(1:samples_to_display);

% presnting the h_total[n] graph
figure; % new figure window
plot(n1, h_total_n_1, 'b', 'LineWidth', 2);
title('Reconstruction System Impulse Response h_{rec}[n]');
xlabel('Time [sec]');
ylabel('Amplitude');
set(gca, 'FontSize', 20); % defining scale of axis size
grid on;

% -------------------------- Question 2E ----------------------------------
% reconstracting x[n] from y[n] recording
L_max = max(length(y), length(y_test)); % maximum length for zero pad
N_FFT_2 = 2 ^ nextpow2(L_max); % defining FFT length for efficiency

% creating DFT for all the signals
Y_z_k_new = fft(y_z, N_FFT_2); % Y_z[k] - noise
Y_test_k_new = fft(y_test, N_FFT_2); % Y_test[k]
X_test_k_new = fft(x_test, N_FFT_2); % X_test[k]

% calculating H_total in Frequency domain
Y_test_clean_k_new = Y_test_k_new - Y_z_k_new; % noise "clean-up"
H_total_k_new = Y_test_clean_k_new ./ (X_test_k_new + eps); % eps = tiny epsilon is..->
% ..-> making sure we are not dividing by zero

% building the reconstruction filter H_rec[k]
H_rec_k_new = 1 ./ (H_total_k_new + eps);

% defining Y_0[k]
Y_k = fft(y,N_FFT_2); % Y^d[k] = DFT{y[n]}
Y_clean_k = Y_k - Y_z_k_new; % noise "clean-up" from y[n]

% reconstracting x[n]
X_rec_k = Y_clean_k .* H_rec_k_new; % reconstruction in frequency domain
x_rec_n = ifft(X_rec_k, 'symmetric'); % reconstruction in time domain

% limiting time to 3.8 sec
time_limit_new = 3.8;
samples_to_display_new = round(time_limit_new * Fs); % number of samplings to show 
time_axis_sec_new = (0:length(y)-1) / Fs; % defining the time axis
n2 = time_axis_sec_new(1:samples_to_display_new);
x_rec_n_2 = x_rec_n(1:samples_to_display_new);

% presnting the h_total[n] graph
figure; % new figure window
plot(n2, x_rec_n_2, 'b', 'LineWidth', 2);
title('Restored Original Signal x[n] from y[n]');
xlabel('Time [sec]');
ylabel('Amplitude');
set(gca, 'FontSize', 20); % defining scale of axis size
grid on;

% -------------------------- Question 2F ----------------------------------
disp('Playing restored signal x_rec[n]...');
soundsc(x_rec_n_2, Fs);


