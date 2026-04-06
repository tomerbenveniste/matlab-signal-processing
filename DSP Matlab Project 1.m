clc;
clear all;
close all;
 
% -------------------------- question 1A ----------------------------------
t = 0.2:0.01:3; %continuous time vector
Wm = 3*pi; 
A = 4./(Wm*pi*(t.^2)); % signal amplitude
x_t = A.*((sin(Wm*t)).^2).*(cos(Wm*t)).*(sin(2*Wm*t)); %continuous signal
x_t_abs = abs(x_t); % absolute value of the signal

% showing the signal |x(t)|
figure; % new figure window
plot(t,x_t_abs,'b','LineWidth',2);  % continuous signal graph
title('Signal |x(t)| - time domain - Question 1A','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('|x(t)| [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;


% -------------------------- question 1B ----------------------------------
w = -17*pi:0.01:17*pi; % frequency vector 
% defining the triangles which assemble X(w)
triangle1 = tripuls(w + 3*Wm, 4*Wm); 
triangle2 = -tripuls(w - 3*Wm, 4*Wm);
triangle3 = tripuls(w + Wm, 4*Wm);
triangle4 = -tripuls(w - Wm, 4*Wm);

%showing all the triangular pulses in one plot
figure; %new figure window
plot(w,triangle1,'r','LineWidth',2);
hold on;
plot(w,triangle2,'g','LineWidth',2);  
hold on;
plot(w,triangle3,'m','LineWidth',2);  
hold on;
plot(w,triangle4,'c','LineWidth',2);  
title('all the triangular pulses','FontSize',18); 
xlabel('w [rad/sec]','FontSize',16); 
ylabel('triangular pulses','FontSize',16);
legend([{'+3Wm'};{'-3Wm'};{'+Wm'};{'-Wm'}],'FontSize',16);
% defining ticks for w axis
xticks([-17*pi, -15*pi, -9*pi, -3*pi, 0, 3*pi, 9*pi, 15*pi, 17*pi]); 
labels = {'-17\pi', '-15\pi', '-9\pi', '-3\pi', '0', '3\pi', '9\pi', '15\pi', '17\pi'};
xticklabels(labels);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;

% X(w) is the sum of the triangles above
X_w = abs(triangle1 + triangle2 + triangle3 + triangle4);

% showing the signal |X(w)|
figure; %new figure window
plot(w,X_w,'r','LineWidth',2);
title('Signal |X(w)| - Frequency domain - Question 1B','FontSize',18); 
xlabel('w [rad/sec]','FontSize',16); 
ylabel('|X(w)|','FontSize',16);
% defining ticks for w axis
xticks([-17*pi, -15*pi, -9*pi, -3*pi, 0, 3*pi, 9*pi, 15*pi, 17*pi]); 
labels = {'-17\pi', '-15\pi', '-9\pi', '-3\pi', '0', '3\pi', '9\pi', '15\pi', '17\pi'};
xticklabels(labels);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;


% -------------------------- question 1C ----------------------------------
ws = 10*Wm; % sampling frequency
Ts = (2*pi)/ws; % cycle/period time
tn = 0.2:Ts:3; % discrete time vector
A = 4./(Wm*pi*(tn.^2)); % amplitude with discrete time
x_tn = A.*((sin(Wm*tn)).^2).*(cos(Wm*tn)).*(sin(2*Wm*tn)); % discrete signal

%showing the signal |x(tn)|
figure; %new figure window
plot(t,x_t,'b','LineWidth',2);   
hold on;
stairs(tn, x_tn, '-*r', 'LineWidth', 1.5); % creating ZOH signal with the samplings
title('Signals x(t) & ZOH[x(tn)] - Question 1C','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('x(t) & ZOH[x(tn)] [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'x(t)'};{'ZOH[x(tn)]'}],'FontSize',16);
grid on;


% -------------------------- question 1D ----------------------------------
k_vector = [-1,0,1]; % these are the only k's in the infinite sum that will be in the domain
% X_sum will be the sum of all triangles in the domain
X_sum = zeros(size(w));
for k = k_vector
    triangle_1k = tripuls(w - k*ws + 3*Wm, 4*Wm);
    triangle_2k = -tripuls(w - k*ws - 3*Wm, 4*Wm);
    triangle_3k = tripuls(w - k*ws + Wm, 4*Wm);
    triangle_4k = -tripuls(w - k*ws - Wm, 4*Wm);
    X_sum = X_sum + triangle_1k + triangle_2k + triangle_3k + triangle_4k;
end
% defining Xzoh(w)
e_w = exp(-1j*w*(Ts/2)); % sub function
X_w_ZOH = 1j*e_w.*sinc(w/ws).*X_sum;
% defining abs of Xzoh(w)
X_w_ZOH_abs = abs(X_w_ZOH);

%showing the signal |X_ZOH(w)|
figure; %new figure window
plot(w,X_w_ZOH_abs,'r','LineWidth',2);
title('Signal |XZOH(w)| - Frequency domain - Question 1D','FontSize',18); 
xlabel('w [rad/sec]','FontSize',16); 
ylabel('|XZOH(w)|','FontSize',16);
% defining ticks for w axis
xticks([-17*pi, -15*pi, -9*pi, -3*pi, 0, 3*pi, 9*pi, 15*pi, 17*pi]); 
labels = {'-17\pi', '-15\pi', '-9\pi', '-3\pi', '0', '3\pi', '9\pi', '15\pi', '17\pi'};
xticklabels(labels);
set(gca, 'FontSize', 14); % defining scale of axis size
grid on;


% -------------------------- question 1E ----------------------------------
Window = abs(w) <= ws/2; % creating a boolean window for H filter
H_w = (exp(1j*pi*(w/ws))./sinc(w/ws)).*Window; % defining H(w)
% Xrec(w) = Xzoh(w)H(w) - convolution in time <-> multiplication in frequency
X_rec_w = X_w_ZOH.*H_w;
% defining inverse fourier transform with integral using trapz function
x_rec_t = zeros(size(t));
for i = 1:length(t)
    integrand = X_rec_w.*exp(1j*w*t(i));
    x_rec_t(i) = (1/(2*pi))*trapz(w,integrand);
end

%showing the signal X(t) & X_rec(t)
figure; %new figure window
plot(t,x_t,'k','LineWidth',4);   
hold on;
plot(t,x_rec_t,'--r','LineWidth',2); 
title('Signals X(t) & Xrec(t) - Question 1E','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('X(t) & Xrec(t) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'X(t)'};{'Xrec(t)'}],'FontSize',16);
grid on;

% calculate error
error_avg = mean(abs(x_t - x_rec_t));


% -------------------------- question 1F ----------------------------------
ws_new = 9*Wm; % new sampling frequency that will cause aliasing
Ts_new = (2*pi)/ws_new; % new cycle time
% X_sum will be the sum of all triangles in the domain
X_sum_new = zeros(size(w));
for k = k_vector
    triangle_1k = tripuls(w - k*ws_new + 3*Wm, 4*Wm);
    triangle_2k = -tripuls(w - k*ws_new - 3*Wm, 4*Wm);
    triangle_3k = tripuls(w - k*ws_new + Wm, 4*Wm);
    triangle_4k = -tripuls(w - k*ws_new - Wm, 4*Wm);
    X_sum_new = X_sum_new + triangle_1k + triangle_2k + triangle_3k + triangle_4k;
end
% defining the new Xzoh(w)
e_w_new = exp(-1j*w*(Ts_new/2));
X_w_ZOH_new = 1j*e_w_new.*sinc(w/ws_new).*X_sum_new;
X_w_ZOH_abs_new = abs(X_w_ZOH_new);

%showing the signals |X_ZOH(w)| for both ws frequncies on the same plot
figure; %new figure window
plot(w,X_w_ZOH_abs,'r','LineWidth',2);    
hold on;
plot(w,X_w_ZOH_abs_new,'c','LineWidth',2);
title('Signal |XZOH(w)| in both sampling frequencies','FontSize',18); 
xlabel('w [rad/sec]','FontSize',16); 
ylabel('|XZOH(w)|','FontSize',16);
% defining ticks for w axis
xticks([-17*pi, -15*pi, -9*pi, -3*pi, 0, 3*pi, 9*pi, 15*pi, 17*pi]); 
labels = {'-17\pi', '-15\pi', '-9\pi', '-3\pi', '0', '3\pi', '9\pi', '15\pi', '17\pi'};
xticklabels(labels);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'Ws = 10Wm'};{'Ws = 9Wm'}],'FontSize',16);
grid on;

% repeating 1E section
% redefining everything - can be seen in section 1E
Window_new = abs(w) <= ws_new/2; 
H_w_new = (exp(1j*pi*(w/ws_new))./sinc(w/ws_new)).*Window_new;
X_rec_w_new = X_w_ZOH_new.*H_w_new;
x_rec_t_new = zeros(size(t));
for i = 1:length(t)
    integrand_new = X_rec_w_new.*exp(1j*w*t(i));
    x_rec_t_new(i) = (1/(2*pi))*trapz(w,integrand_new);
end

%showing the signal X(t) & the new X_rec(t) on the same plot
figure; %new figure window
plot(t,x_t,'k','LineWidth',4);   
hold on;
plot(t,x_rec_t_new,'--r','LineWidth',2); 
title('Signals X(t) & Xrec(t) - Ws = 9Wm - Question 1F','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('X(t) & Xrec(t) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'X(t)'};{'Xrec(t) - Ws = 9Wm'}],'FontSize',16);
grid on;

% calculate error vector
error_avg_new = mean(abs(x_t - x_rec_t_new));


% -------------------------- question 2A ----------------------------------
wA = (2*pi)/5;
wB = pi/7;
t = 0 : 0.01 : 70; % continous time vector
x_t = 5*cos(wA*t) - 3*sin(wB*t); % defining x(t)
T = 70; % cycle time
w0 = (2*pi)/T; % the signal frequency
N = 15; % number of samplings required
dt = 70/15; % diffrence between samplings
tn = 0 : dt : (T-dt); % discrete time vector
x_tn = 5*cos(wA*tn) - 3*sin(wB*tn); % defining x(tn)

%showing the signal x(t) & x(tn)
figure; %new figure window
plot(t,x_t,'c','LineWidth',2);   
hold on;
plot(tn,x_tn,'*r','LineWidth',2); 
title('Signals x(t) & x(tn) - Question 2A','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('x(t) & x(tn) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'x(t)'};{'x(tn)'}],'FontSize',16);
grid on;


% -------------------------- question 2B ----------------------------------
n_vec = (0 : N-1).'; % coloumn vector
k_vec = (-7 : 7); % row vector
nk = (n_vec*k_vec); % creating a matrix in nxk
F_mat = exp(1j*((2*pi)/N)*nk); % exp(matrix) -> matrix of exponents
inv_F = inv(F_mat); % creating the inverse(F) matrix
a_vec = inv_F*(x_tn.'); % calculating the coefficients vector


% -------------------------- question 2C ----------------------------------
% creating Xrec(t) as the sum of harmonics as k runs from -7 to 7
x_rec_t = zeros(size(t));
k_vec(-1+8) = 14;
k_vec(1+8) = -14;
% assembling the fourier sum
for i = 1:length(k_vec)
    a_k = a_vec(i);
    k = k_vec(i);
    x_k = a_k * exp(1j*w0*k*t);
    x_rec_t = x_rec_t + x_k;
end

%showing the signal x(t) & Xrec(t)
figure; %new figure window
plot(t,x_t,'k','LineWidth',4);   
hold on;
plot(t,x_rec_t,'--c','LineWidth',2);
hold on;
plot(tn,x_tn,'*r','LineWidth',2); 
title('Signals X(t) & Xrec(t) - Question 2C','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('X(t) & Xrec(t) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'X(t)'};{'Xrec(t)'};{'X(tn)'}],'FontSize',16);
grid on;


% -------------------------- question 2D ----------------------------------
% repeating 2A
% defining random samplings
tn_random = sort(rand(N,1)*T); % discrete random time coloumn vector
% safety loop that checks that all the samplings are far enough from each other
% and also verifies that we did not sample exactly on t=0 and t=70
while min(diff(tn_random)) < 0.1 || (tn_random(1) == 0 && tn_random(15) == 70)
    tn_random = sort(rand(N,1)*T); % discrete random time coloumn vector
end 
% defining x(tn_random)
x_tn_random = 5*cos(wA*tn_random) - 3*sin(wB*tn_random); 

%showing the signal x(t) & x(tn_random)
figure; %new figure window
plot(t,x_t,'k','LineWidth',2);   
hold on;
plot(tn_random,x_tn_random,'*r','LineWidth',2); 
title('Signals x(t) & x(tn random) - Question 2D','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('x(t) & x(tn random) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'x(t)'};{'x(tn random)'}],'FontSize',16);
grid on;

% repeating 2B
tnk = (tn_random*k_vec); % creating a matrix in nxk
F_mat = exp(1j*w0*tnk); % exp(matrix) -> matrix of exponents
inv_F = inv(F_mat); % creating the inverse(F) matrix
a_vec_random = inv_F*(x_tn_random); % calculating the coefficients vector

% repeating 2C
% creating Xrec(t) with random samplings from the fourier sum
x_rec_t_random = zeros(size(t));
% assembling the fourier sum
for i = 1:length(k_vec)
    a_k = a_vec(i);
    k = k_vec(i);
    x_k = a_k * exp(1j*w0*k*t);
    x_rec_t_random = x_rec_t_random + x_k;
end

%showing the signal x(t) & Xrec(t) after random samplings
figure; %new figure window
plot(t,x_t,'k','LineWidth',4);   
hold on;
plot(t,x_rec_t_random,'--c','LineWidth',2);
hold on;
plot(tn_random,x_tn_random,'*r','LineWidth',2); 
title('Signals X(t) & Xrec(t) after random samplings - Question 2D','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('X(t) & Xrec(t) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'X(t)'};{'Xrec(t)'};{'X(tn random)'}],'FontSize',16);
grid on;


% -------------------------- question 2E ----------------------------------
% repeating section 2B - Uniform sample with uncertainty of tn = tn+0.01*rand(1)
tn_2 = zeros(size(tn));
for i = 1:length(tn)
    tn_2(i) = tn(i) + 0.01*rand(1);
end
% calculating new a vector
tn2_k = (tn_2.'*k_vec); % creating a matrix in tnxk
F_mat_2 = exp(1j*w0*tn2_k); % exp(matrix) -> matrix of exponents
inv_F_2 = inv(F_mat_2); % creating the inverse(F) matrix
a_vec_2 = inv_F_2*(x_tn.'); % calculating the coefficients vector

% calculating K(F) - condition number of F
K_F = cond(F_mat_2);

% repeating section 2C - Uniform sample with uncertainty of tn = tn+0.01*rand(1)
% creating Xrec(t) as the sum of harmonics
x_rec2_t = zeros(size(t));
% assembling the fourier sum
for i = 1:length(k_vec)
    a_k = a_vec_2(i);
    k = k_vec(i);
    x_k = a_k * exp(1j*w0*k*t);
    x_rec2_t = x_rec2_t + x_k;
end

%showing the signal x(t) & Xrec2(t)
figure; %new figure window
plot(t,x_t,'k','LineWidth',4);   
hold on;
plot(t,x_rec2_t,'--c','LineWidth',2);
hold on;
plot(tn,x_tn,'*r','LineWidth',2); 
title('Signals X(t) & Xrec(t) with uncertainty of tn = tn+0.01*rand(1) - Question 2E','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('X(t) & Xrec(t) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'X(t)'};{'Xrec(t)'};{'X(tn)'}],'FontSize',16);
grid on;

% calculate error
error_avg = mean(abs(x_t - x_rec2_t));

% repeating 2D which is repeating again 2A,B,C with random samplings
% repeating section 2B - RANDOM sample with uncertainty of tn = tn+0.01*rand(1)
tn_random_2 = zeros(size(tn));
for i = 1:length(tn)
    tn_random_2(i) = tn_random(i) + 0.01*rand(1);
end
% calculating new a vector
tn_random_2_k = (tn_random_2.'*k_vec); % creating a matrix in tnxk
F_mat_random_2 = exp(1j*w0*tn_random_2_k); % exp(matrix) -> matrix of exponents
inv_F_random_2 = inv(F_mat_random_2); % creating the inverse(F) matrix
a_vec_random_2 = inv_F_random_2*(x_tn_random); % calculating the coefficients vector

% calculating K(F) - condition number of F
K_F_2 = cond(F_mat_random_2);

% repeating section 2C - RANDOM sample with uncertainty of tn = tn+0.01*rand(1)
% creating Xrec(t) as the sum of harmonics
x_rec2_random_t = zeros(size(t));
% assembling the fourier sum
for i = 1:length(k_vec)
    a_k = a_vec_random_2(i);
    k = k_vec(i);
    x_k = a_k * exp(1j*w0*k*t);
    x_rec2_random_t = x_rec2_random_t + x_k;
end

%showing the signal x(t) & Xrec2(t)
figure; %new figure window
plot(t,x_t,'k','LineWidth',2);   
hold on;
plot(t,x_rec2_random_t,'--r','LineWidth',2);
hold on;
plot(tn_random,x_tn_random,'*m','LineWidth',2); 
title('Signals X(t) & Xrec(t) - random tn with uncertainty of tn = tn+0.01*rand(1) - Question 2E','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('X(t) & Xrec(t) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'X(t)'};{'Xrec(t)'};{'X(tn)'}],'FontSize',16);
grid on;

% calculate error
error_avg_2 = mean(abs(x_t - x_rec2_random_t));


% -------------------------- question 2F ----------------------------------
% repeating section 2A - 40 RANDOM samplings
% defining 40 random samplings
N_40 = 40;
tn_random_40 = sort(rand(N_40,1)*T); % discrete random time coloumn vector
% safety loop that checks that all the samplings are far enough from each other
% and also verifies that we did not sample exactly on t=0 and t=70
while min(diff(tn_random_40)) < 0.1 || (tn_random_40(1) == 0 && tn_random_40(40) == 70)
    tn_random_40 = sort(rand(N_40,1)*T); % discrete random time coloumn vector
end 
% defining x(tn_random)
x_tn_random_40 = 5*cos(wA*tn_random_40) - 3*sin(wB*tn_random_40); 

%showing the signal x(t) & x(tn_random_40)
figure; %new figure window
plot(t,x_t,'k','LineWidth',2);   
hold on;
plot(tn_random_40,x_tn_random_40,'*r','LineWidth',2); 
title('Signals x(t) & x(tn random) - 40 samplings Question 2F','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('x(t) & x(tn random) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'x(t)'};{'x(tn random)'}],'FontSize',16);
grid on;

% repeating section 2B - 40 RANDOM samplings with uncertainty of tn = tn+0.01*rand(1)
tn_random_for_F_40 = zeros(size(tn_random_40));
for i = 1:length(tn_random_40)
    tn_random_for_F_40(i) = tn_random_40(i) + 0.01*rand(1);
end
% calculating new a vector
tn_random_40_k = (tn_random_for_F_40*k_vec); % creating a matrix in 40x15 (tnk)
F_40 = exp(1j*w0*tn_random_40_k); % exp(matrix) -> matrix of exponents
F_40_T = transpose(F_40); % transposing F

% a = [[(F^T)*F]^(-1)]*(F^T)*x
a_vec_random_40 = ((inv(F_40_T*F_40))*F_40_T)*x_tn_random_40;

% calculating K(F) - condition number of F
K_F_40 = cond(F_40);

% repeating section 2C - 40 RANDOM samplings with uncertainty of tn = tn+0.01*rand(1)
% creating Xrec(t) as the sum of harmonics
x_rec40_random_t = zeros(size(t));
% assembling the fourier sum
for i = 1:length(k_vec)
    a_k = a_vec_random_40(i);
    k = k_vec(i);
    x_k = a_k * exp(1j*w0*k*t);
    x_rec40_random_t = x_rec40_random_t + x_k;
end

%showing the signal x(t) & Xrec40(t)
figure; %new figure window
plot(t,x_t,'k','LineWidth',2);   
hold on;
plot(t,x_rec40_random_t,'--r','LineWidth',2);
hold on;
plot(tn_random_40,x_tn_random_40,'*m','LineWidth',2); 
title('Signals X(t) & Xrec(t) - 40 random samplings with uncertainty of tn = tn+0.01*rand(1) - Question 2F','FontSize',18); 
xlabel('t [sec]','FontSize',16); 
ylabel('X(t) & Xrec(t) [V]','FontSize',16);
set(gca, 'FontSize', 14); % defining scale of axis size
legend([{'X(t)'};{'Xrec(t)'};{'X(tn)'}],'FontSize',16);
grid on;

% calculate error
error_avg_40 = mean(abs(x_t - x_rec40_random_t));