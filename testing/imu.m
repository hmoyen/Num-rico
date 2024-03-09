%  Daniel Kawano, Rose-Hulman Institute of Technology
%  Last modified:  Dec 19, 2017

clear all
close all
clc

%  Parameter values:

R0 = 5;             %  m
omega0 = pi;        %  rad/s
dt = 0.001;         %  s
N = 1;

%  Exact solution:

te = [0 : dt : 2*N]';           %  s
psie = omega0*te;               %  rad
xA1e = R0*cos(psie);            %  m
xA2e = R0*sin(psie);            %  m

figure
set(gcf, 'color', 'w')
plot(xA1e, xA2e, '-r', 'linewidth', 2)
xlabel('\itx_{A\rm_{1}}\rm (m)')
ylabel('\itx_{A\rm_{2}}\rm (m)')
axis([1.2*R0*[-1, 1, -1, 1]])
axis equal

%  Integration drift:

dt = 0.01;                      %  s
N = 100; 
t = [0 : dt : 2*N]';            %  s

%  (1)  Noise-free IMU "data":

psidote = omega0*ones(length(t),1);             %  rad/s
a1e = (-R0*omega0^2)*ones(length(t),1);         %  m/s^2
a2e = 0*a1e;                                    %  m/s^2

figure
set(gcf, 'color', 'w')
subplot(311)
plot(t, psidote, '-b', 'linewidth', 2)
xlabel('Time (s)')
ylabel({['Angular velocity,'], ['\it\omega\rm_{3} (rad/s)']})
xlim([0, 10])
ylim([2, 4])
subplot(312)
plot(t, a1e, '-b', 'linewidth', 2)
xlabel('Time (s)')
ylabel({['Specific force,'], ['\itf\rm_{1} (m/s^2)']})
xlim([0, 10])
ylim([-51, -48])
subplot(313)
plot(t, a2e, '-b', 'linewidth', 2)
xlabel('Time (s)')
ylabel({['Specific force,'], ['\itf\rm_{2} (m/s^2)']})
xlim([0, 10])
ylim([-1, 1])

psie = 0 + cumtrapz(t, psidote);                                %  rad

xAddot1e = cos(psie).*a1e - sin(psie).*a2e;                     %  m/s^2
xAddot2e = sin(psie).*a1e + cos(psie).*a2e;                     %  m/s^2

xA1id = R0 + cumtrapz(t, 0 + cumtrapz(t, xAddot1e));            %  m
xA2id = 0 + cumtrapz(t, R0*omega0 + cumtrapz(t, xAddot2e));     %  m

figure
set(gcf, 'color', 'w')
plot(xA1id, xA2id, '-g', xA1e, xA2e, '-r', 'linewidth', 2)
xlabel('\itx_{A\rm_{1}}\rm (m)')
ylabel('\itx_{A\rm_{2}}\rm (m)')
axis([1.2*R0*[-1, 1, -1, 1]])
axis equal

%  (2)  IMU "data" with noise:

n = 5;

psidotn = add_awgn_noise(psidote, 30);            %  rad/s
a1n = add_awgn_noise(a1e, 30);                    %  m/s^2
a2n = add_awgn_noise(a2e, 30);                    %  m/s^2

figure
set(gcf, 'color', 'w')
subplot(311)
plot(t(1:n*(length(t)-1)/N), psidotn(1:n*(length(t)-1)/N), '-b', 'linewidth', 2)
xlabel('Time (s)')
ylabel({['Angular velocity,'], ['\it\omega\rm_{3} (rad/s)']})
ylim([2, 4])
subplot(312)
plot(t(1:n*(length(t)-1)/N), a1n(1:n*(length(t)-1)/N), '-b', 'linewidth', 2)
xlabel('Time (s)')
ylabel({['Specific force,'], ['\itf\rm_{1} (m/s^2)']})
ylim([-51, -48])
subplot(313)
plot(t(1:n*(length(t)-1)/N), a2n(1:n*(length(t)-1)/N), '-b', 'linewidth', 2)
xlabel('Time (s)')
ylabel({['Specific force,'], ['\itf\rm_{2} (m/s^2)']})
ylim([-1, 1])

psin = 0 + cumtrapz(t, psidotn);                                %  rad

xAddot1n = cos(psin).*a1n - sin(psin).*a2n;                     %  m/s^2
xAddot2n = sin(psin).*a1n + cos(psin).*a2n;                     %  m/s^2

xA1idn = R0 + cumtrapz(t, 0 + cumtrapz(t, xAddot1n));           %  m
xA2idn = 0 + cumtrapz(t, R0*omega0 + cumtrapz(t, xAddot2n));    %  m

figure
set(gcf, 'color', 'w')
plot(xA1idn(1:n*(length(t)-1)/N), xA2idn(1:n*(length(t)-1)/N), '-g', xA1e, xA2e, '-r', 'linewidth', 2)
xlabel('\itx_{A\rm_{1}}\rm (m)')
ylabel('\itx_{A\rm_{2}}\rm (m)')
axis([1.2*R0*[-1, 1, -1, 1]])
axis equal


function [y, n] = add_awgn_noise(x, SNR_dB)
% [y, n] = add_awgn_noise(x, SNR_dB)
% Adds AWGN noise vector to signal 'x' to generate a resulting signal vector 'y'
% of specified SNR in dB. It also returns the noise vector 'n' that is added to the signal x.

L = length(x);
SNR = 10^(SNR_dB / 10); % SNR to linear scale
Esym = sum(abs(x).^2) / L; % Calculate actual symbol energy
N0 = Esym / SNR; % Find the noise spectral density

if isreal(x)
    noiseSigma = sqrt(N0); % Standard deviation for AWGN Noise when x is real
    n = noiseSigma * randn(1, L); % Computed noise
else
    noiseSigma = sqrt(N0 / 2); % Standard deviation for AWGN Noise when x is complex
    n = noiseSigma * (randn(1, L) + 1i * randn(1, L)); % Computed noise
end

y = x + n; % Received signal
end
