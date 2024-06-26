%derived from https://cecas.clemson.edu/cvel/Reports/CVEL-14-058.pdf
%for plane wave isotropic material etc etc etc
%all material property inputs are their relative values (ie complex
%dielectric constant
clear
constants
atmosphereprops


%% Radome Matl properties (relative)
epsdash = 4.25;
epsddash = -1.25;
mudash =1;
muddash = 0;
% sig_rd = epsddash*eps_0*5e9*2*pi;
sig_rd = 2.27E-07;

t = 5; %mm

eps_rd = eps_0*(epsdash + 1i*epsddash);
mu_rd = mu_0*(mudash + 1i*muddash);
DP_rd = [eps_rd mu_rd sig_rd];
%% PLot SE for fixed thickness, DPs for a range of frequencies,
hzmin = 8e9;
hzmax = 12e9;
shteps = 100;
hzrange = linspace(hzmin, hzmax, shteps);
freqrange  = linspace(hzmin*2*pi, hzmax*2*pi, shteps);

% Initialize arrays to store the results
SE_values = zeros(length(freqrange), 1);
A_values = zeros(length(freqrange), 1);
R_values = zeros(length(freqrange), 1);
B_values = zeros(length(freqrange), 1);
sigma_values = zeros(length(freqrange), 1);

% Evaluate the function for each omega value in the range
for i = 1:shteps
    freq = freqrange(i);
    % sig_rd = epsddash * eps_0 * freq;
    % sigma_values(i) = sig_rd;
    % DP_rd = [eps_rd mu_rd sig_rd];

    [SE, A, R, B] = shieldingeffectiveness(freq, DPs1, DP_rd, DPs2, t);
    SE_values(i) = SE;
    A_values(i) = A;
    R_values(i) = R;
    B_values(i) = B;
end


% figure
% semilogx(hzrange, SE_values);
% hold on
% semilogx(hzrange, A_values);
% semilogx(hzrange, R_values);
% semilogx(hzrange, B_values);
% xlabel("Frequency (Hz)")
% ylabel("SE (dB)")
% legend('SE','A', "R", "B");
% title('')
% hold off


figure
plot(hzrange/1e9, SE_values);
hold on
plot(hzrange/1e9, A_values);
plot(hzrange/1e9, R_values);
plot(hzrange/1e9, B_values);
legend('SE','A', "R", "B");
xlabel("Frequency (GHz)")
ylabel("SE (dB)")
legend('SE','A', "R", "B");
title("Shielding Effectivness over Frequency for \epsilon' = " + num2str(epsdash) + ', \epsilon" = ' + num2str(epsddash) + ', t = ' + num2str(t) + 'mm')
hold off

% plot(hzrange/1e9, sigma_values);




