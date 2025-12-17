% Appendix A: MATLAB Code for Magnus Coefficient and Spin Ratio Curve Fitting
% (Reconstructed from report appendix image)

clear; clc; close all;

% Spin ratio (Sp) and Magnus coefficient (Cm) data
Sp = [0.018 0.027 0.023 0.026 0.036 0.022 0.015 0.037 0.038 0.027 0.028 0.016];
Cm = [0.015 0.022 0.018 0.021 0.029 0.018 0.012 0.030 0.031 0.022 0.022 0.012];

% Quadratic curve fit: Cm = a*Sp^2 + b*Sp + c
p = polyfit(Sp, Cm, 2);

Cm_fit = polyval(p, Sp);

% Plot data and fitted curve
figure; hold on;
scatter(Sp, Cm, 'filled');
plot(Sp, Cm_fit, 'LineWidth', 2);
xlabel('Spin Ratio (Sp)');
ylabel('Magnus Coefficient (Cm)');
grid on;

% Display curve-fit equation in the Command Window
fprintf('Curve fit: Cm = %.4f*Sp^2 + %.4f*Sp + %.4f\n', p(1), p(2), p(3));
