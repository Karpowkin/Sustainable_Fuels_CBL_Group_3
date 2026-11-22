%% SOI_simple.m - find the optimal injection timing (SOI) with an SOI sweep
% SOI in degrees before TDC (DBTDC), positive numbers. In the combustion code Ca = -SOI.
% Run this script on its own (it starts with clear). Put it next to Data_analysis_init.m.
% Do ONE sweep per load (30 / 50 / 70 %): change loadName and run again.
clear; clc; close all;

loadName = '50 % load';          % only used in the plot titles

%% PART 1: list of SOIs to test (give this to the lab)
SOI_min    = 14;     % latest SOI of the sweep [DBTDC]
SOI_max    = 21;     % earliest SOI allowed [DBTDC] - earlier than 21 is NOT allowed
Pmax_limit = 100;    % max allowed peak pressure [bar] (ask the lab for this year's limit)

step       = 2;      % step size 2 or 3 CAD
includeMax = false;  % true = also test the limit itself (21) if the steps do not land on it

SOI_test = SOI_min:step:SOI_max;
if includeMax && SOI_test(end) < SOI_max
    SOI_test(end+1) = SOI_max;
end
fprintf('SOIs to run in the lab [DBTDC]: %s\n', mat2str(SOI_test));

%% PART 2: results from the data analysis (fill in after the lab)
% One value per SOI in SOI_test, SAME ORDER. Example of how to fill in (4 SOIs for 14:2:20):
%   Eff  = [0.30 0.31 0.32 0.31];   % efficiency [-]
%   Pmax = [ 55   60   65   72 ];   % peak pressure [bar]
%   CA50 = [ 16   14   12   10 ];   % CA50 [deg after TDC]
Eff  = NaN(size(SOI_test));
Pmax = NaN(size(SOI_test));
CA50 = NaN(size(SOI_test));
 
% (columns: SOI, CA50, Pmax). Efficiency still has to be typed in above.
useCSV = false;
if useCSV
    T = sortrows(readmatrix('SOI_results.csv'), 1);
    SOI_test = T(:,1)';
    CA50     = T(:,2)';
    Pmax     = T(:,3)';
    Eff      = NaN(size(SOI_test));      % type Eff = [...] below this line, same order as SOI_test
end

if any(isnan(Eff)) || any(isnan(Pmax)) || any(isnan(CA50))
    disp('Fill in Eff, Pmax and CA50 in PART 2 (same length as SOI_test), then run again.');
    return
end
assert(numel(Eff)==numel(SOI_test) && numel(Pmax)==numel(SOI_test) && numel(CA50)==numel(SOI_test), ...
    'Eff, Pmax and CA50 must have the same number of values as SOI_test (%d).', numel(SOI_test));

%% PART 3: find the best SOI
Eff_ok = Eff;
Eff_ok(Pmax > Pmax_limit) = NaN;             % runs above the pressure limit are not allowed

[bestEff, i] = max(Eff_ok);                  % highest efficiency of the allowed runs (NaN is ignored)
if isnan(bestEff)
    error('All runs exceed Pmax_limit (%g bar) - no allowed SOI.', Pmax_limit);
end
SOI_best = SOI_test(i);

fprintf('Best SOI: %.1f DBTDC (efficiency %.3f, Pmax %.0f bar, CA50 %.1f deg ATDC)\n', ...
    SOI_best, bestEff, Pmax(i), CA50(i));

%% PART 4: fine sweep around the best SOI (run this next in the lab)
fineStep = 0.5;                              % check that the ECU can really be set in 0.5 deg steps
SOI_fine = (SOI_best - 1.5):fineStep:(SOI_best + 1.5);
fprintf('Fine sweep SOIs [DBTDC]: %s\n', mat2str(SOI_fine));
% Put SOI_fine in place of SOI_test, fill in the new results, and run again.

%% PART 5: plots
figure;
subplot(1,3,1);
plot(SOI_test, Eff, 'o-'); hold on;
plot(SOI_best, bestEff, 'r*', 'MarkerSize', 12);
xlabel('SOI [DBTDC]'); ylabel('Efficiency [-]'); grid on;
title(['Efficiency vs SOI, ' loadName]);

subplot(1,3,2);
plot(SOI_test, Pmax, 'o-'); hold on;
yline(Pmax_limit, 'r--', 'limit');
xlabel('SOI [DBTDC]'); ylabel('Peak pressure [bar]'); grid on;
title('Peak pressure vs SOI');

subplot(1,3,3);
plot(SOI_test, CA50, 'o-');
xlabel('SOI [DBTDC]'); ylabel('CA50 [deg ATDC]'); grid on;
title('CA50 vs SOI');