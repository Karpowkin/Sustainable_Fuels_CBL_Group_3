%% Info
% Toerental: 1500 RPM
% Resolutie van 0.2º CA
% Data_analysis_s2: same as Data_analysis_init, but runs every operating point
% of one test, pegs to the measured intake pressure and saves the results.
%
%% init
clear; clc;close all;
addpath( "Functions","Nasa");
%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv
Runiv = 8.314;
[SpS,El]        = myload('Nasa\NasaThermalDatabase.mat',{'Diesel','O2','N2','CO2','H2O'});
%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 0;
Vd   = pi/4*Cyl.Bore^2*Cyl.Stroke;   % NEW: displacement volume [m^3]
Nrpm = 1500;                         % NEW: engine speed [rpm]
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355; % Intake valve opening
CaIVC = -135; % Intake valve closing
CaEVO = 149; % Exhaust valve opening
CaEVC = -344; % Exhaust valve closing

%% test settings 
testName = 'Test1';
fuelName = 'B7 diesel';
LHV_J_kg = 43.1e6;      % fuel lower heating value [J/kg] (JEC)
dataDir  = fullfile('Data', testName);
% fast data file, slow data file, load (bar IMEP), SOI (deg bTDC)
points = { ...
  '20261123_0000005_group 3_1.5 bar_14fdaq.txt',     '20261123_0000005_group 3_1.5 bar_sdaq.txt',       1.5, 14; ...
  '20261123_0000006_group 3_1.5 bar_16soi_fdaq.txt', '20261123_0000006_group 3_1.5 bar_16soi_sdaq.txt', 1.5, 16; ...
  '20261123_0000007_group 3_1.5 bar_18soi_fdaq.txt', '20261123_0000007_group 3_1.5 bar_18soi_sdaq.txt', 1.5, 18; ...
  '20261123_0000008_group 3_1.5 bar_20soi_fdaq.txt', '20261123_0000008_group 3_1.5 bar_20soi_sdaq.txt', 1.5, 20; ...
  '20261123_0000009_group 3_3.0 bar_14soi_fdaq.txt', '20261123_0000009_group 3_3.0 bar_14soi_sdaq.txt', 3.0, 14; ...
  '20261123_0000011_group 3_3.0 bar_16soi_fdaq.txt', '20261123_0000011_group 3_3.0 bar_16soi_sdaq.txt', 3.0, 16; ...
  '20261123_0000012_group 3_3.0 bar_18soi_fdaq.txt', '20261123_0000012_group 3_3.0 bar_18soi_sdaq.txt', 3.0, 18; ...
  '20261123_0000014_group 3_4.5 bar_14soi_fdaq.txt', '20261123_0000014_group 3_4.5 bar_14soi_sdaq.txt', 4.5, 14};
iselectPoint = 1;          % point used for the single point plots
showDiagnostics = false;   
maxdIMEP = 0.1;            % allowed deviation of IMEP from the set load (bar)
nRange   = [1.25 1.35];    % expected polytropic exponent of compression 

nP = size(points,1);
names = {'Load','SOI_set','SOI_meas','IMEP','COV_IMEP','n_comp','CA10','CA50','CA50_std', ...
         'CA90','BurnDur','IgnDelay','pmax','pmax_cycmax','CA_pmax','Qtot_frac','Eff_i','Flag'};
R = nan(nP, numel(names));
CCall = cell(nP,1);

for j = 1:nP
fprintf('\n===== Point %d: %.1f bar, SOI %d =====\n', j, points{j,3}, points{j,4});
%% Load data 
FullName = fullfile(dataDir, points{j,1});
MatName  = strrep(FullName, '.txt', '.mat');

txtInfo = dir(FullName);
matInfo = dir(MatName);
if ~isempty(matInfo) && matInfo.datenum >= txtInfo.datenum
    load(MatName,'dataIn');                  % fast path: cached binary
else
    fid = fopen(FullName,'r');
    dataIn = fscanf(fid,'%f',[3 Inf]).';     % fast path for plain numeric text
    fclose(fid);
    save(MatName,'dataIn');                  % cache for next run
end
[Nrows,Ncols] = size(dataIn);                 % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p               = reshape(dataIn(:,2),[],Ncycles)*bara; 
injection_current     = reshape(dataIn(:,3),[],Ncycles)*bara;

%% Load Slowdaq
FullNameSlowdaq = fullfile(dataDir, points{j,2});
MatNameSlowdaq  = strrep(FullNameSlowdaq, '.txt', '.mat');

txtInfoSlowdaq = dir(FullNameSlowdaq);
matInfoSlowdaq = dir(MatNameSlowdaq);

if ~isempty(matInfoSlowdaq) && matInfoSlowdaq.datenum >= txtInfoSlowdaq.datenum
    load(MatNameSlowdaq,'dataSlowdaq');
else
    fid = fopen(FullNameSlowdaq,'r');
    dataSlowdaq = fscanf(fid,'%f',[5 Inf]).';
    fclose(fid);
    save(MatNameSlowdaq,'dataSlowdaq');
end

fuelMassFlow_g_s = dataSlowdaq(:,1);
intakeTemperature_C = dataSlowdaq(:,2);
exhaustTemperature_C = dataSlowdaq(:,3);
intakePressure_Pa = dataSlowdaq(:,4)*bara;
airMassFlow_g_s = dataSlowdaq(:,5);

fuelFlow_kg_s = mean(fuelMassFlow_g_s)*g;   % averaged over the measurement
airFlow_kg_s = mean(airMassFlow_g_s)*g;     % averaged over the measurement

% SOI of this point 
CaSOI = -points{j,4};
inj = mean(injection_current,2);
iw  = find(Ca(:,1) > -40 & Ca(:,1) < 10);              % search window around TDC
iS  = iw(find(inj(iw) > 0.2*max(inj(iw)),1));          % first rise above 20% of the peak
SOI_meas = -Ca(iS,1);                                  % measured SOI (deg bTDC)

% Plotting
if j == iselectPoint
f1=figure(1);
set(f1,'Position',[200 800 1200 400]);
plot(Ca,p/bara,'LineWidth',1);
hold on
plot(Ca,injection_current/bara,'LineWidth',1);     
xlabel('Ca');ylabel('p [bar]');
xlim([-360 360]);ylim([0 80]);
iselect = 10;
plot(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b');
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r');
line([CaSOI CaSOI],YLIM,'LineWidth',1,'Color','y');
set(gca,'XTick',-360:60:360,'XGrid','on','YGrid','on');
title('All cycles in one plot.')
hold off
end
%% pV-diagram
V = CylinderVolume(Ca,Cyl);

% Indicated work (average of all cycles)
VClosed = [V; V(1,:)];
pClosed = [p; p(1,:)];

WiCycles = sum(0.5 .* (pClosed(1:end-1,:) + pClosed(2:end,:)) .* ...
    diff(VClosed,1,1), 1);    % [J]
WiAverage = mean(WiCycles);

fprintf('Average indicated work over %d cycles = %.2f J\n', ...
    Ncycles, WiAverage);
IMEP     = WiAverage/Vd;                          % [Pa]
COV_IMEP = std(WiCycles)/mean(WiCycles)*100;      % cycle-to-cycle variation [%]

% Engine KPI 
brakePower_W = WiAverage*Nrpm/120;   % indicated power [W]
CO2Flow_kg_s = NaN %CO2 mass flow rate [kg/s] 
NOxFlow_kg_s = NaN %NOx mass flow rate [kg/s] 

KPI = CalculateKPI(brakePower_W, fuelFlow_kg_s, ...
    LHV_J_kg, CO2Flow_kg_s, NOxFlow_kg_s);

disp(KPI);
% Maximum pressure across all cycles
[maxPressure_Pa, maxIdx] = max(p(:));
[maxCaIdx, maxCycle] = ind2sub(size(p), maxIdx);

fprintf('Maximum pressure across all cycles = %.2f bar (cycle %d at %.1f CA)\n', ...
    maxPressure_Pa/bara, maxCycle, Ca(maxCaIdx, maxCycle));


% combustion characterisation
% check that the variables used below exist
needed = {'Ca','p','V','SpS','CaIVC','CaEVO','CaSOI','bara','dm','Ncycles','Runiv'};
for k = 1:numel(needed)
    if ~exist(needed{k},'var'), warning('Missing variable: %s', needed{k}); end
end
assert(isequal(size(V), size(Ca)), 'V must have the same size as Ca');

% settings 
SOI_run    = points{j,4};   % (deg bTDC) of this run, from the test list
usePegging = true;
Opt.CaIVC   = CaIVC;
Opt.CaEVO   = CaEVO;
if ~isnan(SOI_run), Opt.CaSOI = -SOI_run; else, Opt.CaSOI = CaSOI; end
Opt.CaPeg   = -180;
Opt.TIVC    = 330;
Opt.FiltWin = 4;
Opt.Gamma   = 'nasa';     % or constant
Opt.pIntake = mean(intakePressure_Pa);  % measured intake pressure avereged
Opt.PegMode = 'ambient';      % peg to the measured intake pressure at BDC
Opt.CaPoly  = [-110 -50];     % only used for polytropic
Opt.kPoly   = 1.32;           % only used for polytropic

Cam = Ca(:,1);  pm = mean(p,2);  Vm = V(:,1);
iB  = abs(Cam + 180) <= 3;
idx = Cam >= -100 & Cam <= -30;
pPeg = pm - mean(pm(iB)) + Opt.pIntake;   % pegged to the measured intake pressure

%% ===== DIAGNOSTICS =====
if showDiagnostics
% 1) Volume function check
fprintf('\n--- DIAGNOSTICS ---\n');
fprintf('V(0)    = %.1f cm^3 (expect 35.2)\n', interp1(Cam,Vm,0)*1e6);
fprintf('V(-180) = %.1f cm^3 (expect 757)\n',  interp1(Cam,Vm,-180)*1e6);

% 2) Raw pressure at BDC and position of the pressure peak
fprintf('raw p at BDC = %.2f bar (expect about 1.0 absolute)\n', mean(pm(iB))/bara);
[~,ip] = max(pm);  fprintf('pmax at %.1f deg\n', Cam(ip));

% 3) Compression exponent, raw and pegged
c = polyfit(log(Vm(idx)), log(pm(idx)),   1);  n_raw = -c(1);
c = polyfit(log(Vm(idx)), log(pPeg(idx)), 1);  n_peg = -c(1);
fprintf('n raw = %.2f | n pegged = %.2f (expect 1.25-1.35)\n', n_raw, n_peg);

% 4) Compare: no pegging, ambient pegging, polytropic pegging
OptT = Opt;
OptT.PegMode = 'none';       CC0 = CombustionAnalysis(Cam,pm,Vm,SpS,OptT);
OptT.PegMode = 'ambient';
CC1 = CombustionAnalysis(Cam,pm,Vm,SpS,OptT);
OptT.PegMode = 'polytropic'; CC2 = CombustionAnalysis(Cam,pm,Vm,SpS,OptT);
base = @(CC) mean(CC.aROHR(CC.Ca >= -60 & CC.Ca <= -20));
fprintf('%-11s | aROHR0 | CA50 | CA90 | Qtot\n','peg');
fprintf('%-11s | %6.2f | %4.1f | %4.1f | %4.0f\n','none',      base(CC0),CC0.CA50,CC0.CA90,CC0.Qtotal);
fprintf('%-11s | %6.2f | %4.1f | %4.1f | %4.0f\n','ambient',   base(CC1),CC1.CA50,CC1.CA90,CC1.Qtotal);
fprintf('%-11s | %6.2f | %4.1f | %4.1f | %4.0f\n','polytropic',base(CC2),CC2.CA50,CC2.CA90,CC2.Qtotal);
fprintf('Rair = %g | gamma %.3f..%.3f\n', CC2.Rair, min(CC2.gamma), max(CC2.gamma));

% 5) Log-log pV against an ideal polytrope (n = 1.30, intake pressure at BDC)
figure(6); loglog(Vm/dm^3, pPeg/bara,'b'); hold on
i0 = find(Cam >= -180,1);
loglog(Vm/dm^3, Opt.pIntake/bara*(Vm(i0)./Vm).^1.30, 'k--'); grid on
xlabel('V [dm^3]'); ylabel('p [bar]'); legend('measured (pegged)','n = 1.30');
title('Log-log pV: compare slope with the ideal line'); hold off

% 6) SOI from the injector signal vs the value used
fprintf('SOI from injector signal = %.1f deg | Opt.CaSOI = %.1f deg\n', -SOI_meas, Opt.CaSOI);

% 7) TDC phasing test
fprintf('\n--- TDC shift test (ambient peg) ---\n');
OptS = Opt; OptS.PegMode = 'ambient';
for sh = -2:1:2
    Vs  = CylinderVolume(Cam + sh, Cyl);
    cc  = polyfit(log(Vs(idx)), log(pPeg(idx)), 1);
    Cs  = CombustionAnalysis(Cam,pm,Vs,SpS,OptS);
    fprintf('shift %+d deg: n = %.2f | pre-ign aROHR %.2f | CA50 %.1f | CA90 %.1f | Qtot %.0f J\n', ...
        sh, -cc(1), mean(Cs.aROHR(Cs.Ca >= -60 & Cs.Ca <= -20)), Cs.CA50, Cs.CA90, Cs.Qtotal);
end
fprintf('--- END DIAGNOSTICS ---\n\n');
end

% Mean cycle
pMean = mean(p, 2);
CC = CombustionAnalysis(Ca(:,1), pMean, V(:,1), SpS, Opt);
fprintf('\nMean cycle: CA10 = %.1f, CA50 = %.1f, CA90 = %.1f deg ATDC\n', CC.CA10, CC.CA50, CC.CA90);
fprintf('Burn duration = %.1f deg, ignition delay = %.1f deg, Qtotal = %.0f J\n', ...
        CC.BurnDuration, CC.IgnitionDelay, CC.Qtotal);
fprintf('pmax = %.1f bar at %.1f deg ATDC\n', CC.pMax/bara, CC.CApMax);

% Per cycle
CA50cyc = nan(Ncycles,1);  pMaxcyc = nan(Ncycles,1);
for i = 1:Ncycles
    Ci = CombustionAnalysis(Ca(:,i), p(:,i), V(:,i), SpS, Opt);
    CA50cyc(i) = Ci.CA50;
    pMaxcyc(i) = Ci.pMax;
end
fprintf('CA50 per cycle: mean %.2f, std %.2f deg | pmax per cycle: mean %.1f, std %.2f bar\n', ...
        mean(CA50cyc), std(CA50cyc), mean(pMaxcyc)/bara, std(pMaxcyc)/bara);

fprintf('Ambient peg offset = %.3f bar | p at BDC = %.2f bar\n', ...
    CC.pegOffset/bara, (interp1(Cam,pm,-180) + CC.pegOffset)/bara);

%% data quality checks and results of this point
c = polyfit(log(Vm(idx)), log(pPeg(idx)), 1);  n_peg = -c(1);      % polytropic exponent of compression
Qfrac = CC.Qtotal/(fuelFlow_kg_s/(Nrpm/120)*LHV_J_kg);              % heat released
eta = KPI.BrakeEfficiency;                                          % indicated efficiency here

flag = 0;     % 1 = load deviates from set value, 2 = n outside expected range
if abs(IMEP/bara - points{j,3}) > maxdIMEP, flag = 1; end
if n_peg < nRange(1) || n_peg > nRange(2),  flag = 2; end

R(j,:) = [points{j,3}, points{j,4}, SOI_meas, IMEP/bara, COV_IMEP, n_peg, CC.CA10, CC.CA50, std(CA50cyc), ...
          CC.CA90, CC.BurnDuration, CC.IgnitionDelay, CC.pMax/bara, max(pMaxcyc)/bara, ...
          CC.CApMax, Qfrac, eta, flag];
CCall{j} = CC;
fprintf('IMEP %.2f bar, COV %.1f%%, n %.3f, eta %.3f, flag %d\n', IMEP/bara, COV_IMEP, n_peg, eta, flag);

% plots of the selected point
if j == iselectPoint
figure(3)
plot(CC.Ca, CC.aROHR)
xlabel('Ca [deg]');
ylabel('aROHR [J/deg]')
xlim([-30 90]);
grid on
title('Apparent rate of heat release')

figure(4)
plot(CC.Ca, CC.xb)
xlabel('Ca [deg]');
ylabel('Mass fraction burned [-]')
xlim([-30 90]);
grid on
title('Cumulative heat release')

figure(5)
plot(1:Ncycles, CA50cyc, 'o')
xlabel('Cycle [-]');
ylabel('CA50 [deg ATDC]');
grid on
title('CA50 per cycle')
end
end   % end of the loop over the points

% saving results
if ~exist('Results','dir'), mkdir('Results'); end
writetable(array2table(R, 'VariableNames', names), fullfile('Results', ['Results_' testName '.csv']));

% best SOI per load
loads = unique(R(:,1))';
fprintf('\nBest SOI per load (%s):\n', fuelName);
for L = loads
    k  = find(R(:,1) == L);
    ok = k(R(k,18) == 0);
    if isempty(ok), fprintf('%.1f bar: no valid points\n', L); continue; end
    [~, b] = max(R(ok,17));  b = ok(b);
    fprintf('%.1f bar: SOI %d deg bTDC (eta %.3f, CA50 %.1f deg ATDC, pmax %.1f bar)\n', ...
        L, R(b,2), R(b,17), R(b,8), R(b,13));
    % file for SOI_simple.m (columns: SOI, CA50, pmax, efficiency)
    writematrix(R(k,[2 8 13 17]), fullfile('Results', sprintf('SOI_results_%s_%.1fbar.csv', testName, L)));
end

%plots over all points
mk = {'o-','s-','^-'};
lg = arrayfun(@(L) sprintf('%.1f bar', L), loads, 'UniformOutput', false);

figure(7); hold on
for q = 1:numel(loads), k = R(:,1)==loads(q); plot(R(k,2), R(k,8), mk{q}); end
xlabel('SOI [deg bTDC]'); ylabel('CA50 [deg ATDC]'); grid on; legend(lg, 'Location','best')
title(['CA50 vs SOI, ' fuelName]); hold off

figure(8); hold on
for q = 1:numel(loads), k = R(:,1)==loads(q); plot(R(k,2), R(k,13), mk{q}); end
xlabel('SOI [deg bTDC]'); ylabel('Peak pressure [bar]'); grid on; legend(lg, 'Location','best')
title(['Peak pressure vs SOI, ' fuelName]); hold off

figure(9); hold on
for q = 1:numel(loads), k = R(:,1)==loads(q); plot(R(k,2), 100*R(k,17), mk{q}); end
xlabel('SOI [deg bTDC]'); ylabel('Indicated efficiency [%]'); grid on; legend(lg, 'Location','best')
title(['Indicated efficiency vs SOI, ' fuelName]); hold off

figure(10); hold on
for q = 1:numel(loads), k = R(:,1)==loads(q); plot(R(k,8), 100*R(k,17), mk{q}); end
xlabel('CA50 [deg ATDC]'); ylabel('Indicated efficiency [%]'); grid on; legend(lg, 'Location','best')
title(['Indicated efficiency vs CA50, ' fuelName]); hold off

figure(11); hold on
k = find(R(:,1)==loads(1))';
for j = k, plot(CCall{j}.Ca, CCall{j}.aROHR); end
xlim([-30 60]); xlabel('Ca [deg]'); ylabel('aROHR [J/deg]'); grid on
legend(arrayfun(@(s) sprintf('SOI %d', s), R(k,2), 'UniformOutput', false), 'Location','best')
title(sprintf('aROHR at %.1f bar IMEP, %s', loads(1), fuelName)); hold off

figure(12); hold on
for q = 1:numel(loads), k = R(:,1)==loads(q); plot(R(k,2), R(k,11), mk{q}); end
xlabel('SOI [deg bTDC]'); ylabel('Burn duration CA10-CA90 [deg]'); grid on; legend(lg, 'Location','best')
title(['Burn duration vs SOI, ' fuelName]); hold off
