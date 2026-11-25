%% Info
% 4CBLB10 Group 3 - engine analysis for ALL test points in one run
% Engine : Hatz 1D90E, 1500 RPM
% Fuel   : B7 diesel (Test 1, 23 Nov 2026)  -> change Fuel.* below for other fuels
% Data   : 0.2 deg CA resolution, 200 cycles per test point
%          all *fdaq.txt / *sdaq.txt files in the folder Actual_Data_set
%
% Sections
%   Settings / units / fuel / geometry / emissions table
%   1. Find all data files (fdaq + matching sdaq by run number)
%   2. Loop over all test points
%        2a load fast + slow data (cached as .mat after the first run)
%        2b indicated work
%        2c combustion characterisation (pegging, CA10/50/90)   [teammate]
%        2d diagnostics (only for one point, if switched on)    [teammate]
%        2e KPIs: efficiency, isFC, isCO2, isNOx via CalculateKPI,
%           plus isCO, isHC, soot and analyser-CO2 check           [NEW]
%   3. Summary table + plots of all points
%   4. Optimal SOI                                               [TODO]
%
% Needs the folders "Functions" (CylinderVolume, CombustionAnalysis,
% CalculateKPI, myload, myfind)
% and "Nasa" next to this script.

%% init
clear; clc; close all;
% Work relative to THIS file, so it runs no matter which folder MATLAB is in
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(fullfile(scriptDir,'Functions'), fullfile(scriptDir,'Nasa'));

%% Units
mm      = 1e-3; dm = 0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;

%% Settings
% Folder with all fdaq/sdaq files. The script looks in these places, in order:
dataDirCandidates = {fullfile(scriptDir,'Actual_Data_set'), fullfile(scriptDir,'..','Actual_Data_set')};
RPM             = 1500;
NdatapointsperCycle = 720/0.2;
nSkip           = 20;     % slow-daq samples skipped at the start (engine still settling)
showDiagnostics = false;  % true = run the diagnostics block for ONE point (prints a lot)
diagPoint       = 5;      % which point (index after sorting, see printed list)

%% Fuel properties (B7, handbook Table 2.1)
Fuel.name = 'B7';
Fuel.LHV  = 43*MJ;   % [J/kg]
Fuel.wC   = 0.858;   % carbon mass fraction [-]
Fuel.wH   = 0.136;
Fuel.wO   = 0.006;

% Molar masses [g/mol]
Mol.C   = 12.011;
Mol.CO2 = 44.01;
Mol.NOx = 46.01;     % NOx reported as NO2 (convention)
Mol.CO  = 28.01;
Mol.HC  = 13.9;      % HC as C1 (CH1.85) - check analyser basis
Mol.exh = 28.8;      % lean diesel exhaust, close to air
rhoExhStd = 1.293;   % exhaust density at 0 degC, 1 atm [kg/m^3] (for soot)

%% Emission analyser values (already averaged, Test 1 sheet)
%            load  SOI   CO2[%] NOx[ppm] Soot[FSN] CO[%] HC[ppm] O2[%]
gasTable = [ 1.5   14    1.86   121      0.15      0.03  9       18.07
             1.5   16    1.78   156      0.13      0.03  9       18.06
             1.5   18    1.77   185      0.12      0.02  9       18.13
             1.5   20    1.83   228      0.10      0.02  5       17.90
             3.0   14    3.14   280      0.12      0.02  8       16.02
             3.0   16    3.21   350      0.11      0.02  8       15.99
             3.0   18    3.25   440      0.07      0.01  8       15.92
             4.5   14    4.79   463      0.14      0.01  8       13.83 ];

%% Load NASA maybe you need it at some point?
global Runiv %#ok<GVMIS> (needed by the NASA functions)
Runiv = 8.314;
[SpS,El] = myload(fullfile(scriptDir,'Nasa','NasaThermalDatabase.mat'),{'Diesel','O2','N2','CO2','H2O'});

%% Engine geom data (check if these are correct)
Cyl.Bore             = 104*mm;
Cyl.Stroke           = 85*mm;
Cyl.CompressionRatio = 21.5;
Cyl.ConRod           = 136.5*mm;
Cyl.TDCangle         = 0;
Vd = pi/4*Cyl.Bore^2*Cyl.Stroke;   % displacement volume [m^3]
% -- Valve events
CaIVO = -355; % Intake valve opening
CaIVC = -135; % Intake valve closing
CaEVO = 149;  % Exhaust valve opening
CaEVC = -344; % Exhaust valve closing

%% 1. Find all data files
dataDir = '';
for c = 1:numel(dataDirCandidates)
    if exist(dataDirCandidates{c},'dir'), dataDir = dataDirCandidates{c}; break; end
end
if isempty(dataDir)
    error('Folder Actual_Data_set not found. Set dataDirCandidates at the top.');
end

fList = dir(fullfile(dataDir,'*fdaq.txt'));
if isempty(fList), error('No fdaq files in %s', dataDir); end

pts = struct([]);
for k = 1:numel(fList)
    fName = fList(k).name;
    % e.g. 20261123_0000009_group 3_3.0 bar_14soi_fdaq.txt
    %      20261123_0000005_group 3_1.5 bar_14fdaq.txt
    tokRun  = regexp(fName,'^\d+_(\d+)_','tokens','once');
    tokLoad = regexp(fName,'_(\d+\.\d+) bar','tokens','once');
    tokSOI  = regexp(fName,'_(\d+)(?:soi_)?fdaq','tokens','once');
    sMatch  = dir(fullfile(dataDir,['*_' tokRun{1} '_*sdaq.txt']));
    if numel(sMatch) ~= 1
        error('Expected 1 sdaq file for run %s, found %d', tokRun{1}, numel(sMatch));
    end
    pts(k).run   = str2double(tokRun{1});
    pts(k).load  = str2double(tokLoad{1});   % IMEP set point [bar]
    pts(k).SOI   = str2double(tokSOI{1});    % [deg bTDC]
    pts(k).fName = fName;
    pts(k).sName = sMatch(1).name;
end
[~,order] = sortrows([[pts.load]' [pts.SOI]']);
pts  = pts(order);
nPts = numel(pts);

fprintf('Found %d test points in %s:\n', nPts, dataDir);
for k = 1:nPts
    fprintf('  %d: run %d, %.1f bar, SOI %d\n', k, pts(k).run, pts(k).load, pts(k).SOI);
end

%% 2. Loop over all test points
R = struct([]);
for k = 1:nPts
    pt = pts(k);
    fprintf('\n===== Point %d/%d: %.1f bar, SOI %d deg bTDC =====\n', k, nPts, pt.load, pt.SOI);

    %% 2a. Load fast daq (txt -> cached .mat)
    FullName = fullfile(dataDir, pt.fName);
    MatName  = strrep(FullName, '.txt', '.mat');
    txtInfo = dir(FullName); matInfo = dir(MatName);
    if ~isempty(matInfo) && matInfo.datenum >= txtInfo.datenum
        load(MatName,'dataIn');
    else
        fid = fopen(FullName,'r');
        dataIn = fscanf(fid,'%f',[3 Inf]).';
        fclose(fid);
        save(MatName,'dataIn');
    end
    Nrows = size(dataIn,1);
    if mod(Nrows, NdatapointsperCycle) ~= 0
        error('%s: rows not a multiple of %d', pt.fName, NdatapointsperCycle);
    end
    Ncycles = Nrows/NdatapointsperCycle;
    Ca      = reshape(dataIn(:,1),[],Ncycles);
    p       = reshape(dataIn(:,2),[],Ncycles)*bara;
    injection_current = reshape(dataIn(:,3),[],Ncycles);   % raw injector signal

    % Slow daq (txt -> cached .mat)
    FullNameSlowdaq = fullfile(dataDir, pt.sName);
    MatNameSlowdaq  = strrep(FullNameSlowdaq, '.txt', '.mat');
    txtInfoS = dir(FullNameSlowdaq); matInfoS = dir(MatNameSlowdaq);
    if ~isempty(matInfoS) && matInfoS.datenum >= txtInfoS.datenum
        load(MatNameSlowdaq,'dataSlowdaq');
    else
        fid = fopen(FullNameSlowdaq,'r');
        dataSlowdaq = fscanf(fid,'%f',[5 Inf]).';
        fclose(fid);
        save(MatNameSlowdaq,'dataSlowdaq');
    end
    dataSlowdaq = dataSlowdaq(nSkip+1:end,:);   % drop settling samples

    fuelMassFlow_g_s     = dataSlowdaq(:,1);
    intakeTemperature_C  = dataSlowdaq(:,2);
    exhaustTemperature_C = dataSlowdaq(:,3);
    intakePressure_Pa    = dataSlowdaq(:,4)*bara;
    airMassFlow_g_s      = dataSlowdaq(:,5);
    fuelFlow_kg_s = fuelMassFlow_g_s*g;
    airFlow_kg_s  = airMassFlow_g_s*g;

    %% 2b. Indicated work (average of all cycles)
    V = CylinderVolume(Ca,Cyl);
    VClosed = [V; V(1,:)];
    pClosed = [p; p(1,:)];
    WiCycles = sum(0.5 .* (pClosed(1:end-1,:) + pClosed(2:end,:)) .* ...
        diff(VClosed,1,1), 1);    % [J]
    WiAverage = mean(WiCycles);
    fprintf('Average indicated work over %d cycles = %.2f J\n', Ncycles, WiAverage);

    %% 2c. Combustion characterisation
    CaSOI = -pt.SOI;
    Opt.CaIVC   = CaIVC;
    Opt.CaEVO   = CaEVO;
    Opt.CaSOI   = CaSOI;
    Opt.CaPeg   = -180;
    Opt.TIVC    = 330;
    Opt.FiltWin = 4;
    Opt.Gamma   = 'nasa';                    % or constant
    Opt.pIntake = mean(intakePressure_Pa);   % ambient (intake p not measured)
    Opt.PegMode = 'polytropic';              % 'polytropic' or 'ambient'
    Opt.CaPoly  = [-110 -50];
    Opt.kPoly   = 1.32;

    Cam = Ca(:,1); Vm = V(:,1);
    pMean = mean(p,2);
    CC = CombustionAnalysis(Cam, pMean, Vm, SpS, Opt);
    fprintf('Mean cycle: CA10 = %.1f, CA50 = %.1f, CA90 = %.1f deg ATDC\n', CC.CA10, CC.CA50, CC.CA90);
    fprintf('Burn duration = %.1f deg, ignition delay = %.1f deg, Qtotal = %.0f J\n', ...
        CC.BurnDuration, CC.IgnitionDelay, CC.Qtotal);
    fprintf('pmax = %.1f bar at %.1f deg ATDC\n', CC.pMax/bara, CC.CApMax);

    % Per cycle
    CA50cyc = nan(Ncycles,1); pMaxcyc = nan(Ncycles,1);
    for i = 1:Ncycles
        Ci = CombustionAnalysis(Ca(:,i), p(:,i), V(:,i), SpS, Opt);
        CA50cyc(i) = Ci.CA50;
        pMaxcyc(i) = Ci.pMax;
    end
    fprintf('CA50 per cycle: mean %.2f, std %.2f deg | pmax per cycle: mean %.1f, std %.2f bar\n', ...
        mean(CA50cyc), std(CA50cyc), mean(pMaxcyc)/bara, std(pMaxcyc)/bara);

    %% 2d. Diagnostics (one point only)
    if showDiagnostics && k == diagPoint %#ok<UNRCH>
        fprintf('\n--- DIAGNOSTICS (point %d) ---\n', k);
        fprintf('V(0)    = %.1f cm^3 (expect 35.2)\n', interp1(Cam,Vm,0)*1e6);
        fprintf('V(-180) = %.1f cm^3 (expect 757)\n',  interp1(Cam,Vm,-180)*1e6);

        iB = abs(Cam + 180) <= 3;
        fprintf('raw p at BDC = %.2f bar (expect about 1.0 absolute)\n', mean(pMean(iB))/bara);
        [~,ip] = max(pMean); fprintf('pmax at %.1f deg\n', Cam(ip));

        idx  = Cam >= -100 & Cam <= -30;
        pPeg = pMean - mean(pMean(iB)) + 1.0e5;
        c = polyfit(log(Vm(idx)), log(pMean(idx)), 1); n_raw = -c(1);
        c = polyfit(log(Vm(idx)), log(pPeg(idx)),  1); n_peg = -c(1);
        fprintf('n raw = %.2f | n pegged = %.2f (expect 1.25-1.35)\n', n_raw, n_peg);

        OptT = Opt;
        OptT.PegMode = 'none';       CC0 = CombustionAnalysis(Cam,pMean,Vm,SpS,OptT);
        OptT.PegMode = 'ambient';    OptT.pIntake = 1.0e5;
        CC1 = CombustionAnalysis(Cam,pMean,Vm,SpS,OptT);
        OptT.PegMode = 'polytropic'; CC2 = CombustionAnalysis(Cam,pMean,Vm,SpS,OptT);
        base = @(C) mean(C.aROHR(C.Ca >= -60 & C.Ca <= -20));
        fprintf('%-11s | aROHR0 | CA50 | CA90 | Qtot\n','peg');
        fprintf('%-11s | %6.2f | %4.1f | %4.1f | %4.0f\n','none',      base(CC0),CC0.CA50,CC0.CA90,CC0.Qtotal);
        fprintf('%-11s | %6.2f | %4.1f | %4.1f | %4.0f\n','ambient',   base(CC1),CC1.CA50,CC1.CA90,CC1.Qtotal);
        fprintf('%-11s | %6.2f | %4.1f | %4.1f | %4.0f\n','polytropic',base(CC2),CC2.CA50,CC2.CA90,CC2.Qtotal);
        fprintf('Rair = %g | gamma %.3f..%.3f\n', CC2.Rair, min(CC2.gamma), max(CC2.gamma));

        figure(100); loglog(Vm/dm^3, pPeg/bara,'b'); hold on
        i0 = find(Cam >= -180,1);
        loglog(Vm/dm^3, (Vm(i0)./Vm).^1.30, 'k--'); grid on; hold off
        xlabel('V [dm^3]'); ylabel('p [bar]'); legend('measured (pegged)','n = 1.30');
        title('Log-log pV: compare slope with the ideal line');

        inj = mean(injection_current,2); iS = find(inj > 0.5*max(inj),1);
        fprintf('SOI from injector signal = %.1f deg | Opt.CaSOI = %.1f deg\n', Cam(iS), Opt.CaSOI);

        fprintf('\n--- TDC shift test (ambient peg) ---\n');
        OptS = Opt; OptS.PegMode = 'ambient'; OptS.pIntake = 1.0e5;
        for sh = -2:1:2
            Vs = CylinderVolume(Cam + sh, Cyl);
            cc = polyfit(log(Vs(idx)), log(pPeg(idx)), 1);
            Cs = CombustionAnalysis(Cam,pMean,Vs,SpS,OptS);
            fprintf('shift %+d deg: n = %.2f | pre-ign aROHR %.2f | CA50 %.1f | CA90 %.1f | Qtot %.0f J\n', ...
                sh, -cc(1), mean(Cs.aROHR(Cs.Ca >= -60 & Cs.Ca <= -20)), Cs.CA50, Cs.CA90, Cs.Qtotal);
        end

        fprintf('Polytropic peg offset = %.3f bar | implied p at BDC = %.2f bar\n', ...
            CC.pegOffset/bara, (interp1(Cam,pMean,-180) + CC.pegOffset)/bara);
        for kk = [1.28 1.32 1.36]
            OptK = Opt; OptK.kPoly = kk;
            Ck = CombustionAnalysis(Cam,pMean,Vm,SpS,OptK);
            fprintf('k=%.2f: offset %.3f bar | CA50 %.1f | CA90 %.1f | Qtot %.0f J | pre-ign aROHR %.2f\n', ...
                kk, Ck.pegOffset/bara, Ck.CA50, Ck.CA90, Ck.Qtotal, ...
                mean(Ck.aROHR(Ck.Ca >= -60 & Ck.Ca <= -20)));
        end

        figure(101); set(gcf,'Position',[200 500 1200 400]);
        plot(Ca,p/bara,'LineWidth',1); hold on
        plot(Cam,pMean/bara,'k','LineWidth',2);
        YLIM = ylim;
        line([CaIVC CaIVC],YLIM,'Color','b'); line([CaEVO CaEVO],YLIM,'Color','r');
        line([CaSOI CaSOI],YLIM,'Color','y'); hold off
        xlabel('Ca'); ylabel('p [bar]'); xlim([-360 360]); grid on
        title(sprintf('All cycles, %.1f bar, SOI %d', pt.load, pt.SOI));
        fprintf('--- END DIAGNOSTICS ---\n');
    end

    %% 2e. KPIs  [NEW]
    % Only cylinder pressure is measured, no shaft torque, so these are
    % INDICATED values (is = indicated specific), not brake values.
    Pi_W     = WiAverage * RPM/60/2;     % indicated power [W] (4-stroke)
    IMEP_bar = WiAverage / Vd / bara;

    mf   = mean(fuelFlow_kg_s);          % [kg/s]
    ma   = mean(airFlow_kg_s);           % [kg/s]
    mexh = mf + ma;                      % exhaust mass flow [kg/s]
    nexh = mexh / (Mol.exh*g);           % exhaust molar flow [mol/s]
    lambda = (ma/mf) / 14.50;            % AFRst B7 = 14.50

    toSpec = @(mdot) mdot / Pi_W * kWhr / g;   % [kg/s] -> [g/kWh]

    % Emission values for this point
    row = gasTable(gasTable(:,1) == pt.load & gasTable(:,2) == pt.SOI, :);
    if isempty(row), warning('No gas data for %.1f bar / SOI %d', pt.load, pt.SOI); row = nan(1,8); end
    gas.CO2 = row(3); gas.NOx = row(4); gas.FSN = row(5);
    gas.CO  = row(6); gas.HC  = row(7); gas.O2  = row(8);

    % CO2, method a: fuel carbon balance (all C -> CO2, slide 43)
    mCO2_fuel = mf * Fuel.wC * Mol.CO2/Mol.C;
    % CO2, method b: analyser (vol fraction x molar exhaust flow)
    mCO2_meas = gas.CO2*volperc * nexh * Mol.CO2*g;
    % NOx from analyser
    mNOx = gas.NOx*ppm * nexh * Mol.NOx*g;

    % --- Main KPIs with the existing CalculateKPI function ---
    % NOTE: Pi_W is INDICATED power (no torque sensor), so the outputs
    % named "Brake..." / "BS..." are really indicated values here.
    KPI = CalculateKPI(Pi_W, mf, Fuel.LHV, mCO2_fuel, mNOx);
    eta   = KPI.BrakeEfficiency;   % indicated efficiency [-]
    isFC  = KPI.BSFC_g_kWh;        % [g/kWh]
    isCO2 = KPI.BSCO2_g_kWh;       % [g/kWh]
    isNOx = KPI.BSNOx_g_kWh;       % [g/kWh]

    % --- Extra emissions, not in CalculateKPI yet [NEW] ---
    % CO and HC from analyser
    mCO  = gas.CO*volperc * nexh * Mol.CO*g;
    mHC  = gas.HC*ppm     * nexh * Mol.HC*g;
    % Soot: FSN -> mg/m^3 (AVL correlation), times exhaust volume flow
    Csoot_mg_m3 = (1/0.405) * 4.95 * gas.FSN * exp(0.38*gas.FSN);
    mSoot = Csoot_mg_m3*1e-6 * (mexh/rhoExhStd);   % [kg/s]

    fprintf('IMEP %.2f bar | Pi %.0f W | eta %.1f %% | isFC %.0f | isCO2 %.0f (analyser %.0f) | isNOx %.2f | soot %.3f g/kWh\n', ...
        IMEP_bar, Pi_W, eta*100, isFC, isCO2, toSpec(mCO2_meas), isNOx, toSpec(mSoot));

    %% Store results of this point
    R(k).run       = pt.run;
    R(k).load      = pt.load;
    R(k).SOI       = pt.SOI;
    R(k).Ncycles   = Ncycles;
    R(k).IMEP_bar  = IMEP_bar;
    R(k).Pi_W      = Pi_W;
    R(k).lambda    = lambda;
    R(k).T_exh_C   = mean(exhaustTemperature_C);
    R(k).eta       = eta;
    R(k).isFC      = isFC;
    R(k).isCO2     = isCO2;
    R(k).isCO2_meas= toSpec(mCO2_meas);
    R(k).isNOx     = isNOx;
    R(k).isCO      = toSpec(mCO);
    R(k).isHC      = toSpec(mHC);
    R(k).isSoot    = toSpec(mSoot);
    R(k).CA10      = CC.CA10;
    R(k).CA50      = CC.CA50;
    R(k).CA90      = CC.CA90;
    R(k).CA50_std  = std(CA50cyc);
    R(k).pMax_bar  = CC.pMax/bara;
    R(k).CApMax    = CC.CApMax;
    R(k).Ca        = Cam;
    R(k).pMean_bar = pMean/bara;
    R(k).V_dm3     = Vm/dm^3;
    R(k).CC        = CC;
    R(k).CA50cyc   = CA50cyc;
end

save(fullfile(scriptDir,'Results_AllPoints.mat'),'R');

%% 3. Summary table + plots
fprintf('\n%5s %4s %6s %6s %6s %7s %7s %7s %6s %7s %6s\n', 'Load','SOI','IMEP','eta%', ...
    'isFC','isCO2','CO2meas','isNOx','isCO','isSoot','CA50');
for k = 1:nPts
    fprintf('%5.1f %4d %6.2f %6.1f %6.0f %7.0f %7.0f %7.2f %6.2f %7.3f %6.1f\n', R(k).load, R(k).SOI, ...
        R(k).IMEP_bar, R(k).eta*100, R(k).isFC, R(k).isCO2, R(k).isCO2_meas, ...
        R(k).isNOx, R(k).isCO, R(k).isSoot, R(k).CA50);
end

labels = arrayfun(@(r) sprintf('%.1f bar, SOI %d', r.load, r.SOI), R, 'UniformOutput', false);

figure(1)
hold on
for k = 1:nPts, plot(R(k).Ca, R(k).pMean_bar, 'LineWidth', 1); end
hold off
xlabel('Ca [deg]'); ylabel('p [bar]'); xlim([-60 90]); grid on
legend(labels,'Location','northeast'); title('Mean cylinder pressure, all points')

figure(2); set(gcf,'Position',[200 200 1100 450]);
subplot(1,2,1); hold on
for k = 1:nPts, plot(R(k).V_dm3, R(k).pMean_bar); end
hold off; grid on; xlabel('V [dm^3]'); ylabel('p [bar]'); xlim([0 0.8])
title('pV-diagram (mean cycle), linear')
subplot(1,2,2)
for k = 1:nPts, loglog(R(k).V_dm3, R(k).pMean_bar); hold on; end
hold off; grid on; xlabel('V [dm^3]'); ylabel('p [bar]'); xlim([0.03 0.8])
title('pV-diagram (mean cycle), log-log'); legend(labels,'Location','southwest')

figure(3)
hold on
for k = 1:nPts, plot(R(k).CC.Ca, R(k).CC.aROHR); end
hold off
xlabel('Ca [deg]'); ylabel('aROHR [J/deg]'); xlim([-30 90]); grid on
legend(labels); title('Apparent rate of heat release')

figure(4)
hold on
for k = 1:nPts, plot(R(k).CC.Ca, R(k).CC.xb); end
hold off
xlabel('Ca [deg]'); ylabel('Mass fraction burned [-]'); xlim([-30 90]); grid on
legend(labels,'Location','southeast'); title('Cumulative heat release')

figure(5)
hold on
for k = 1:nPts, plot(1:R(k).Ncycles, R(k).CA50cyc, '.'); end
hold off
xlabel('Cycle [-]'); ylabel('CA50 [deg ATDC]'); grid on
legend(labels); title('CA50 per cycle')

% KPIs vs SOI, one line per load
loads = unique([R.load]);
kpiList = {'eta','Indicated efficiency [-]'; 'isCO2','isCO2 [g/kWh]'; ...
           'isNOx','isNOx [g/kWh]'; 'isSoot','Soot [g/kWh]'; ...
           'CA50','CA50 [deg ATDC]'; 'isCO','isCO [g/kWh]'};
figure(6); set(gcf,'Position',[100 100 1100 650]);
for q = 1:size(kpiList,1)
    subplot(2,3,q); hold on
    for L = loads
        sel = [R.load] == L;
        plot([R(sel).SOI], [R(sel).(kpiList{q,1})], 'o-', 'LineWidth', 1.5);
    end
    hold off; grid on
    xlabel('SOI [deg bTDC]'); ylabel(kpiList{q,2});
end
subplot(2,3,1);
legend(arrayfun(@(L) sprintf('%.1f bar', L), loads, 'UniformOutput', false), 'Location','best');

%% 4. Optimal SOI  [TODO - teammate]
% Starting point: per load, SOI with the highest efficiency and the NOx
% at that SOI. Replace with the agreed criterion (e.g. CA50 target,
% NOx limit, trade-off weighting).
fprintf('\nOptimal SOI (max efficiency, first guess):\n');
for L = loads
    sel = find([R.load] == L);
    [~,ib] = max([R(sel).eta]);
    b = sel(ib);
    fprintf('  %.1f bar: SOI %d deg bTDC | eta %.1f %% | CA50 %.1f | isNOx %.2f g/kWh\n', ...
        L, R(b).SOI, R(b).eta*100, R(b).CA50, R(b).isNOx);
end