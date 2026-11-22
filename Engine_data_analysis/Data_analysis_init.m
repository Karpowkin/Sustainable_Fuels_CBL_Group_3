%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
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
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355; % Intake valve opening
CaIVC = -135; % Intake valve closing
CaEVO = 149; % Exhaust valve opening
CaEVC = -344; % Exhaust valve closing
CaSOI = -3.2; % Start of Injection
% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams
%% Load data (if txt file)  
FullName = fullfile('Data','ExampleDataSet.txt');
MatName  = fullfile('Data','ExampleDataSet.mat');

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
p               = reshape(dataIn(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
p_injection     = reshape(dataIn(:,3),[],Ncycles)*bara;
%% Plotting 
f1=figure(1);
set(f1,'Position',[200 800 1200 400]);
plot(Ca,p/bara,'LineWidth',1);
hold on
plot(Ca,p_injection/bara,'LineWidth',1);
xlabel('Ca');ylabel('p [bar]');
xlim([-360 360]);ylim([0 50]);
iselect = 10;
plot(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
% YLIM = ylim;
% line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b');
% line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r');
% line([CaSOI CaSOI],YLIM,'LineWidth',1,'Color','y');
set(gca,'XTick',-360:60:360,'XGrid','on','YGrid','on');
title('All cycles in one plot.')
hold off
%% pV-diagram
V = CylinderVolume(Ca,Cyl);
Vselect = V(:,iselect);
f2 = figure(2);
set(f2,'Position',[200 400 600 800]);
subplot(2,1,1)
plot(Vselect/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');
xlim([0 0.8]);ylim([0.5 50]);
set(gca,'XTick',0:0.1:0.8,'XGrid','on','YGrid','on');
title({'pV-diagram',"linear scale"})
subplot(2,1,2)
loglog(Vselect/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');
xlim([0.02 0.8]);ylim([0 50]);
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');
title({'pV-diagram',"log-log scale"})
drawnow;

% Engine KPI

% Replace NaN when data will be obtained 
brakePower_W = NaN;   % Mechanical shaft power [W]
fuelFlow_kg_s = NaN;  % Fuel mass flow rate [kg/s]
LHV_J_kg = NaN;       % Fuel lower heating value [J/kg]
CO2Flow_kg_s = NaN;   % CO2 mass flow rate [kg/s]
NOxFlow_kg_s = NaN;   % NOx mass flow rate [kg/s]

KPI = CalculateKPI(brakePower_W, fuelFlow_kg_s, ...
    LHV_J_kg, CO2Flow_kg_s, NOxFlow_kg_s);

disp(KPI);
% Indicated work (average of all cycles)
VClosed = [V; V(1,:)];
pClosed = [p; p(1,:)];

WiCycles = sum(0.5 .* (pClosed(1:end-1,:) + pClosed(2:end,:)) .* ...
    diff(VClosed,1,1), 1);    % [J]
WiAverage = mean(WiCycles);

fprintf('Average indicated work over %d cycles = %.2f J\n', ...
    Ncycles, WiAverage);
% Maximum pressure across all cycles
[maxPressure_Pa, maxIdx] = max(p(:));
[maxCaIdx, maxCycle] = ind2sub(size(p), maxIdx);

fprintf('Maximum pressure across all cycles = %.2f bar (cycle %d at %.1f CA)\n', ...
    maxPressure_Pa/bara, maxCycle, Ca(maxCaIdx, maxCycle));
