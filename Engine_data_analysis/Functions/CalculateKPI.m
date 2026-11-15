function KPI = CalculateKPI(brakePower_W, fuelFlow_kg_s, ...
    LHV_J_kg, CO2Flow_kg_s, NOxFlow_kg_s)
% Calculate KPIs for one engine operating point.
%
% Inputs:
% brakePower_W : mechanical shaft power [W]
% fuelFlow_kg_s: fuel mass flow rate [kg/s]
% LHV_J_kg     : lower heating value [J/kg]
% CO2Flow_kg_s : CO2 mass flow rate [kg/s]
% NOxFlow_kg_s : NOx mass flow rate [kg/s]

% Energy supplied by the fuel per second [W]
fuelPower_W = fuelFlow_kg_s * LHV_J_kg;

% Brake thermal efficiency [-]
brakeEfficiency = brakePower_W / fuelPower_W;

% Conversion from kg/J to g/kWh
conversionFactor = 1000 * 3.6e6;

% Brake-specific fuel consumption and emissions [g/kWh]
BSFC_g_kWh = fuelFlow_kg_s / brakePower_W * conversionFactor;
BSCO2_g_kWh = CO2Flow_kg_s / brakePower_W * conversionFactor;
BSNOx_g_kWh = NOxFlow_kg_s / brakePower_W * conversionFactor;

% Results table
KPI = table(brakeEfficiency, BSFC_g_kWh, ...
    BSCO2_g_kWh, BSNOx_g_kWh, ...
    'VariableNames', {'BrakeEfficiency', 'BSFC_g_kWh', ...
    'BSCO2_g_kWh', 'BSNOx_g_kWh'});

end