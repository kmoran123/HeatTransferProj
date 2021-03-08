clear
clc
format compact

% Transport 3 Hot Dog Project Code
% for Katy Moran and Anna Hartig


% list of material properties and design constraints
densityHotDog = 880;
kHotDog = 0.52;
cHotDog = 3350;
diameterHotDog = 0.0254;

initalTemp = 15 + 273.15;
maxTs = 100 + 273.15;
finalTemp = 60 + 273.15;
gasTemp = 295 + 273.15;
coalsTemp = 425 + 273.15;
emissivityCoals = 0.8;
emissivityHotDog = 0.45;
diameterCoals = 0.045;

deltar = 0.001; %**delta r** can be adjusted if needed
alphaHotDog = kHotDog/(densityHotDog*cHotDog);

% h rad calculation
hRad = 28.7;

% h convec calculation
hConv = 15;

%total h
hTot = hRad + hConv;

%stability criteria calculation for delta t
deltatSurface = (deltar^2 * diameterHotDog * kHotDog)/(alphaHotDog*(2*diameterHotDog*kHotDog - kHotDog + 2*hTot*diameterHotDog*deltar));
deltatInterior = deltar^2/(2*alphaHotDog);
deltat = 0;

if deltatInterior > deltatSurface 
    deltat = deltatInterior;
    
else
    deltat = deltatSurface;
end

    