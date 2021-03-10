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
radiusHotDog = diameterHotDog/2;

initialTemp = 15 + 273.15;
maxTs = 100 + 273.15;
finalTemp = 60 + 273.15;
gasTemp = 295 + 273.15;
coalsTemp = 425 + 273.15;
emissivityCoals = 0.8;
emissivityHotDog = 0.45;
diameterCoals = 0.045;

deltar = 0.0001; %**delta r** can be adjusted if needed
M = radiusHotDog/deltar;
alphaHotDog = kHotDog/(densityHotDog*cHotDog);

% h rad calculation
hRad = 0;

% h convec calculation
hConv = 0.00001;

%total h
hTot = hRad + hConv;

%stability criteria calculation for delta t
deltatSurface = (deltar^2 * M * kHotDog)/(alphaHotDog*(2*M*kHotDog - kHotDog + 2*hTot*M*deltar));
deltatInterior = deltar^2/(2*alphaHotDog);
deltatguess = 0.1; %preferred value of delta t that may or may not fit stability criteria

%here this is determining which delta t value to use. 0.03 is the default,
%but if the stability criteria not fulfilled, will use the other values
if (deltatguess > deltatSurface && deltatguess > deltatInterior) 
    deltat = deltatguess;
elseif deltatInterior > deltatSurface 
    deltat = deltatInterior;
else
    deltat = deltatSurface;
end

%initialize the temp matrix
temperature = zeros(1000, ceil(M)+1);
temperature(1,:) = initialTemp;

% numerical solution
% loop until centerline Temp = 60 degrees
i=1;
% while temperature(i,0) < finalTemp
%     %codey code code
%     i=i+1;
% end