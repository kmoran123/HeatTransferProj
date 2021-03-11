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
hRad = 15;
%hello
% h convec calculation
hConv = 15;

%total h
hTot = hRad + hConv;


%stability criteria calculation for delta t
deltatSurface = (deltar^2 * M * kHotDog)/(alphaHotDog*(2*M*kHotDog - kHotDog + 2*hTot*M*deltar));
deltatCenter = (deltar^2)/(4*alphaHotDog);
deltatInterior = (deltar^2)/(2*alphaHotDog);

%here this is determining which delta t value to use. 0.03 is the default,
%but if the stability criteria not fulfilled, will use the other values
if (deltatInterior < deltatSurface && deltatInterior < deltatCenter) 
    deltat = deltatInterior;
elseif deltatCenter < deltatSurface 
    deltat = deltatCenter;
else
    deltat = deltatSurface;
    fprintf('Surface Stability')
end

%Calculating Fo and Bi
Fo = (deltat*alphaHotDog)/deltar^2;
%Bi = hTot*ro/kHotDog;

%initializing the matrices
temperature = zeros(16798, ceil(M)+1);
time=zeros(16798,1);
temperature(1,:) = initialTemp;

% numerical solution
% loop until centerline Temp = 60 degrees
i=1;

while temperature(i,1) < finalTemp
    for j=1:ceil(M)+1
        if j==1
            temperature(i+1,j) = 4*Fo*temperature(i, j+1) - (4*Fo-1)*temperature(i,j);
        elseif j==ceil(M)+1
            temperature(i+1, j) = (1 - 2*Fo*(1-(1/(2*M))) - 2*(hTot/kHotDog)*Fo*deltar)*temperature(i,j) + ...
                2*Fo*(1-(1/(2*M)))*temperature(i, j-1) + 2*(hTot/kHotDog)*Fo*deltar*gasTemp;
        else
            temperature(i+1, j) = (1-2*Fo)*temperature(i,j) + Fo*(1-(1/(2*(j-1))))*temperature(i, j-1) + ...
                Fo*(1+(1/(2*(j-1))))*temperature(i, j+1);
        end
            
    end
    time(i+1) = deltat*i;
    i=i+1;
    
end

xmatrix = zeros(128,1);
for k=1:M+1
    xmatrix(k+1) = k*deltar;
end

figure
plot(time, temperature(:,1), time, temperature(:,128));
title('Numerical Solution');
xlabel('time, [seconds]')
ylabel('Temperature, [K]')
legend('Centerline', 'Surface');

figure
plot(xmatrix, temperature(22,:), xmatrix, temperature(8468,:), xmatrix, temperature(16798,:));
title('Numerical Solution');
xlabel('Distribution, [cm]')
ylabel('Temperature, [K]')
legend('t = 0.2976 sec', 't = 120.0034 seconds', 't = 238.0652 seconds');


