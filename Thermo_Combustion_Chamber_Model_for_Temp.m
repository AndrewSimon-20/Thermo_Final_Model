% SFRJ Fuel consumption for uneven cylindrical fuel grain.

% Assumptions Made
% 1) Constant radial burn rate of Fuel at all points along fuel grain.
% 2) Insulation does not react when exposed. Or does so with insignificant gas
% and thermal generation.
% 3) Assumes slope of fuel grain from bow to aft to be linear and constant.
% 4) Assumes uniform fuel composition.

% Parameters  The parmeters that say 'don't touch this please' have been
% input from the real test data.

% Radial Burn rate of Fuel (regression rate) (meters/sec)
br = 0.00023368; % (meters/sec) % Don't Touch This Please.
% Thickness at thin side(meters)
FGThin = (0.00635 + 0.007874)/2; % (meters) % Don't Touch This Please 
% Thickness at thick side (meters)
FGThick = (0.007366 + 0.009652)/2; % (meters) %  Don't Touch This Please
% Outter diameter of fuel grain (meters)
oDFG = 0.066675; % (meters) % Don't Touch This Please
% Length of Fuel grain
lenFG = 0.508; % (meters) % Don't Touch This Please
% Density of Fuel
denFuel = 941; % (kg/m^3) % Don't Touch This Please
% Mass air In (kg/sec)
massAirIn = 0.5619; % (kg/sec) %Don't touch this please.

% Mols of HTPB per kg of Fuel
% Note this is 1000g/kg * wt % of HTPB in fuel / molar mass of "one unit of
% HTPB"
molHTPBpKG = (1000)*(0.66)/(125.188); % (mols/kg)
% Mols Styrene per kg of Fuel
molStyrenepKG = (1000)*(0.25)/(104.15); % (mols/kg)
% mols O2 consumed per kg of fuel
molO2pKG = 10*molStyrenepKG + 10.25*molHTPBpKG; % (mols/kg)
% mols CO2 produced per kg of fuel burned
molCO2pKG = 8*molStyrenepKG + 8*molHTPBpKG; % (mols/kg)
% mols H2O produced per kg of fuel burned
molH2OpKG = 4*molStyrenepKG + 6.5*molHTPBpKG; % (mols/kg)
% the initial temperature of the air entering the combustion chamber
T0 = 667 ;% (kelvin)
% Time interval of the burn from our data set.
tInitial = 0; % (seconds) % Please Don't Touch
tFinal = 3.8; % (seconds)  % Please Don't Touch



% Accuracy Modifiers

% Time interval between calculation of layers. (Smaller = More accurate).
tLayer = 0.01;
% dz (chunks of fuel being burned dz tall) (Smaller = More accurate) (meters).
dz = 0.001;





% Calculations

% Determines the slope of the fuel grain
slopeFG = abs((FGThin - FGThick)/lenFG) ;

% Breaks the possible situations into 2,fuel grain has one side that is thicker
% and Fuel grain is the same thickness throughout.
if FGThin <  FGThick
    % Calculates the max amount of time that there will still be fuel
    % burning.
    maxBurnTime = FGThick/br ;
    disp("This is the max burn Time:")
    disp(maxBurnTime)
    % Creates a vector for the inner diameter of the fuel at each point dz
    % along the length of the fuel grain. Uses an average fuel thickness
    % between the thicknesses on either side of dz.
    iDFG = zeros(1,ceil(lenFG/dz));
    for j = 1:ceil(lenFG/dz)
        iDFG(1,j) = oDFG - 2*FGThin - 2*(j-1)*dz*slopeFG;
    end
    % Creates a vector for the volume of fuel burned at each point dz at
    % each time point dt (dt = tLayer).
    volFuelBurned = zeros(ceil(maxBurnTime/tLayer),ceil(lenFG/dz));
    volFuelTotal = zeros(1,ceil(lenFG/dz));
    for j = 1:ceil(lenFG/dz)
        volFuelTotal(1,j) = ((pi/4)*(oDFG^2 - (iDFG(1,j))^2))*dz;
    end
    totalVolFuel = sum(volFuelTotal);
    
    for i = 1:ceil(maxBurnTime/tLayer)
        for k = 1:ceil(lenFG/dz)
            if iDFG(1,k)< oDFG
                volFuelBurned(i,k) = dz*(pi/4)*(((iDFG(1,k)+(br*tLayer))^2)-(iDFG(1,k)^2));
                iDFG(1,k) = iDFG(1,k)+ 2*br*tLayer;
            end
        end
    end
    disp("The volume of fuel burned in this period is:")
    if tInitial == 0
        disp(sum(sum(volFuelBurned(ceil(tInitial/tLayer)+1:ceil(tFinal/tLayer),:))))
    else
        disp(sum(sum(volFuelBurned(ceil(tInitial/tLayer):ceil(tFinal/tLayer),:))))
    end
    
elseif FGThin == FGThick
    % Calculates the max amount of time that there will still be fuel
    % burning.
    maxBurnTime = FGThin/br ;
    disp("This is the max burn Time:")
    disp(maxBurnTime)
    % Creates vector for volume of Fuel burned at each tLayer interval.
    volBurned = zeros(1,ceil(maxBurnTime/tLayer));
    % Creates the inner diameter of the fuel grain
    iDFG = oDFG - 2*FGThin;
    % Calculates volume of fuel at t = 0
    volFuelTotal = ((pi/4)*(oDFG^2 - (oDFG-FGThin)^2))*lenFG;
    for i = 1:ceil(maxBurnTime/tLayer)
        volBurned(1,i) = lenFG*(pi/4)*(((iDFG+(br*tLayer))^2)-(iDFG^2));
        iDFG = iDFG + 2*br*tLayer;
    end
    disp("The volume of fuel burned in this period is:")
    if tInitial == 0
        disp(sum(volBurned(1,ceil(tInitial/tLayer)+1:ceil(tFinal/tLayer))))
    else
        disp(sum(volBurned(1,ceil(tInitial/tLayer):ceil(tFinal/tLayer))))
    end
    
end




% Calculates mols of components consumed and produced.
if FGThin <  FGThick
    molsHTPBConsumed = zeros(1,ceil(maxBurnTime/tLayer));% Initializes a vector for mols of HTPB consumed at each time interval.
    for i = 1:ceil(maxBurnTime/tLayer)
        molsHTPBConsumed(1,i) = molHTPBpKG*sum(volFuelBurned(i,:))*denFuel;
    end
elseif FGThin ==  FGThick
    molsHTPBConsumed = zeros(1,ceil(maxBurnTime/tLayer));% Initializes a vector for mols of HTPB consumed at each time interval.
    for i = 1:ceil(maxBurnTime/tLayer)
        molsHTPBConsumed(1,i) = molHTPBpKG*volBurned(1,i)*denFuel;
    end
else
    disp("Please input the thicker end parmeter in for FGThick not for FGThin")
end

molsSTYRENEConsumed = molsHTPBConsumed*(molStyrenepKG/molHTPBpKG);
molsO2Consumed = molsHTPBConsumed*(molO2pKG/molHTPBpKG);
molsCO2Produced = molsHTPBConsumed*(molCO2pKG/molHTPBpKG);
molsH2OProduced = molsHTPBConsumed*(molH2OpKG/molHTPBpKG);

massO2Consumed = molsO2Consumed * 0.0319988; % (kg)
massH2OProduced = molsH2OProduced * 0.01801528;  % (kg)
massCO2Produced = molsCO2Produced * 0.04401; % (kg)

massFuelBurned = zeros(1,ceil(tFinal/tLayer));
massExit= zeros(1,ceil(tFinal/tLayer));
for ii = 1:ceil(tFinal/tLayer)
    massFuelBurned(1,ii) = sum(volFuelBurned(ii,:))*denFuel; % (kg/sec)
    massExit(1,ii) = massAirIn + massFuelBurned(1,ii); % (kg/sec)
end

massO2In = 0.2314*massAirIn; % (kg/sec)
massARIn = 0.0129*massAirIn;% (kg/sec)
massN2In = 0.7552*massAirIn; % (kg/sec)

massPctO2 = ((massO2In*tLayer)-massO2Consumed)/(massAirIn*tLayer); % time intervals = tLayer
massPctCO2 = (massCO2Produced)/(massAirIn*tLayer);  % time intervals = tLayer
massPctH2O = massH2OProduced/(massAirIn*tLayer); % time intervals = tLayer
massPctN2 = (massN2In*tLayer)/(massAirIn*tLayer); % time intervals = tLayer
massPctAR = (massARIn*tLayer)/(massAirIn*tLayer); % time intervals = tLayer
% Explanation of the + 0.0718: This gives < 1, the missing mass percent is from IPDI which we assume
% doesnt react. We are assuming that it has the same heat capacity of
% nitrogen so we will add 0.0718 to massPctN2 to account for it roughly.


tmid = ceil((tFinal/2)/tLayer) ; % Finds the time at the middle of the burn and calculates the mass percents at that time.
mPctO2 = massPctO2(tmid);
mPctCO2 = massPctCO2(tmid);
mPctH2O = massPctH2O(tmid);
mPctAr = massPctAR;
mPctN2 = massPctN2;
% NOTE these do not add to 100%, or 1, because we made the assumption mass
% flow rate in is equivilent to mass flow rate out because so little mass
% is burned. 

% stores the shomate equations for our components as functions
cpO2 = @(T) (30.03235 + 8.772972*T - 3.988133*(T^2) + 0.0788313*(T^3) - 0.741599/(T^2))/(31.998/1000); % range 700-2000K
cpCO2 = @(T) (58.16639 + 2.720074*T - 0.0492289*(T^2) + 0.038844*(T^3) - 6.447293/(T^2))/(44.01/1000); % range 1200-6000K
cpH2O = @(T) (30.092 + 6.832514*T + 6.793435*(T^2) - 2.53448*(T^3) + 0.082139/(T^2))/(18.01528/1000); % range 500 -1700k
cpAr = 0.52*1000 ; % the heat capacity of argon is constant out to the 4th sigfig from 298k to 6000k
cpN2 = @(T) (19.50583 + 19.88705*T - 8.598535*(T^2) + 1.369784*(T^3) + 0.527601/(T^2))/(28.0134/1000); % range 500 -2000k

% creates the weighted cp function
cpWtd= @(T) ((cpO2(T)*mPctO2)+(cpCO2(T)*mPctCO2)+(cpH2O(T)*mPctH2O)+(cpAr*mPctAr)+(cpN2(T)*mPctN2))/1000; % range 1200-1700k
% the division by 1000 above is to convert to KJ from Joules.
fuelMassFlowRate = sum(massFuelBurned(tmid-ceil(0.5/tLayer):tmid+ceil(0.5/tLayer))); % kg/sec

% heat of combustion of one kg of polystyrene is 42632 kj/kg
% heat of combustion of one kg of butadiene 47020 kj/kg
% Calculates the weighted heat of combustion
heatComb = 0.66*47020 + 0.25*42632; % KJ/KG fuel

bestT = 0 ;
error = 1;
for ij = 1.2:0.001:1.7
    if  error > (massAirIn*(ij*1000-T0)*cpWtd(ij))/(fuelMassFlowRate*heatComb)-1 && (massAirIn*(ij*1000-T0)*cpWtd(ij))/(fuelMassFlowRate*heatComb)-1 > 0  
        error = (massAirIn*(ij*1000-T0)*cpWtd(ij))/(fuelMassFlowRate*heatComb)-1;
        bestT = ij*1000;
    end
end


% The actual combustion chamber temp was 2292 Rankine or 1273 kelvin

clc
prctError = 100*abs(bestT-1273)/1273;
disp('The percent error of this calculation is:')
disp(prctError)
disp(' %')
