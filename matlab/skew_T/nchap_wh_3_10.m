% -*- coding: utf-8 -*-
% $$$ Exercise 3.10 A parcel of air with an initial temperature
% $$$ of 14 °C and dew point 2 °C is lifted adiabatically
% $$$ from the 1000-hPa level. Determine its LCL
% $$$ and temperature at that level. If the air parcel is
% $$$ lifted a further 200 hPa above its LCL, what is its
% $$$ final temperature and how much liquid water is condensed
% $$$ during this rise?

c=constants;
press0 = 95.e3;
Temp0 = 14 + c.Tc;
wv0 = .008; %g/kg
wt0 = .008;
refpress = 1.e5;

% if unstaturated follow dry adiabat to LCl, then follow moist adiabat at
% LCL to refpress
% if saturated follow moist adiabat to refpress
% return wet bulb potential temperature by using rootfinder

Tdew = findTdwv(wv0, press0);
thetae = thetaep(Tdew, Temp0, press0);
thetaw = findTmoist(thetae, refpress);

fprintf('found Wet Bulb Potential Temperature=%8.2f (deg C)\n',...
        thetaw - c.Tc);
    
%raise to 700 hPa and return the new total water mixing ratio
%broken into vapor and liquid ratios, as well as the temperature at this
%level.  The new water vapor ratio is wsat at this Temperature and Pressure
%level.

[Temp1, wv1, wl1]=tinvert_thetae(thetae,wv0,7.e4);

fprintf('found temperature =%8.2f (deg C) at 700 \n',...
        Temp1 - c.Tc);
fprintf('found new water mixing ratios at 700 Pa=%8.2f, %8.2f (g/kg)\n',...
        wv1*1000, wl1*1000);
% after precipitation, there is .3 of the liquid water left in the air
% parcel.

wt1 = wv1 + .3*wl1;

fprintf('After precipitation the new water mixing ratio at 700 Pa=%8.2f (g/kg)\n',...
        wt1*1000);
    
%now bring it down to original pressure level noting that wsat, Td and the
%LCL are now different but the adiabat(s) and wet bulb potential temperature
%should be the same

press1 = 7.e4;
Tdew1 = findTdwv(wv1, press1);
thetae1 = thetaep(Tdew1, Temp1, press1);
thetaw1 = findTmoist(thetae1, refpress);
[Temp2,wv2,wl2]=tinvert_thetae(thetae, wt1, press0);

fprintf('New temperature=%8.2f (deg C)\n',...
        Temp2 - c.Tc);

fprintf('potential temperature=%8.2f (deg C)\n',...
        thetae1-c.Tc);
fprintf('mixing ratio=%8.2f (g/kg)\n',...
        (wv2 + wl2)*1000);
fprintf('wet bulb pot temp= %8.2f (deg C)\n',...
        thetaw1-c.Tc);
%8.2f (deg C)