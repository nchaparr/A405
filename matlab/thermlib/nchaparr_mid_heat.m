c=constants;
wtA=14.e-3;
pressA=900.e2;
tempA=25 + c.Tc;

TdA=findTdwv(wtA,pressA);
thetaeA=thetaep(TdA,tempA,pressA);
wtB=wtA;
pressB=700.e2;
TdB=findTdwv(wtB,pressB);
thetaeB=thetaeA;
[tempB,wvB,wlB]=tinvert_thetae(thetaeB, wtB, pressB);

wtC=wtA;
pressC=900.e2;
TdC=findTdwv(wtC,pressC);
tempC=tempB;
thetaeC=thetaep(TdC,tempC,pressC);
skew=30.;
figHandle=figure(1);
[figureHandle,outputws,handlews]=makeSkew(figHandle,skew);
xtempA=convertTempToSkew(tempA - c.Tc,pressA*0.01,skew);
xtempB=convertTempToSkew(tempB - c.Tc,pressB*0.01,skew);
xtempC=convertTempToSkew(tempC - c.Tc,pressC*0.01,skew);
text(xtempA,pressA*0.01,'A',...,
            'edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempB,pressB*0.01,'B',...
            'edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempC,pressC*0.01,'C',...
            'edgecolor','b','fontweight','bold','fontsize',22,'color','b');
pressLevs=linspace(700,900,60)*100.;
for i=1:numel(pressLevs)
    thePress=pressLevs(i);
    [temp,wv,wl]=tinvert_thetae(thetaeA, wtA, thePress);
    lineAB(i)=temp;
    rho=thePress/(c.Rd*temp);
    rhoAB(i)=rho;
end
tempCA=linspace(tempC,tempA,100);
rhoCA=pressA./(c.Rd*tempCA);
press900Vec=NaN(size(rhoCA));
press900Vec(:)=pressA;
rhoBC=pressLevs/(c.Rd*tempB);
xtemp=convertTempToSkew(lineAB - c.Tc,pressLevs*0.01,skew);    
semilogy(xtemp,pressLevs*0.01,'k-.','linewidth',2);
semilogy([xtempB,xtempC],[700.,900.],'b-.','linewidth',2);
semilogy([xtempC,xtempA],[900.,900.],'r-.','linewidth',2);
title('heat engine problem');
xleft=convertTempToSkew(10,1000.,skew);
xright=convertTempToSkew(30,1000.,skew);
axis([xleft,xright,650,1000.]);
print -depsc prob1.eps

figure(2)
clf;
plot(1./rhoCA,press900Vec*1.e-2);
ylim([700,1000.]);
hold on;
plot(1./rhoAB,pressLevs*1.e-2);
plot(1./rhoBC,pressLevs*1.e-2);
title('volume - temperature plot')
hold off;

%Calculate work by integrating over area between pressure vs alpha curves:

%Set up arrays of specific voloumes and corresponding pressures 
%for each transition, as represented on the plot:

alphaAB = 1./rhoAB;
pressAB = pressLevs;

alphaCA = 1./rhoCA;
pressCA = press900Vec;

alphaBC = 1./rhoBC;
pressBC = pressLevs;

%Calculate the area under constant pressure line:

disp(alphaCA(end) - alphaCA(1));
disp(pressA - pressB);
workCA = (alphaCA(end) - alphaCA(1))*(pressA - pressB);
disp('workCA');
disp(workCA);

%For the curves representing the other two transitions
%et up arrays of specific volumes, for which pressure values
%will be calculated using the interpolated functions:
alphas = linspace(alphaAB(1), alphaAB(end), 100);
alphas1 = linspace(alphaBC(1), alphaBC(end), 100);

%Create the interpolated functions:
alphaf1 = @(alphas) interp1(alphaAB, pressAB, alphas);
alphaf2 = @(alphas1) interp1(alphaBC, pressBC, alphas1);

%Combine them to get the difference in area between the two curves:
work = workCA - quad(alphaf1, alphaAB(1), alphaAB(end)) + quad(alphaf2, alphaBC(1), alphaBC(end));
disp('The three areas: C to A, A to B and B to C ')
disp(workCA); 
disp(quad(alphaf1, alphaAB(1), alphaAB(end)));
disp(quad(alphaf2, alphaBC(1), alphaBC(end)));
disp('work obtained by difference of pressure-alpha curve integrals');
disp(work);
alphaf2 = @(alphas) interp1(alphaBC, pressBC, alphas1);

%Tests
%Do rough estimate of integral under 1st curve
disp('triangle one');
disp(0.5*(alphas(end) - alphas(1))*20000);

%Do rough estimate of integral under 2nd curve
disp('triangle two');
disp(0.5*(alphas(end) - alphas(1))*(pressBC(end) - pressBC(1)));

%test against values calculated directly from formulae
[tempC,wvC,wlC]=tinvert_thetae(thetaeC, wtC, pressC);
wvA = wtA;
Qin = c.cpd*(tempA - tempC) + c.lv0*(wvA-wvC);
Qout = tempB*c.cpd*(log(thetaeC) - log(thetaeA));
disp('work calculated analytically');
disp(Qin + Qout);






