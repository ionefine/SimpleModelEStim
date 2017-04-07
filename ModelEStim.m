%   =======================================================================
%   RetinalModel.m
%   Model designed to replicate the basic principles of estim
%
%   Date: 04.02.2017
%   =======================================================================
%   Some useful numbers
% cap=1-2.5 microF/cm2 capacity of cell membrane
% r_m resistance of resting cell membrane, ~1 kOhms
% rho_i - resistivity of axoplasm, 50-200 Ohms.cm
% rho_e - resistivity of extracellular fluid, 200 Ohms.cm

clear all; close all
figure(1); clf; set(gcf, 'Name', 'Time')
figure(2); clf; set(gcf, 'Name', 'Space')

stim.tsample=0.005; % sampling in ms
%% stimulus pulse train
stim.pulsedur=.1; %pulseduration in ms
stim.freq=50; % Hz
stim.dur=200; % duration of the train in ms
stim.amp=.6 / 1000; % stimulus current in  amps
stim.tsform=Istim(stim); % create pulse train

stim.t=0:stim.tsample:stim.dur;
stim.t=stim.t(1:end-1);
figure(1); clf; subplot(2, 2, 1); plot(stim.t, stim.tsform.*stim.amp.*1000); xlabel('time (ms)'); ylabel('uAmps'); 

% spatial position of the electrode
e.radius=5; % radius in microns
e.box=200; % length and width of box in which electrode is simulated
e.ssample=.5; % distance sampling in microns
[e.xgrid,e.ygrid] = meshgrid(-e.box:e.ssample:e.box);  %microns
e.rgrid =sqrt(e.xgrid.^2+e.ygrid.^2);

e.currentdensity=ones(size(e.rgrid)).*stim.amp./(pi*((e.radius/1000).^2));
e.currentdensity(e.rgrid>e.radius)=0;

figure(2)
subplot(2,2, 1); title('current density'); colormap(gray(256));
imagesc(e.currentdensity)
cb=colorbar; cb.Label.String=('Current Density (Amps/cm2)');

%% now model electrical response at a particular cell
% we're assuming this is the same tissue as the tissue directly under the
% electrode
cl.xloc=0 ;  cl.yloc=0 ; cl.zloc=40/1000; % cell distance in cm
cl.rho=300; %resistance extracellular medium, 0.3kOhm.cm is a reasonable number
% Frankenhauser and Huxley model for current/voltage drop off
s=0;
for x=-1:1
    for y=-1:1
        xloc=cl.xloc*1000/e.ssample;yloc=cl.yloc*1000/e.ssample;zloc=cl.zloc*1000/e.ssample;
        dgrid=sqrt((zloc.^2) + ((e.xgrid-xloc+x).^2) + ((e.ygrid-yloc+y).^2)); % distance from cell 2 the electrode array, in pix
        cl.dgrid=(dgrid*e.ssample)/1000; % convert to cm
        Vfield=(cl.rho./(4*pi*cl.dgrid))*(e.ssample.^2); % for this cell, this is the amount of current attentuation for each point in the electrode array
        % now add up the currents
        V=e.currentdensity.*Vfield;
        cl.V=sum(sum(V));
        s=s+cl.V;
        if x==0 && y==0
            c=cl.V;
            figure(2); subplot(2,2, 2);colormap(gray(256));
            imagesc(cl.dgrid);cb=colorbar; cb.Label.String=('Distance (cm)');
            title('distance 2 e @ cell');
            subplot(2, 2,3); 
            imagesc(Vfield*1000);cb=colorbar; cb.Label.String=('V Att');
            title('Att @ cell'); colormap(gray(256));
            subplot(2, 2,4); imagesc(V*1000); cb=colorbar; cb.Label.String=('Voltage&cell (mV)');
            title('V @ cell');      
        end
    end
end


Iin=(s-[9*c]).*stim.tsform; % convert to milivolts
Hodgkin_Huxley;

return
%% a digression for heating, not correct right now so IGNORE
% Joules law is the simple equation
% H=(I.^2) x resistivity x duration

% but a more accurate way of describing heating is by using the Pennes Bioheat equation. This includes
% the initial body temperature, and heat removal due to blood
% so now we need to define some properties of the tissue being stimulated
% directly under the electrode
figure(3); clf; set(gcf,'Name', 'Heat');
cl.TB=37; % initial blood temperature
cl.pB=1060; % mass density of blood, kh/m3
cl.cB=3960; % specific heat of blood, J/kg.K
cl.WB=0.0085; % blood perfusion rate, s-1

cl.limit=44; % safe limit for thermal damage
a=(max(e.currentdensity(:))).^2; % amplitude of stimulation squared
JL=cumsum(a.*(stim.tsample/1000).*abs(stim.tsform)).*(cl.rho*(e.ssample/1000)); % Joules Law heating
cl.heat=cl.TB+(JL./(cl.pB.*cl.cB.*cl.WB)); % Pennes Bioheat
figure(1); subplot(2,2,2); plot(stim.t, cl.heat);
title('heat, time'); xlabel('time'); ylabel('Temperature');

%% more of this digression, showing that Heat rises with the square of current level 
ct=1;
cdensities=0:.1:3;
for i=cdensities
    a=i.^2; % amplitude of stimulation
    JL=cumsum(abs(stim.tsform*(stim.tsample/1000)).*a).*(cl.rho*(e.ssample/1000)); % Joules Law heating
    heat(ct)=max(cl.TB+(JL./(cl.pB.*cl.cB.*cl.WB)));
    ct=ct+1;
end
figure(3);clf
subplot(2,2,1); plot(cdensities, heat); hold on
plot(max(e.currentdensity(:)), max(cl.heat), 'r*')
title('heat, current'); xlabel('current density'); ylabel('Temperature');

