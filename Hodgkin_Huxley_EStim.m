% Hodgkin_Huxley.m
%
% Matlab implementation of the Hodgkin_Huxley model
% Uses the simplest (lowest order) approximation of the differential
% equation, so it's easy to read, but falls appart for large time steps
%
% Written on 11/20/06 by G.M. Boynton at the Salk Institute, SNL-B
% Modified for CSHL 2014, 7/18/14
%Generate time-course for input current.


dt=stim.tsample;
t = linspace(0, length(Iin)*dt, length(Iin));

% % Parameters for voltage clamp
gClamp = 4;

%Equilibrium potentials
VNa = 115;      %sodium
VK = -12;       %potassium
Vl = 10.613;    %extracellular

C=1;            %membrane capacitance
gL = 0.3;       %leak conductance
T=6.3;          %temperature (C)
k=3^(.1*T-.63); %constant (function of T)

%zero out vectors
V = zeros(1,length(t));  %membrane potential
n = zeros(1,length(t));  %'n' gate
m = zeros(1,length(t));  %'m' gate
h = zeros(1,length(t));  %'h' gate
gNa = zeros(1,length(t));  %Sodium conductance
gK  = zeros(1,length(t));  %Potassium conductance

%nonzero initial parameters (important!)
m(1)=0.0505;
n(1)=0.3835;
h(1)=0.4782;

for i=1:length(t)-1
    %voltage dependent alpha values for m,n and h
    am = (2.5-.1*V(i))/(exp(2.5-.1*V(i))-1);
    an = (1-.1*V(i))/(10*(exp(1-.1*V(i))-1));
    ah = .07*exp(-V(i)/20);
    
    %voltage dependent beta values for m,n and h
    bm = 4*exp(-V(i)/18);
    bn = 0.125*exp(-V(i)/80);
    bh = 1/(exp(3-.1*V(i))+1);
    
    %differential equations for m,n and h
    m(i+1) = m(i)+ dt*k*(am*(1-m(i))-bm*m(i));
    n(i+1) = n(i)+ dt*k*(an*(1-n(i))-bn*n(i));
    h(i+1) = h(i)+ dt*k*(ah*(1-h(i))-bh*h(i));
    
    %sodium and potassium conductances (functions of m,n and h)
    gNa(i+1) = 120*m(i+1)^3*h(i+1);
    gK(i+1) = 36*n(i+1)^4;
    
    %Differential equation for membrane potential
    V(i+1)=V(i)+dt/C*(Iin(i)-gK(i+1)*(V(i)-VK)-gNa(i+1)*(V(i)-VNa)-gL*(V(i)-Vl));
end

%The rest is stuff for plotting.
fontSize = 10;
figure(11)  %Membrane potential over time
clf
plot(t,V,'LineWidth',2);
xlabel('Time (ms)','FontSize',fontSize);
ylabel('Voltage (mv)','FontSize',fontSize);
set(gca,'YLim',[-20,150]);
title(sprintf('%5.2g mV injected',max(Iin)),'FontSize',fontSize);
hold on
%plot(t, Iin,'r-','LineWidth',.5);
set(gcf,'PaperPosition',[1,1,11,2]);
set(gca,'FontSize',fontSize');

figure(12) %Na and K conductance over time
clf
hold on
plot(t,gNa,'r-');
plot(t,gK,'g-');
legend({'gNa','gK'},'FontSize',fontSize);
xlabel('Time (ms)','FontSize',fontSize);
ylabel('Conductance','FontSize',fontSize);
set(gca,'FontSize',fontSize');

figure(13)  %m,n and h gates over time
clf
hold on
plot(t,m,'r-');
plot(t,n,'g-');
plot(t,h,'b-');
legend({'m','n','h'},'FontSize',fontSize);
xlabel('Time (ms)','FontSize',fontSize);
ylabel('p','FontSize',fontSize);
set(gca,'FontSize',fontSize');

figure(14)  %phase plots for m,n and h vs membrane potential
clf
subplot(1,3,1)
plot(m,V);
xlabel('m','FontSize',fontSize);
ylabel('V (mV)','FontSize',fontSize);
set(gca,'FontSize',fontSize');

subplot(1,3,2)
plot(n,V);
xlabel('n','FontSize',fontSize);
set(gca,'FontSize',fontSize');

subplot(1,3,3)
plot(h,V);
xlabel('h','FontSize',fontSize);
set(gca,'FontSize',fontSize');

figure(15)
%voltage dependent alpha values for m,n and h
v=linspace(min(V), max(V), 100);
am = (2.5-.1*v)./(exp(2.5-.1.*v)-1);
an = (1-.1*v)./(10*(exp(1-.1.*v)-1));
ah = .07*exp(-v./20);

%voltage dependent beta values for m,n and h
bm = 4*exp(-v./18);
bn = 0.125*exp(-v./80);
bh = 1./(exp(3-.1*v)+1);

list={'am', 'an', 'ah', 'bm', 'bn', 'bh'};
clist=['r', 'g', 'b', 'r', 'g', 'b'];
for i=1:6
    subplot(2, 3, i)
    plot(v, eval(list{i}), clist(i)); hold on
    xlabel('Voltage');
    ylabel(list{i})
    title(list{i});
    set(gca, 'XLim', [min(v) max(v)]);
end

%% plotting over time
return
figure(16); clf
for i=1:100:length(t)
    subplot(2,4 ,1)
    plot(t(1:i), n(1:i), 'g-'); title('n'); xlabel('time'); ylabel('n')
    set(gca, 'XLim', [0 max(t)]);
    set(gca, 'YLim', [0 1]);
    subplot(2,4 ,5); hold on
    plot(t(1:i), m(1:i), 'r-'); title('m'); xlabel('time'); ylabel('m / h')
    plot(t(1:i), h(1:i), 'b-');
    set(gca, 'XLim', [0 max(t)]);
    set(gca, 'YLim', [0 1]); hold on
    
    
    subplot(2, 4, 2)
    plot(t(1:i),gK(1:i),'g-'); title('gk'); xlabel('time');  ylabel('gK');
    set(gca, 'XLim', [0 max(t)]);
    set(gca, 'YLim', [-1 30]);
    subplot(2, 4, 6)
    plot(t(1:i),gNa(1:i),'m-'); title('gNa'); xlabel('time');  ylabel('gNa');
    set(gca, 'XLim', [0 max(t)]);
    set(gca, 'YLim', [-1 30]);
    
    subplot(2,4 ,3)
    plot(t(1:i), V(1:i), 'k-'); title('voltage'); xlabel('time');  ylabel('V')
    set(gca, 'XLim', [0 max(t)]);
    set(gca, 'YLim', [-25 150]);
    
    subplot(2, 4, 4)
    plot(n(1:i), V(1:i), 'g-'); title('n vs. V'); xlabel('n'); ylabel('V')
    set(gca, 'XLim', [.3 .8]);
    set(gca, 'YLim', [-25 150]);
    
    subplot(2, 4, 8)
    plot(m(1:i), V(1:i), 'r-'); title('m & h vs. V'); xlabel('m/h'); ylabel('V')
    plot(h(1:i), V(1:i), 'r-');
    set(gca, 'XLim', [0 1]);
    set(gca, 'YLim', [-25 150]);
    
    drawnow;
end

