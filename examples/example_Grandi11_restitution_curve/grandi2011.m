% Human Atrial Action Potential and Ca2+ Model: Sinus Rhythm and Chronic Atrial Fibrillation.
% Grandi E, Pandit SV, Voigt N, Workman AJ, Dobrev D, Jalife J, Bers DM.
% Circ Res. 2011 Oct 14;109(9):1055-66. Epub 2011 Sep 15.
% PMID:21921263

%11/11/09 added ik1
%11/12/09 changed Itof, remove Itos
%11/12/09 add ikur
%11/18/09 corrected equations for Ikur as in Maleckar et al. 2009
%11/18/09 changed ICal - data from Li & Nattel AJP 1997
%12/04/09 SR Vmax *1.5
%12/09/09 IbarNaK = 0.7*
%12/09/09 INCX = 0.8*

function output = human_aclamp()
%
%SVP 11/11/09
clear all;
close all;

p = 0;%par;  % Parameter array for passing nondefault conditions
%% Initial conditions
mo=1.405627e-3;
ho= 9.867005e-1;
jo=9.915620e-1; 
do=7.175662e-6; 
fo=1.000681; 
fcaBjo=2.421991e-2;
fcaBslo=1.452605e-2;
xtoso=4.051574e-3;
ytoso=9.945511e-1; 
xtofo=4.051574e-3; 
ytofo= 9.945511e-1; 
xkro=8.641386e-3; 
xkso= 5.412034e-3;
RyRro=8.884332e-1;
RyRoo=8.156628e-7; 
RyRio=1.024274e-7; 
NaBjo=3.539892;
NaBslo=7.720854e-1; 
TnCLo=8.773191e-3; 
TnCHco=1.078283e-1; 
TnCHmo=1.524002e-2; 
CaMo=2.911916e-4; 
Myoco=1.298754e-3; 
Myomo=1.381982e-1;
SRBo=2.143165e-3; 
SLLjo=9.566355e-3; 
SLLslo=1.110363e-1; 
SLHjo=7.347888e-3; 
SLHslo=7.297378e-2; 
Csqnbo= 1.242988;
Ca_sro=0.1e-1; %5.545201e-1; 
Najo=9.136;%8.80329; 
Naslo=9.136;%8.80733; 
Naio=9.136;%8.80853; 
Kio=120; 
Cajo=1.737475e-4; 
Caslo= 1.031812e-4; 
Caio=8.597401e-5; 
Vmo=-8.09763e+1; 
rtoso=0.9946; 
ICajuncinto=1; 
ICaslinto=0;
C1o=0.0015;       % [] 
C2o=0.0244;       % [] 
C3o=0.1494;       % [] 
C4o=0.4071;       % [] 
C5o=0.4161;       % [] 
C7o=0.0001;       % [] 
C8o=0.0006;       % [] 
C9o=0.0008;       % [] 
C10o=0;           % [] 
C11o=0;           % [] 
C12o=0;           % [] 
C13o=0;           % [] 
C14o=0;           % [] 
C15o=0;           % [] 
O1o=0;            % [] 
O2o=0;            % [] 
C6o=1-(C1o+C2o+C3o+C4o+C5o+C7o+C8o+C9o+C10o+C11o+C12o+C13o+C14o+C15o+O1o+O2o);       % []


%initial conditions for IKur
rkuro = 0;
skuro = 1.0;

% Gating variables      
%   1       2       3       4       5       6       7       8       9       10      11      12      13
%%   m       h       j       d       f       fcaBj   fcaBsl   xtos    ytos    xtof    ytof    xkr     xks   
%y10=[1.2e-3;0.99;   0.99;   0.0;    1.0;    0.0141; 0.0141;     0;      1;      0.0;    1.0;    0.0;    6e-3;];
y10=[mo; ho; jo; do; fo; fcaBjo; fcaBslo; xtoso; ytoso; xtofo; ytofo; xkro; xkso;];   
% RyR and Buffering variables
%   14      15      16      17      18      19      20      21      22      23      24
%%   RyRr    RyRo    RyRi    NaBj    NaBsl   TnCL    TnCHc   TnCHm   CaM     Myoc    Myom  
y20=[RyRro; RyRoo; RyRio; NaBjo; NaBslo; TnCLo; TnCHco; TnCHmo; CaMo; Myoco; Myomo;];           
%y20=[1;     0;      0;      1.8;   0.8;    0.012;   0.112;  0.01;   0.4e-3; 1.9e-3; 0.135;];
% More buffering variables
%   25      26      27      28      29      30
%%   SRB     SLLj   SLLsl    SLHj    SLHsl  Csqnb
y30=[SRBo; SLLjo; SLLslo; SLHjo; SLHslo; Csqnbo];
%y30=[3.3e-3; 0.012; 0.012; 0.13;  0.13;  1.5;];
%   Intracellular concentrations/ Membrane voltage
%    31      32      33      34      35      36      37     38     39    40   41
%%    Ca_sr   Naj     Nasl    Nai     Ki      Caj    Casl    Cai   Vm  rtos ?
y40=[Ca_sro; Najo; Naslo; Naio; Kio; Cajo; Caslo; Caio; Vmo; rtoso; 1]; 
y50=[C1o; C2o; C3o; C4o; C5o; C6o; C7o; C8o; C9o; C10o; C11o; C12o; C13o; C14o; C15o; O1o];
% y40=[0.9;    8.8;    8.8;    8.8;    135;    0.1e-3; 0.1e-3; 0.1e-3; -88;  0.89; 0;          0;];
% y50=[UIC3o; UIC2o; UIFo; UIM1o; UC3o; UC2o; UC1o; UOo; UIM2o; LC3o; LC2o; LC1o; LOo ];    

% SVP put initial variables for IKur; 11/12/09
%(y(58); y(59));
y60=[rkuro; skuro];

% EG put initial variables for INaL; 
%(y(60); y(61));
y70=[1; 0; 0];

% Put everything together
%y0  = [y10;y20;y30;y40;y50]    
 y0  = [y10;y20;y30;y40;y50;y60;y70];
%load('yfinal'); % load output of previous simulation saved as yfinal.mat
%y0 = yfinal;

%%
%% Single Run Simulation
tspan = [0;1e3];
%%
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 
[t,y] = ode15s(@f,tspan,y0,options,p);

yfinal = y(end,:);
output = yfinal;
save 'yfinal'

currents = calcCurrents(t,y,p);

figure(1);
subplot(4,1,1); hold on,plot(t,y(:,39)); title('Voltage');
subplot(4,1,2); hold on,plot(t,y(:,38)); title('Cyto calcium');
subplot(4,1,3); hold on,plot(t,y(:,32)'); title('Na conc'); %,t,y(:,33),t,y(:,34)
subplot(4,1,4); hold on,plot(t,y(:,31)); title('Ca conc');%t,y(:,36),t,y(:,37),


function output = f(t,y,p,runType)

ydot = zeros(size(y));

%% Model Parameters
%% EPI or ENDO?
epi=1;
%% AF
AF=0;
%% ISO
ISO=0;
%% Right ATRIUM
RA=0;


% Constants
R = 8314;       % [J/kmol*K]  
Frdy = 96485;   % [C/mol]  
Temp = 310;     % [K]
FoRT = Frdy/R/Temp;
Cmem = 1.1e-10;   % [F] membrane capacitance 1.3810e-10;%
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100;     % cell length [um]113;%100
cellRadius = 10.25;   % cell radius [um]12;%10.25
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]
distSLcyto = 0.45;    % dist. SL to cytosol [um]
distJuncSL = 0.5;  % dist. junc to SL [um]
DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec] 
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; Vsr = 0.035*Vcell; Vsl = 0.02*Vcell; Vjunc = 1*0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius;  % [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2]
%J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
%J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 9.9664e-014
%J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.7460e-012
%J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 6.6240e-013
%J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.5618e-011
% tau's from c-code, not used here
J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments
Fjunc = 0.11;   Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations     
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 1;    % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

% Na transport parameters
      
%%
GNa=23*(1-0.1*AF);  % [mS/uF]
GNaB = 0.597e-3;    % [mS/uF] 
IbarNaK = 1.26;     % [uA/uF]
KmNaip = 11*(1-0.25*ISO);         % [mM]11
KmKo =1.5;         % [mM]1.5
Q10NaK = 1.63;  
Q10KmNai = 1.39;

%% K current parameters
pNaK = 0.01833;      
gkp = 0.002;

% Cl current parameters
GClCa =0.0548;   % [mS/uF]
GClB = 9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]

% I_Ca parameters
pNa = (1+0.5*ISO)*(1-0.5*AF)*0.75e-8;       % [cm/sec]
pCa = (1+0.5*ISO)*(1-0.5*AF)*2.7e-4;       % [cm/sec]
pK = (1+0.5*ISO)*(1-0.5*AF)*1.35e-7;        % [cm/sec]
Q10CaL = 1.8;       

%% Ca transport parameters
IbarNCX = (1+0.4*AF)*3.15;      % [uA/uF]5.5 before - 9 in rabbit
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact =0.384e-3;   % [mM] 0.256 rabbit
Q10NCX = 1.57;      % [none]
IbarSLCaP =  0.0471; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
KmPCa =0.5e-3;     % [mM] 
GCaB = 6.0643e-4;    % [uA/uF] 3
Q10SLCaP = 2.35;    % [none]

% SR flux parameters
Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = (2.5-1.25*ISO)*0.246e-3;          % [mM] default
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25;                 % [1/ms]      
koCa = 10+20*AF+10*ISO*(1-AF);               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]

% Buffering parameters
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
Bmax_Naj = 7.561;       % [mM] % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = (1+0.5*ISO)*19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 


%% Membrane Currents
mss = 1 / ((1 + exp( -(56.86 + y(39)) / 9.03 ))^2);
taum = 0.1292 * exp(-((y(39)+45.79)/15.54)^2) + 0.06487 * exp(-((y(39)-4.823)/51.12)^2);                 
 
ah = (y(39) >= -40) * (0)... 
   + (y(39) < -40) * (0.057 * exp( -(y(39) + 80) / 6.8 )); 
bh = (y(39) >= -40) * (0.77 / (0.13*(1 + exp( -(y(39) + 10.66) / 11.1 )))) ...
   + (y(39) < -40) * ((2.7 * exp( 0.079 * y(39)) + 3.1*10^5 * exp(0.3485 * y(39)))); 
tauh = 1 / (ah + bh); 
hss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);
 
aj = (y(39) >= -40) * (0) ...
    +(y(39) < -40) * (((-2.5428 * 10^4*exp(0.2444*y(39)) - 6.948*10^-6 * exp(-0.04391*y(39))) * (y(39) + 37.78)) / ...
                     (1 + exp( 0.311 * (y(39) + 79.23) )));
bj = (y(39) >= -40) * ((0.6 * exp( 0.057 * y(39))) / (1 + exp( -0.1 * (y(39) + 32) ))) ...
   + (y(39) < -40) * ((0.02424 * exp( -0.01052 * y(39) )) / (1 + exp( -0.1378 * (y(39) + 40.14) ))); 
tauj = 1 / (aj + bj);
jss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);         
 
ydot(1) = (mss - y(1)) / taum;
ydot(2) = (hss - y(2)) / tauh;
ydot(3) = (jss - y(3)) / tauj;
    
I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc);
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl);
I_Na = I_Na_junc+I_Na_sl;


% Late I_Na
GNaL=0.0025*AF;
aml = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bml = 0.08*exp(-y(39)/11);
hlinf = 1/(1+exp((y(39)+91)/6.1));
tauhl=600;
ydot(60) = aml*(1-y(60))-bml*y(60);
ydot(61) = (hlinf-y(61))/tauhl;

I_NaL_junc = Fjunc*GNaL*y(60)^3*y(61)*(y(39)-ena_junc);
I_NaL_sl = Fsl*GNaL*y(60)^3*y(61)*(y(39)-ena_sl);
I_NaL = I_NaL_junc + I_NaL_sl;
if t<9050
    ydot(62)=0;
else
    ydot(62)=I_NaL;
end
% I_nabk: Na Background Current
I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = 1*Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = 1*Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;

%% I_kr: Rapidly Activating K Current
gkr =0.035*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+10)/5));
tauxr = 550/(1+exp((-22-y(39))/9))*6/(1+exp((y(39)-(-11))/9))+230/(1+exp((y(39)-(-40))/20));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+74)/24));
I_kr = gkr*y(12)*rkr*(y(39)-ek);


%% I_ks: Slowly Activating K Current
markov_iks=0;
% pcaks_junc = -log10(y(36))+3.0; 
% pcaks_sl = -log10(y(37))+3.0;  
% gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
% gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));     

eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));

if markov_iks==0;
gks_junc=1*(1+1*AF+2*ISO)*0.0035*1;
gks_sl=1*(1+1*AF+2*ISO)*0.0035*1; %FRA
xsss = 1 / (1+exp(-(y(39)+40*ISO + 3.8)/14.25)); % fitting Fra
tauxs=990.1/(1+exp(-(y(39)+40*ISO+2.436)/14.12));
ydot(13) = (xsss-y(13))/tauxs;
I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
I_ks = I_ks_junc+I_ks_sl;
else
    gks_junc=1*0.0065;
    gks_sl=1*0.0065; %FRA
    alpha=3.98e-4*exp(3.61e-1*y(39)*FoRT);
    beta=5.74e-5*exp(-9.23e-2*y(39)*FoRT);
    gamma=3.41e-3*exp(8.68e-1*y(39)*FoRT);
    delta=1.2e-3*exp(-3.3e-1*y(39)*FoRT);
    teta=6.47e-3;
    eta=1.25e-2*exp(-4.81e-1*y(39)*FoRT);
    psi=6.33e-3*exp(1.27*y(39)*FoRT);
    omega=4.91e-3*exp(-6.79e-1*y(39)*FoRT);
    
ydot(42)=-4*alpha*y(42)+beta*y(43);
ydot(43)=4*alpha*y(42)-(beta+gamma+3*alpha)*y(43)+2*beta*y(44);
ydot(44)=3*alpha*y(43)-(2*beta+2*gamma+2*alpha)*y(44)+3*beta*y(45);
ydot(45)=2*alpha*y(44)-(3*beta+3*gamma+alpha)*y(45)+4*beta*y(46);
ydot(46)=1*alpha*y(44)-(4*beta+4*gamma)*y(46)+delta*y(50);    
ydot(47)=gamma*y(43)-(delta+3*alpha)*y(47)+beta*y(48);   
ydot(48)=2*gamma*y(44)+3*alpha*y(47)-(delta+beta+2*alpha+gamma)*y(48)+2*beta*y(49)+2*delta*y(51);
ydot(49)=3*gamma*y(45)+2*alpha*y(48)-(delta+2*beta+1*alpha+2*gamma)*y(49)+3*beta*y(50)+2*delta*y(52);
ydot(50)=4*gamma*y(46)+1*alpha*y(49)-(delta+3*beta+0*alpha+3*gamma)*y(50)+2*delta*y(53);
ydot(51)=1*gamma*y(48)-(2*delta+2*alpha)*y(51)+beta*y(52);  
ydot(52)=2*gamma*y(49)+2*alpha*y(51)-(2*delta+beta+1*alpha+gamma)*y(52)+2*beta*y(53)+3*delta*y(54);
ydot(53)=3*gamma*y(50)+1*alpha*y(52)-(2*delta+2*beta+2*gamma)*y(53)+3*delta*y(55);
ydot(54)=1*gamma*y(52)-(3*delta+1*alpha)*y(54)+beta*y(55);  
ydot(55)=2*gamma*y(53)+1*alpha*y(54)-(3*delta+1*beta+1*gamma)*y(55)+4*delta*y(56);
ydot(56)=1*gamma*y(55)-(4*delta+teta)*y(56)+eta*y(57);
O2=1-(y(42)+y(43)+y(44)+y(45)+y(46)+y(47)+y(49)+y(48)+y(50)+y(51)+y(52)+y(53)+y(54)+y(55)+y(56)+y(57));
ydot(57)=1*teta*y(56)-(eta+psi)*y(57)+omega*O2;
I_ks_junc = Fjunc*gks_junc*(y(57)+O2)*(y(39)-eks);
I_ks_sl = Fsl*gks_sl*(y(57)+O2)*(y(39)-eks);                                                                                                                                   
I_ks = I_ks_junc+I_ks_sl;
end
%I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(y(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(y(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

%% I_to: Transient Outward K Current (slow and fast components)
% modified for human myocytes

GtoFast=(1.0-0.7*AF)*0.165*1.0; %nS/pF maleckar; %human atrium

%11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
%atrium
%equations for activation; 
xtoss = ( (1)./ ( 1 + exp( -(y(39)+1.0)/11.0 ) ) );
tauxtof = 3.5*exp(-((y(39)/30.0)^2.0))+1.5;

%equations for inactivation;
ytoss = ( (1.0)./ ( 1 + exp( (y(39)+40.5)/11.5) ) ) ;
tauytof =25.635*exp(-(((y(39)+52.45)/15.8827)^2.0))+24.14;%14.14

ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = 1.0*GtoFast*y(10)*y(11)*(y(39)-ek);

I_to = 1*I_tof;

%% I_kur: Ultra rapid delayed rectifier Outward K Current
%Equation for IKur; from Maleckar et al. 2009 - EG
%atrium
%equations for activation;
Gkur = 1*(1.0-0.5*AF)*(1+2*ISO)* 0.045*(1+0.2*RA); %nS/pF maleckar 0.045
xkurss = ( (1)./ ( 1 + exp( (y(39)+6)/-8.6 ) ) );
tauxkur = 9/(1+exp((y(39)+5)/12.0))+0.5;

%equations for inactivation;
ykurss = ( (1)./ ( 1 + exp( (y(39)+7.5)/10 ) ) );
tauykur = 590/(1+exp((y(39)+60)/10.0))+3050;

ydot(58) = (xkurss-y(58))/tauxkur;
ydot(59) = (ykurss-y(59))/tauykur;
I_kur = 1*Gkur*y(58)*y(59)*(y(39)-ek);

%% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss = aki/(aki+bki);

%I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*(y(39)-ek);
%SVP 11/11/09
%multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
%resting potential
I_ki =(1+1*AF)*0.0525*sqrt(Ko/5.4)*kiss*(y(39)-ek);

% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;
I_Clbk = GClB*(y(39)-ecl);

GClCFTR=0;%4.9e-3*ISO;     % [mS/uF]
I_ClCFTR = GClCFTR*(y(39)-ecl);

%% I_Ca: L-type Calcium Current
dss = 1/(1+exp(-(y(39)+3*ISO+9)/6)); %in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
taud = 1*dss*(1-exp(-(y(39)+3*ISO+9)/6))/(0.035*(y(39)+3*ISO+9)); 
fss = 1/(1+exp((y(39)+3*ISO+30)/7))+0.2/(1+exp((50-y(39)-3*ISO)/20)); % in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+3*ISO+25))^2 )+0.02);
ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-1*11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-1*11.9e-3*y(7); % fCa_sl
fcaCaMSL= 0.1/(1+(0.01/y(37)));
fcaCaj= 0.1/(1+(0.01/y(36)));
fcaCaMSL=0;
fcaCaj= 0;
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
I_Ca_junc = (Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_Ca_sl = (Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45;
I_CaNa_junc = (Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_CaNa_sl = (Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(Kdact/y(36))^2);
Ka_sl = 1/(1+(Kdact/y(37))^2);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37);


I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx = I_ncx_junc+I_ncx_sl;

% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

% I_cabk: Ca Background Current
I_cabk_junc = Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = (1)*koCa/kCaSR;%
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36));          % [mM/ms]

J_serca = 1.0*Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);
J_SRleak = (1)*(1.0+0.25*AF)*5.348e-6*(y(31)-y(36));           %   [mM/ms]


%% Sodium and Calcium Buffering
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
J_CaB_cytosol = sum(ydot(19:25));

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current
% ydot(31)=0;

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc+I_NaL_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl+I_NaL_sl;   % [uA/uF]
I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
   +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
%FluxNaSL=ydot(33);
% ydot(32) = 0;
% ydot(33) = 0;
ydot(34) = J_na_slmyo/Vmyo*(y(33)-y(34));             % [mM/msec] 
% ydot(34)=0;

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur;     % [uA/uF] %SVP: added IKur
% ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]
ydot(35) =0; % -I_K_tot*Cmem/(Vmyo*Frdy);

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
% ydot(36)=0;
% ydot(37)=0;
% ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;%+J_ca_slmyo/Vmyo*(y(37)-y(38));    % [mM/msec]
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y(37)-y(38));
% ydot(38)=0;


%% Simulation type
protocol = 'pace1';

switch lower(protocol)
    case {'none',''},
        I_app = 0;
    case 'pace1',        % pace w/ current injection at rate 'rate'
		%rate =  1e-3;
        rate =1e-3;  %SVP 11/12/09 in 4Hz?
        
       		if mod(t,1/rate) <= 5
            I_app = 12.5;
        else
            I_app = 0.0;
        end
%     case 'pace2',        % pace w/ current injection at rate 'rate'
% 		factor = 2;
%         rate = factor*1e-3;
% 		if (mod(t+900,1/rate) <= 5) & ((t > 5000) & (t < 10000))
%             I_app = 12.0;
%         elseif (mod(t+900,1/rate*factor*2) <= 5)  & ((t <= 5000) | (t >= 10000))
%             I_app = 12.0;
%         else
%             I_app = 0;
%         end    
    case 'pacenew',        % pace w/ current injection at rate 'rate'
		if (t>=0)& (t <= 5)
            I_app = 9.5;
        else
            if (t>=310)& (t <= 315)
            I_app = 9.5;
            else
            I_app = 0.0;
            end
        end
    case 'pace4',
        rate = 40e-3;
        
        
        if (t > 1000) & (t < 2000) & (mod(t+900,1/rate) <= 5) 
                    I_app = 10;    
        elseif (t > 2000) & (mod(t+900,1/rate*10) <= 5)
                    I_app = 10;
                    else
            I_app = 0;
        end
    case 'vclamp',      
		V_hold = -80;           
        V_test = 60;
		if (t > 1000 & t < 1500)
		    V_clamp = V_test;
		else
		    V_clamp = V_hold;
		end
		R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
%           case 'caffeine',      
% 		V_clamp = -80;           
%         R_clamp = 0.02;
% 		I_app = (V_clamp-y(39))/R_clamp;
    case 'iv',
        rate = 0.5e-3; %frequency of protocol sweep in kHz
        
        if mod(t+1900,1/rate)<=100 %duration of the voltage step
            V_clamp = 10;%+floor(t*rate)*10;
		else
		   if mod(t,1/rate)<=100 %duration of the voltage step
            V_clamp = -80+40/100*(t-floor(t*rate)*2000);
        else
		    V_clamp = -80;
           end
        end
          R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
        case 'iv2',
        rate = 1e-3; %frequency of protocol sweep in kHz
        
        if mod(t,1/rate)<=300 %duration of the voltage step
            V_clamp = -30;%+floor(t*rate)*10;
		
        else
		    V_clamp = -120;
           
        end
% V_clamp = -80;
        R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
    case 'recovery',
        rate = 0.1e-3;
        %tdelay=[7.5 15.1 22.6 37.7 67.9 128.3 188.7 301.9 400.0 505.7 603.8 709.4 784.9 890.6 1003.8 1200 1403.8 1607.5 1796.2 2015.1];
      %  tdelay=[10 30 50 80 100 200 300 500 800 1200 2000 5000]; %interpulse intervals
        tdelay=[10 1 5 80];
        if (mod(t,1/rate)<=1000 | mod(t-1000+1/rate-tdelay(floor(t*rate)+1),1/rate)<=1000)
            V_clamp = -20;
        else
		    V_clamp = -80;
        end
        R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
        case 'envelope',
        rate = 0.5e-3; %frequency of protocol sweep in kHz
        %if mod(t,1/rate)<=10+floor(t*rate)*10 %duration of the voltage step
        if mod(t+1800,1/rate)<=p
            V_clamp = 30;
		else
		    V_clamp = -40;
        end
        R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
            
end  

%% Membrane Potential
%%
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk++I_ClCFTR;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
%ydot(39) = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);
ydot(39) = -(I_tot-I_app);
vmax = ydot(39);
% ----- END EC COUPLING MODEL ---------------
% adjust output depending on the function call

if (nargin == 3)
    output = ydot;
elseif (nargin == 4) & strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) & strcmp(runType,'rates')
    output = r;
elseif (nargin == 4) & strcmp(runType,'currents')
    %currents = [I_Na I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_Catot I_ncx I_pca I_cabk J_serca*Vmyo/Vsr];
    %currents = [I_Na I_tof I_tos I_kr I_ks I_ClCa I_Catot J_SRCarel*Vsr/Vmyo J_SRleak RI I_ncx]; 
%     currents= [I_Catot J_SRCarel*2*Vsr*Frdy/Cmem];
% currents=[I_Na I_Catot I_ncx I_nak I_kr I_ks I_tof I_tos I_ki];
currents=[];
    output = currents;
end

%% Calculate timecourse for currents and other intermediates
function currents = calcCurrents(t,y,p)
% After running a simulation, feed the time vector and state variables into
% this function to compute ionic currents, etc.
% currents: [I_Na,I_Catot];
currents=[];
for i=1:size(t)
    if ceil(i/1000)==i/1000
%         disp(['t = ',num2str(ceil(t(i)))]);
    end
    currents=[currents;f(t(i),y(i,:),p,'currents')];
end

% end calcCurrents