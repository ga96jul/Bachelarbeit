function Rf=Gen_Rayleigh_Ch(NofSym,M,Speed)
%% Rayleigh Fading Simulator as the Jakes’ method
Rc = 3.84e6;
Rb = Rc/M;
fD = Doppler(Speed);
%==========================================================================
% fading génération
N = 34;
M  = (N-2)/4;
xc = zeros(NofSym,1);
xs = zeros(NofSym,1);
t = 0:(NofSym-1);
t = t'/(Rb)+ 5000/Rb; %randint(1,1,NofSym)+100/Rb;;%
%------------------------
B0 = pi/4;
a0 = sqrt(2)*cos(B0);
b0 = sqrt(2)*sin(B0);
w0 = 2*pi*fD;
%------------------------
xc = a0.*cos(w0*t);
xs = b0.*cos(w0*t);
%------------------------
phi = 0;
for n = 1:M
    %------------------------
    Bn = n*pi/M;
    an = 2*cos(Bn);
    bn = 2*sin(Bn);
    wn = w0*cos((2*pi*n)/N);
    %------------------------
    xc = xc + an.*cos(wn*t);
    xs = xs + bn.*cos(wn*t);
end;
Rf = (xc +1i.*xs).*(2/sqrt(N));
Rf = Rf.' ;
% figure;
% plot(t,10*log10(abs(Rf)));grid;
%==========================================================================
function FD = Doppler(v)
c = 3e8;      % is the speed of light (m/s).
f = 1.8e9;
v_ms = v.*(1000/3600);
FD =  v_ms*f/c;
%==========================================================================




