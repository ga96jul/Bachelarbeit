varx = 1:0.01:10;

%entropy Y
entY = log(pi*exp(1+varx));

%entropy N
entN = log(pi*exp(1));

%Capacity

CapB = entY - entN;
%CapC = log10(exp(varx)); 



figure;
hold on;
%plot(varx, CapB);
grid on;
legend;

%% QPSK

x = [1 1; -1 1; -1 -1; -1 1];
SNR = 10;
N = 1000;
symbols = N/2;
data = [1,1,-1,1,-1,-1,-1,1];
fData = N/length(data);
Data = repmat(data,1,fData);
br = symbols*8;
f = N;
T = 1/(br);
t = 1/br:1/br:1;

s_p_data = reshape(Data,2,length(Data)/2);

inphase_data = s_p_data(1,:);
str_inphase_data = imresize(inphase_data, [1 br], 'nearest');
Qphase_data = s_p_data(2,:);
str_Qphase_data = imresize(Qphase_data, [1 br], 'nearest');

%oddBits = Data(1:2:end);
%evenBits = Data(2:2:end);
%plot(1:length(oddBits),oddBits);
% multiply with cosine
y = [];
calcY = [];


    y1=str_inphase_data.*cos((pi/4)+2*(pi)*f*t); % inphase component
    y2=str_Qphase_data.*sin((pi/4)+2*(pi)*f*t);% Quadrature component
    y=[y1+j*y2]; %modulated bitdata 

y = awgn(y,SNR,'measured');
plot(t,y1);
plot(t,y2);
plot(t,y);
H = scatterplot(y);
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
grid on;

%% Capacity of QPSK

entropy_noise = entN; %h(N)

%entropy of signal
for k = 1:length(Data)/2
    for p = 1:4
        entY2 = 0.25*(1/pi)*exp(-(y(k)-x(i,:)));
    end
end
entFinal = entY2 - entropy_noise;


plot(10,entFinal);
    
    
    
    
