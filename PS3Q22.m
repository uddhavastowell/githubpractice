%Problem Set 3, Question 2
%PS3_Q2.m
%Vincent U Stowell
%09/19/2019
%Code for solving mRNA dynamic equations (using the Newton Method)
%Here the initial concentration is 0nM, 
%The decay constant is 0.15 /min
%The rate of transcription a, is oscillating between 15 and 0 nM/min
%And t(1/2)=ln(2)/k4

clear
close all

save;

%Input Section:
%Model Parameters 

%Initial Condition
m0=0; % Initial mRNA Concentration in nM

a0=15; %Initial rate of mRNA transcription in nM/min

k4=0.15; %Decay constant of mRNA transcript mi per min

%thalf=log(2)/k4; 
thalf=4.62;

%Integration Parameters
dt=1/60; % Length of time step in minutes
T=200; % Length of time over which to integrate
N=T/dt+1; % Number of time steps


%%%%%%%%%%%%

%Allocate output matrices
m=zeros(1,N); %mRNA concentration
m(1)=m0; % <Input Initial concentration
t=[0:dt:T]; %Create time vector

a=zeros(1,N);% Decay constant of mRNA in minutes

c=ceil(t/thalf);% Define an array by vector t, in which each point is rounded up to the nearest
%multiple of thalf, so that the dataset proceeds in whole steps of thalf.

d=c/2;% Divide dataset c so that half of the points have a decimal remainder.

f= (d~=round(d));% Answer 1 if the datapoint from d does not equal itself when rounded; 
% Converts odd data points from c (that is odd multiples of thalf) into 1s,
% and even points into 0s.

a=f*a0;% Apply 15nM rate of transcription when transcription is active


for n=1:N-1
    dm=a(n)-k4*m(n)*dt;% find change in mRNA concentration during dt, based on current mRNA concentration
    m(n+1)=m(n)+dm;% add change in mRNA concentration to current mRNA concentration
    if c(n)==20
x=m(n);% Define x as the mRNA concentration at the onset of the 10th cycle
elseif c(n)==21
y=m(n);% Define y as the mRNA concentration at the close of the 10th cycle
    end
   
end



v=zeros(1,2);
v(1)=x;
v(2)=y;% Define v as a vector of points x and y

b=ones(1,N);
z=b*mean(v);% Find the mean of the mRNA concentration at onset and close of the 10th cycle
% and define vector z as equal to this mean at all points from 1 to N
    


%Plot the Results
figure
hold on
plot(t,m)
plot(t,z)
hold off
xlabel('time (min)')
ylabel('mRNA Concentration (nM)')


