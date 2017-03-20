%C4 key parameters
f1= 262; %fundamental frequency
Vh0= 2; %initial hammer velocity  

%String parameters
L= [0.62,0.62,0.62]; %string length
Ms= [3.93/1000,(3.93/0.62)*L(2)/1000,(3.93/0.62)*L(3)/1000]; %string mass
T= [670,((f1*2*L(2))^2)*(Ms(2)/L(2)),((f1*2*L(3))^2)*(Ms(3)/L(3))]; %tension 
E= 2e11; %young's modulus
e= 3.82e-5; %stiffness
b1= 0.5; %damping coef (air)
b3= 6.25e-9; %damping coef (str)

%Hammer parameters
a= 0.12; %relative striking position
Mh= 2.97/1000; %hammer mass
p= 2.5; %stiffness exponent
K= 4.5e9; %stiffness coef
m_r= [0.75,Mh/Ms(2),Mh/Ms(3)]; %hammer string mass ratio

%Sampling parameters
Fs= 32000; %sampling frequency
N= 50; %number of string segments

dt= 1/Fs; %delta t
dx= [L(1)/N,L(2)/N,L(3)/N];
c=[sqrt(T(1)/(Ms(1)/L(1))),sqrt(T(2)/(Ms(2)/L(2))),sqrt(T(3)/(Ms(3)/L(3)))]; % Wave speed
i0= round(a*N);
%-------------------------------------------------%
%Defining the coefficients of the wave equation 
D=1+b1*dt+2*b3/dt;
r=[c*dt/dx(1),c*dt/dx(2),c*dt/dx(3)];
    
for k=1:3
a1(k)=(2-2*r(k)^2+b3/dt-6*e*N^2*r(k)^2)/D;
a2(k)=(-1+b1*dt+2*b3/dt)/D;
a3(k)=((r(k)^2)*(1+4*e*N^2))/D;
a4(k)=(b3/dt-e*N^2*r(k)^2)/D;
a5(k)=(-b3/dt)/D;
end

input_length=200000;
played_wave=zeros(1,input_length);
for k=1:3
%Variables initialization
input_length=200000;
y=zeros(N,input_length); % Displacement of the string
yh=zeros(1,input_length); % Displacement of the hammer
F=zeros(1,input_length); % Force signal output

%For first few timesteps
yh(1)=0;
y(:,1)=0;
F(1)= 0; 
yh(2)=Vh0*dt;
y(1,2)=0;
y(N,2)=0;
y(2:N-1,2)=(y(3:N,1)+y(1:N-2,1))/2;
F(2)=K*abs(yh(2)-y(i0,2))^p;
y(1,3)=0;
y(N,3)=0;
y(2:N-1,3)=y(3:N,2)+y(1:N-2,2)-y(2:N-1,1);
y(i0,3)=y(i0+1,2)+y(i0-1,2)-y(i0,1)+((dt)^2*N*F(2))/Ms(k);
yh(3)=2*yh(2)-yh(1)-(((dt)^2)*F(2))/Mh;
F(3)=K*abs(yh(3)-y(i0,3))^p;

%For the remaining time steps
for n=4:input_length
y(1,n)=0;
y(N,n)=0;
y(2,n)= a1(k)*y(2,n-1)+a2(k)*y(2,n-2)+...
a3(k)*(y(3,n-1)+y(1,n-1))+...
a4(k)*(y(4,n-1)-y(2,n-1))+...
a5(k)*(y(3,n-2)+y(1,n-2)+y(2,n-3));
y(N-1,n)= a1(k)*y(N-1,n-1)+a2(k)*y(N-1,n-2)+...
a3(k)*(y(N,n-1)+y(N-2,n-1))+...
a4(k)*(y(N-3,n-1)-y(N-1,n-1))+...
a5(k)*(y(N,n-2)+y(N-2,n-2)+y(N-1,n-3));
y(3:N-2,n)= a1(k)*y(3:N-2,n-1)+a2(k)*y(3:N-2,n-2)+...
a3(k)*(y(4:N-1,n-1)+y(2:N-3,n-1))+...
a4(k)*(y(5:N,n-1)+y(1:N-4,n-1))+...
a5(k)*(y(4:N-1,n-2)+y(2:N-3,n-2)+y(3:N-2,n-3));
y(i0,n)= a1(k)*y(i0,n-1)+a2(k)*y(i0,n-2)+...
a3(k)*(y(i0+1,n-1)+y(i0-1,n-1))+...
a4(k)*(y(i0+2,n-1)+y(i0-2,n-1))+...
a5(k)*(y(i0+1,n-2)+y(i0-1,n-2)+y(i0,n-3))+...
((1/Fs)^2*N*F(n-1))/Ms(k);
yh(n)=2*yh(n-1)-yh(n-2)-((1/Fs)^2*F(n-1))/Mh;
% Check for when the hammer is no longer in contact with the string %
if (yh(n)-y(i0,n))>0
F(n)=K*abs(yh(n)-y(i0,n))^p;
else
F(n)=0;
end
end

s_seg= 2;
sig(k,1:1:input_length)=y(s_seg,:);
if k == 3
figure;
hold on;
plot((1:1:input_length)*dt,sig*1000);
%sound(sig*1000,Fs)
end
played_wave= played_wave + y(s_seg,:); 
ypp= y;
%clear y yh F
end
figure;
plot((1:1:input_length)*dt, played_wave*1000);
sound(played_wave*1000,Fs)
%Need to check 


%allframes= cell(input_length);
%mov(1:input_length) = struct('cdata',[], 'colormap',[]);
%for i=1:input_length
%    fig= figure;
%    plot((1:1:50),y(:,50*i));
%    mo20v(i)= getframe(gcf);
%    close(fig);
%    sound(y(:,50*i))
    %allframes(i)= {F.cdata};
%end
%mov(1:input_length) = struct('cdata',[], 'colormap',[]);
%movie2avi(mov, 'myPeaks1.avi', 'compression','None', 'fps',100);
%winopen('myPeaks1.avi')

