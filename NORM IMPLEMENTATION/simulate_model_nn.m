clc
clear all
%%%%%%%%%%%%%%%--The below code generates time series for a given interaction matrix at a fixed value of driver parameter in presence and absence of social norm---%%%%%%%%%%%%%%%%%
% param a: intrinsic growth rate
% param h: handling time
% param k1: degree of pollinator
% param k2: degree of plant
% param p: mutualisitc trade-off
% param mu: immigration term
% param g: interaction strength
% param d: strength of norm
% param m3: rarity term
% param k3: learning rate
% param l: cost of conservation
a=0.1;mu=0.0001;h=0.4;p=0.5;
q=[];q1=[];
c1=[];c2=[];

 for pp=1
  
 filename=append('Network',int2str(pp),'.mat');
 load(filename)
 B=gammaAP;
[n m]=size(B);                                        
 for i=1:n
     for j=1:m
 if B(i,j)>0
     B(i,j)=1;
 else B(i,j)=0;
 end
     end
 end

b=eye(m);
b1=eye(n);
b=eye(m);
b1=eye(n);
for i=1:m
    for j=1:m
        if i==j
        b(i,j)=1;
        else
            b(i,j)=0.0;
        end
    end
end
for i=1:n
    for j=1:n
        if i==j
        b1(i,j)=1;
        else
            b1(i,j)=0.0;
        end
    end
end
k1=sum(B,1);
k2=sum(B,2);
g=1;
for ii=1:m
    B1(:,ii)=(B(:,ii)./(k1(ii)^p))*g;
end
for ii=1:n
    B2(ii,:)=(B(ii,:)./(k2(ii)^p))*g;
end

% ts=0:0.05:3;
y0=[];
 y0=[rand(m,1); rand(n,1)];
 y0=y0';
 y0=reshape(y0,[1,m+n]);
 p0=y0(1:m);
q0=y0(m+1:m+n);
x=p0;       % initial pollinator abundance
y=q0;       % initial plant abundance
z=0.0001*rand(m,1)';   %initial norm abundance
k3=0.18;
d=0.5;
m3=0.5;
l=0.14;
t=0;t_max =500; dt=0.05;
m1=t_max/dt;
 
k=0;
k_max=1;
k_n=10;
dk=(k_max-k)/k_n;
%np=m/2;  the value of np may be changed depending on the fraction of species we want to apply norm for conservation
% [r I1]=mink(k1,np);  % for applying norm at a fraction of specialist, for generalist replace min by max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for jj=1:k_n+1
  t=0.0;
  while t<t_max
         for j=1:m
%          for i=1:n
    c1(j)=B1(:,j)'*y';
    c1(j)=c1(j)/(1+h*c1(j));   % growth due to mutualism for pollinators
 end
for j=1:n
    c2(j)=B2(j,:)*x';
    c2(j)=c2(j)/(1+h*c2(j));   % growth due to mutualism for plants
end
         B3=b*(x'.*x);       % competition faced by pollinators
         B4=b1*(y'.*y);      % competition faced by plants

 % ================================================================================================= 
                                 %%%---without norm---%%%
 % =================================================================================================

% %        
%  x=x+ (a.*x-k.*x-diag(B3)'+mu+c1.*x)*dt;
%  y=y+(a.*y-diag(B4)'+mu+c2.*y)*dt;
% % % %      

 % ================================================================================================= 
                                 %%%---norm at all pollinator nodes---%%%
 % =================================================================================================


B5=diag(B3)';
B6=diag(B4)';
x=x+ (a*x-k*(1-z).*x-diag(B3)'+mu+c1.*x)*dt;
y=y+(a.*y-diag(B4)'+mu+c2.*y)*dt;
z=z+(k3*z.*(1-z).*((d*(2*z-1))+(1./(x+m3))-l))*dt;
 % ================================================================================================= 
            %%%---fraction of specialist orgeneralist---%%%
           %%%---for running the section below please uncomment lines 87-88---%%%
 % =================================================================================================


% for ii=1:m
%   
%     if find(any(ii==(I1(1,:))))==1
% 
%     x(ii)=x(ii)+ (a*x(ii)-(k*(1-z(ii)).*x(ii))-B5(ii)+mu+(c1(ii).*x(ii)))*dt;
%     z(ii)=z(ii)+(k3*z(ii).*(1-z(ii)).*(d*(2*z(ii)-1)+(1./(x(ii)+m3))-l))*dt;
%     else 
%         x(ii)=x(ii)+(a.*x(ii)-k.*x(ii)-B5(ii)+mu+(c1(ii).*x(ii)))*dt;
%         z(ii)=0.0;
% 
%     end
%  
%   end
%       y=y+(a.*y-B6+mu+c2.*y)*dt;
% 
%  
%    x=abs(x);
%   y=abs(y);
  t=t+dt;
  end
  q=[q;k x y z ];
  
 k=k+dk;    

   save(append('data_norm_opt_pol',int2str(pp)),'q','-v7.3');
 end
 end
figure
 plot(q(:,1),q(:,2:26));
 hold on
 plot(q(:,1),q(:,52:76));
