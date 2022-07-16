clc
clear all

%%%%%%%%%----The below code calculates stability of the reduced socio-mutualisitc network-----%%%%%%%%%%%%%
% param a: intrinsic growth rate
% param h: handling time
% param k1: degree of pollinator
% param k2: degree of plant
% param p: mutualisitc trade-off
% param mu: immigration term
% param g: interaction strength
% param d: strength of norm
% param c: rarity term
% param sigma: learning rate
% param l: cost of conservation
load Network4.mat
B=gammaAP;
d=0.5;c=0.5;h=0.4;alpha=0.1;beta=1;p=0.5;

% % %%%%%%%% eigenvector weighted%%%%%%%%%%%%
% [VP DP]=eigs(B*B');  
% [VA DA]=eigs(B'*B); 
% 
% gammaA=sum(g*k1.^(1-p).*VP(:,1))/sum(VP(:,1));
% gammaP=sum(g*k2.^(1-p).*VA(:,1))/sum(VA(:,1));

gammaA=2.5107;gammaP=2.6396; % calculated from eigen-vector weighted method

A6=[];A7=[];
k1=0:0.1:1;  
l1=0:0.1:1.5;
sigma=0.18;
rt=[];rt2=[];
for ii=1:length(k1)
    for iii=1:length(l1)
        l=l1(iii);
    k=k1(ii);
q1=2*d*beta*h*gammaA*gammaP*(1+alpha*h);
q2=(c*2*d*beta*h*gammaA*gammaP*(1+alpha*h))+((2*d*beta^2)+(2*h*gammaA*alpha*d*beta))-((h*gammaA*gammaP+alpha*(h^2)*gammaA*gammaP)*(d*(2*alpha-k)+k*l))-(2*d*gammaA*gammaP*(1+alpha*h));
q3=c*(2*d*beta^2+h*gammaP*alpha*2*d*beta-(h*gammaA*gammaP+alpha*h^2*gammaA*gammaP)*(d*(2*alpha-k)+k*l)-2*d*gammaA*gammaP*(1+alpha*h))-2*d*gammaA*alpha-(d*(2*alpha-k)+k*l)*(beta+h*gammaA*alpha)-(h*gammaA*gammaP*(1+alpha*h));
q4=-((c*(d*(2*alpha-k)+k*l)*(beta+h*gammaA*alpha))+(2*d*c*gammaA*alpha)-(beta+h*gammaA*alpha));
%     function sols = solve_cubic(q1,q2,q3,q4)
%     syms x
%     sols = solve(q1*x^3 + q2*x^2 + q3*x + q4);
%     end
    
    x = roots([q1 q2 q3 q4]);
    rt=[rt;x'];
    rt1=[];
    for i=1:3
        if imag(rt(i))==0 && rt(i)>0
           rt1=[rt1;rt(i)];
        
        end
    end
%     rt2=[rt2; k rt1'];
% end
        A=max(rt1);
    P=(alpha+((gammaP*A)/(1+h*gammaP*A)))/beta;

x=((A+c)*(d+l)-1)/((A+c)*d*2);
a11=alpha-(2*beta*P)+(gammaP*A)/(1+h*gammaP*A);
a12=(gammaP*P)/((1+h*gammaP*A).^2);
a13=0;
a21=(gammaP*A)/((1+h*gammaA*P)^2);
a22=(alpha-k*(1-x))-(2*beta*A)+((gammaA*P)/(1+h*gammaA*P));
a23=k*A;
a31=0;
a32=-sigma*x*(1-x)/((A+c)^2);
a33=-d*(6*x^2-6*x+1)+(sigma*(1-2*x)*(1-(l*(A+c)))/(A+c));
jacob=[a11 a12 a13;a21 a22 a23;a31 a32 a33];

jacob1=eig(jacob);
%      A6=[A6;  k real(jacob1(1)) real(jacob1(2)) real(jacob1(3)) max(real(jacob1))];
%      v1(iI,j)=real(jacob1(1));
%      v2(i,j)=real(jacob1(2));
     if max(real(jacob1))<0
    
    A6=[A6;k l A P x];
     else 
         A7=[A7;k l A P x]
     end
    
    
end

     
end
     
     
     
     
     
     
     
