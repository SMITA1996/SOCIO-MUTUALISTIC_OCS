clc
clear all
%%%%%---This code finds the optimal conservation strategy for each characteristic network---%%%%%
%%%%%---This code requires function network_detail.m for compilation---%%%%%
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

load Network1.mat
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
com2 = BipartiteModularity.LEADING_EIGENVECTOR(B);
car_M=com2.N;
M1=com2.row_modules;
M2=com2.col_modules;
C1=com2.index_rows;
C2=com2.index_cols;
k1=sum(B,1); 
k2=sum(B,2); 
I1=[];I2=[];
 k_M=zeros(1,m);
 k_mean=zeros(1,m);
 k_std=zeros(1,m);
 kk_mean=[];kk_std=[];
for j=1:length(k1)
    I1=M2(j,1);
    I11=find(M2(:,1)==I1);  %%%%no of pollinators in the same module as j
    I2=find(M1(:,1)==I1); %%%%no of plants in the same module as j 
    k=zeros(1,length(I11)); 
    for ii=1:length(I11)   %pollinator
        for i=1:length(I2) %plant
        if B(I2(i),I11(ii))==1
            k(ii)=k(ii)+1;
        end
        end
        k_mean=mean(k);
        k_sd=std(k);
    end
    
    kk_mean=[kk_mean;k_mean];
    kk_std=[kk_std;k_sd];
    for i=1:length(I2)
        if B(I2(i),j)==1
            k_M(j)=k_M(j)+1;   %%each pollinators in module degree
        end
    end
        
end
%%%%%%%%%%%%%%z score%%%%%%%%%%%%%
% z=[];
% k_M1=k_M';
% for i=1:length(k1)
%     z(i)=(k_M1(i)-kk_mean(i))/(kk_std(i));
% end
%%%%%%%%%%%%%participation coefficient%%%%%%%%%%%%%%
% IND1=[];
% deg=zeros(car_M,m);
% for j=1:length(k1)
%     
%     for i=1:car_M
%         
%         IND1=find(M1(:,1)==i);
%         for l=1:length(IND1)
%             if B(IND1(l),j)==1
%                 deg(i,j)=deg(i,j)+1;
%     end
%     
%         end
%     end
% end
% q=zeros(1,car_M);
% for j=1:length(k1)
%     q=deg(:,j).^2;
%     q1=k1(j)^2;
%     p(j)=1-(sum(q)./q1);
%    
% end
%%%%%%%%---find pollinators in each module---%%%%%%%%%
IND2=[];IND3=[];in_mod_deg=[];tot_deg=[];IND4=[];IND5=[];IND6=[];IND7=[];
F=[];
set1=find(k1(1,:)==1)
A1=zeros(m,car_M);
A2=zeros(m,car_M);
A3=zeros(m,car_M);
G=[];set_new=[];
for i=1:car_M
    IND2=find(M2(:,1)==i); %pollinator in ith module
%     IND3=find(M1(:,1)==i); %plant in ith module
    in_mod_deg=k_M(IND2);  %degree of pollinators in the ith module
    tot_deg=k1(IND2);
    for ii=1:length(IND2)
         A1(ii,i)=IND2(ii);
         A2(ii,i)=in_mod_deg(ii);
         A3(ii,i)=tot_deg(ii);
    end
    A4=A2(:,i)
    A5=A1(:,i);
    A6=A3(:,1);
     [x,y]=find(A4==max(A4));
     if length(x)>1
         [x1,y1]=(max(A6(x,:)));
         y1=x(y1);
     else
         x1=x;
         y1=y;
   
     end
        
 F=[F;A1(y1,i)];
end
 set2=set1';
     
set_new=[set_new;set1,F']; 
S_new=unique(set_new);
I1=S_new;
ext=1;
%%%%%%%%%---model simulation using set of nodes obtained by computation till step 5---%%%%%%%%%%%
a=0.1;mu=0.0001;h=0.4;p=0.5;
q=[];q1=[];
c1=[];c2=[];data=[];
entry=[];
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
g1=1;  
q=[];
g=1;  
for ii=1:m
    B1(:,ii)=(B(:,ii)./(k1(ii)^p))*g;
end
for ii=1:n
    B2(ii,:)=(B(ii,:)./(k2(ii)^p))*g;
end
 ext=1;
for iii=1:25
y0=[];
 y0=[rand(m,1); rand(n,1)];
 y0=y0';
 y0=reshape(y0,[1,m+n]);
 p0=y0(1:m);
q0=y0(m+1:m+n);
x=p0;
y=q0;
z=0.0001*rand(m,1)';
k3=0.18;
d=0.5;
m3=0.5;
l=0.14;
t=0;t_max =500; dt=0.05;
m1=t_max/dt;
 
k=0; % vary epsilon1/epsilon2 in the range [0.5 2.6] and y variable is v
k_max=1;
k_n=10;
dk=(k_max-k)/k_n;


 if length(ext)>0
for jj=1:k_n+1
  t=0.0;
  while t<t_max
         for j=1:m
%          for i=1:n
    c1(j)=B1(:,j)'*y';
    c1(j)=c1(j)/(1+h*c1(j));
 end
for j=1:n
    c2(j)=B2(j,:)*x';
    c2(j)=c2(j)/(1+h*c2(j));
end
         B3=b*(x'.*x);
         B4=b1*(y'.*y);

B5=diag(B3)'; 
B6=diag(B4)'; 
% ================================================================================================= 
                    %%%---Model Equation with norm at nodes obtained using OCS---%%%
% ================================================================================================= 


for ii=1:m
   %  for iii=1:length(I1)
%     if ii==I1(1,:)
    if find(any(ii==(I1(1,:))))==1

    x(ii)=x(ii)+(a*x(ii)-(k*(1-z(ii)).*x(ii))-B5(ii)+mu+(c1(ii).*x(ii)))*dt;
    z(ii)=z(ii)+(k3*z(ii).*(1-z(ii)).*(d*(2*z(ii)-1)+(1./(x(ii)+m3))-l))*dt;
    else
% %         3
        x(ii)=x(ii)+(a.*x(ii)-k.*x(ii)-B5(ii)+mu+(c1(ii).*x(ii)))*dt;
        z(ii)=0.0;
%         y(ii)=y(ii)+(a.*y(ii)-B6(ii)+mu+c2.*y(ii))*dt;

    end
 
  end
      y=y+(a.*y-B6+mu+c2.*y)*dt;

 
   x=abs(x);
  y=abs(y);
  t=t+dt;
  end
  q=[q;k x y z];


 k=k+dk;   
end

     xm=[];
        ext=find(q(end,2:m+1)<0.001); 
        for i=1:length(ext)
            [row(i),column(i)] = find(A1==ext(i));
             xn=A2(row(i),column(i));
             xm=[xm;xn]; 
        end          %in module degree of extinct species
             [m_p n_p]=max(xm);   % max inmodule degree of ext species
             I1=[I1,ext(n_p)];
       
        
           
             
        
        iii
 else 
      iii=25;
  end
end
save(append('optimal_set',int2str(1)),'I1','-v7.3');

