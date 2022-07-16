clc
clear all

%%%%---This code generate nested matrices of desired dimension and connectance---%%%%%%
%%%%%---This code requires the function network_detail.m for compilation---%%%%%
SA = 25;
SP= 25;
cM=0.15;
B=[];
NODF=0.2;

         for ii=1:5000000
         
         
       gammaAP=zeros(SA,SP);
       gammaPA=zeros(SP,SA);  
         for i=1:SA
            for j=1:SP  
                 if rand<cM;
                    gammaAP(i,j)=1;   
                    gammaPA(j,i)=gammaAP(i,j);   %create random adjacency matrix of desired connectance cM and dimension

                 end
            end
         end
        ind = find(sum(gammaAP,1)==0);
        ind1 = find(sum(gammaAP,2)==0);
        if length(ind)==0 && length(ind1)==0
            
         k1=sum(gammaAP,1);  %%%SP
         k2=sum(gammaAP,2);  %%%SA
         A=network_detail(gammaAP);
             
             idX=randi(SA);
            
             idX_z=find(gammaAP(idX,:)==0);
             idX_z_new=idX_z(randi(length(idX_z)));
             idX_z1=find(gammaAP(idX,:)>0);
               
             idX_z_old=idX_z1(randi(length(idX_z1)));
             if k1(idX_z_new)> k1(idX_z_old) && k1(idX_z_old)>1
                 gammaAP(idX,idX_z_new)=1;
                 gammaAP(idX,idX_z_old)=0;
             end
             A=network_detail(gammaAP);    
             if abs(NODF-A)<0.05
                save(append('Network',int2str(ii)),'gammaAP','-v7.3');
             %break;
             end
         end
%         end
%          gammaAP1=gammaAP;
%          gammaAP=gammaAP1;
        end

