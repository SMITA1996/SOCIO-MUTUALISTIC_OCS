function [nodf,qb,Nm] = network_detail(gammaAP)

% Makes use of matlab library "BiMat" available via the link: https://bimat.github.io/. 

fp = Bipartite(gammaAP);

% Nestedness 
fp.nestedness = NestednessNODF(fp.matrix);
fp.nestedness.Detect();
nodf = fp.nestedness.N;

% Modularity
fp.community = LeadingEigenvector(fp.matrix);  % Newman
fp.community.DoKernighanLinTunning = true;  
fp.community.Detect();
qb = fp.community.Qb;
Nm = fp.community.N;




