function [RvalueDetached,RvalueAttached,AreaDetached,AreaAttached]= Rcalc(UvalueWall,UvalueWindow,AreaDetached,AreaAttached,n1)

oneStoryHeight = 3; 
aspectRatio =1;
numberOfStories=2;
%Area = repmat(200,length(UvalueWindow),1);
AreaDetached = trirnd(.9*AreaDetached,1.1*AreaDetached,n1,1);
AreaAttached = trirnd(.9*AreaAttached,1.1*AreaAttached,n1,1);

AwDetached = 2*oneStoryHeight*(aspectRatio + 1)*sqrt(numberOfStories*AreaDetached/aspectRatio);
AwAttached= 2*oneStoryHeight*(aspectRatio + 1)*sqrt(numberOfStories*AreaAttached/(aspectRatio))/2;

lambda = trirnd(0.2,0.3,n1,1); %repmat(0.25,n1,1);


AreaRoofDetached = AreaDetached/numberOfStories;
AreaRoofAttached = AreaAttached/numberOfStories;
Ur = repmat(0.36,n1,1);

density = 1.293; %kg/m3
Cp = 1005; % J/kgK
volume = AreaDetached*oneStoryHeight; %m3
v = trirnd(1,2,n1,1) ;%1.5; % air change per hour

mdotCp = v.*Cp*density.*volume/3600;
%Uwall = (lambda.*UvalueWindow + (1-lambda).*UvalueWall);
Uwall = (lambda.*trirnd(.8*UvalueWindow,1.2*UvalueWindow,n1,1) + (1-lambda).*trirnd(.8*UvalueWall,1.2*UvalueWall,n1,1));

h_outer = repmat(22.7,n1,1)/1000;
h_inner = repmat(8.29,n1,1)/1000;



RvalueDetached = 1./(mdotCp./1000 + Uwall.*AwDetached/1000 + Ur.*AreaRoofDetached/1000) + 1./(h_inner.*AwDetached) + 1./(h_outer.*AwDetached);
RvalueAttached = 1./(mdotCp./1000 + Uwall.*AwAttached/1000 + Ur.*AreaRoofAttached/1000) + 1./(h_inner.*AwAttached) + 1./(h_outer.*AwAttached);

 end