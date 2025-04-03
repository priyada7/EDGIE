function [RvalueDetached,RvalueAttached,AreaDetached,AreaAttached]= Rcalc(UvalueWall,UvalueWindow,AreaDetached,AreaAttached,n1)

oneStoryHeight = 3; 
aspectRatio =1;
numberOfStories=2;
AreaDetached = trirnd(.9*AreaDetached,1.1*AreaDetached,n1,1);
AreaAttached = trirnd(.9*AreaAttached,1.1*AreaAttached,n1,1);

AwDetached = 2*oneStoryHeight*(aspectRatio + 1)*sqrt(numberOfStories*AreaDetached/aspectRatio);
AwAttached= 2*oneStoryHeight*(aspectRatio + 1)*sqrt(numberOfStories*AreaAttached/(aspectRatio))/2;

lambda = trirnd(0.2,0.3,n1,1); 


AreaRoofDetached = AreaDetached/numberOfStories;
AreaRoofAttached = AreaAttached/numberOfStories;
Ur = repmat(0.36,n1,1);

density = 1.293; %kg/m3
Cp = 1005; % J/kgK
volume = AreaDetached*oneStoryHeight; %m3
v = trirnd(0.2,1.5,n1,1) ; % air change per hour

mdotCp = v.*Cp*density.*volume/3600;
Uwall = (lambda.*trirnd(.8*UvalueWindow,1.2*UvalueWindow,n1,1) + (1-lambda).*trirnd(.8*UvalueWall,1.2*UvalueWall,n1,1));





RvalueDetached = 1./(mdotCp./1000 + Uwall.*AwDetached/1000 + Ur.*AreaRoofDetached/1000) ;
RvalueAttached = 1./(mdotCp./1000 + Uwall.*AwAttached/1000 + Ur.*AreaRoofAttached/1000) ;

 end
