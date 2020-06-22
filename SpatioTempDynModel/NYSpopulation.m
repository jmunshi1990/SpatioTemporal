function [populationNorm,C] = NYSpopulation(N)
% Function to calculate normalized population density and intial  %
% cell states                                                     %
% --------------------------------------------------------------- % 
% Description of arguments:										  %
% N          - Number of cell grids                               %
% --------------------------------------------------------------- %

NYS = shaperead('NYS_Civil_Boundaries.shp/NYS_Civil_Boundaries_SHP/Counties_Shoreline.shp');
NY = shaperead('NYS_Civil_Boundaries.shp/NYS_Civil_Boundaries_SHP/State_Shoreline.shp');

X = linspace(NY.BoundingBox(1,1),NY.BoundingBox(2,1),N);
Y = linspace(NY.BoundingBox(1,2),NY.BoundingBox(2,2),N);
population = zeros(N);
inAlbany = zeros(N);
inErie= zeros(N);
inMonroe= zeros(N);
inOnondoga= zeros(N);
[x,y]=meshgrid(X,Y);

for i=1:62
temp1 = NYS(i).X;
temp2 = NYS(i).Y;
name = NYS(i).NAME;
    
a(i) = sum(temp1(:,1:(length(NYS(i).X')-1)))/(length(NYS(i).X')-1);
b(i) = sum(temp2(:,1:(length(NYS(i).Y')-1)))/(length(NYS(i).Y')-1);
[in, on] = inpolygon(x,y,NYS(i).X,NYS(i).Y);
onPop =  (NYS(i).POP2010)/(NYS(i).CALC_SQ_MI)* on;
inPop =  (NYS(i).POP2010)/(NYS(i).CALC_SQ_MI)* in;
population = population + inPop +onPop;
temp1=0;
temp2=0;
inPop = 0;
onPop = 0;
in = 0;
on = 0;
end

populationNorm = population/max(max(population));

%[in, on] = inpolygon(x,y,NY.X,NY.Y);

C=ones(N);
C(populationNorm==0)=-10;
end

