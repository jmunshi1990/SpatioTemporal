function [populationNorm,C] = PASpopulation(N)
% Function to calculate normalized population density and intial  %
% cell states                                                     %
% --------------------------------------------------------------- % 
% Description of arguments:										  %
% N          - Number of cell grids                               %
% --------------------------------------------------------------- %

X = linspace(PAS.BoundingBox(1,1),PAS.BoundingBox(2,1),N);
Y = linspace(PAS.BoundingBox(1,2),PAS.BoundingBox(2,2),N);
population = zeros(N);
[x,y]=meshgrid(X,Y);

for i=1:67
temp1 = PA(i).X;
temp2 = PA(i).Y;
name = PA(i).COUNTY_NAM;
    
a(i) = sum(temp1(:,1:(length(PA(i).X')-1)))/(length(PA(i).X')-1);
b(i) = sum(temp2(:,1:(length(PA(i).Y')-1)))/(length(PA(i).Y')-1);
[in, on] = inpolygon(x,y,PA(i).X,PA(i).Y);
onPop =  (PA(i).POP2020)/(PA(i).AREA_SQ_MI)* on;
inPop =  (PA(i).POP2020)/(PA(i).AREA_SQ_MI)* in;
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

