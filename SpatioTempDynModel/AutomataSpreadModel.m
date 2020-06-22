% -----------------------------------------------------------------
% Cellular Automata (CA) algorithm by Munshi et. al.              %
% Programmed by Joydeep Munshi and Indranil Roy at Lehigh 	  	  %
% University
% Programming dates: Mar 2020 - May 2020                          %
% Last revised: May  2020    								      %
% -----------------------------------------------------------------
% Papers -- Citation Details:
% J. Munshi, I. Roy, G. Balasubramanian, Spatiotemporal           %
% dynamics in demography-sensitive disease transmission:          %
% COVID-19 spread in NY as a case study.						  %
% https://arxiv.org/ftp/arxiv/papers/2005/2005.01001.pdf
% ----------------------------------------------------------------%
% This program implements a version of Cellular Automata framework%
% (CA) to investigate spatio-temporal dynamics of disease outbreak% 
% sensitive to demographic features such as population density,   %
% mobility and employment status. The code implement here uses    %
% Matlab's vectorization and mapping toolbox.                     %

% --------------------------------------------------------------- %
% =============================================================== %
% Notes:                                                          %
% Different implementations may lead to slightly different        %
% behavour and/or results, but there is nothing wrong with it,    %
% as this is the nature of random walks and all metaheuristics.   %
% -----------------------------------------------------------------


function [TotalSuc,TotalLat,TotalQuar,TotalInfect,TotalIso,TotalRecov,TotalDeath] = AutomataSpreadModel(lambda,t1,t2,t3,t4,timePeriod,maxMov,radius)
% --------------------------------------------------------------- %
%Author : Joydeep Munshi, Indranil Roy							  %
%Cellular Automata to model Covid 19 spread						  %
%   Detailed explanation goes here								  %
% S(t), L(t), Q(t), I(t), J(t), R(t), D(t)                        %
% --------------------------------------------------------------- %
% Description of arguments:										  %
% lambda      - Probability of death, morbidity                   %
% t1          - Incubation/Latency period		      			  %
% t2          - Immunization/Recovery period					  %
% t3          - Loss of immunity period							  %
% t4     	  - Lockdown period									  %
% timePeriod  - Simulation time									  %
% maxMov      - Maximum randomWalk								  %
% radius      - Initial radius of infection						  %
% --------------------------------------------------------------- %
if nargin<8
    disp('Please provide 8 arguments');
    return
end

%Initialization
N = 500;						% Number of initial cells NxN
alpha = 0.0001;                 % percentage of quarentined getting infected (function of age group and health)
beta = 0.75;                    % percentage of infected to isolation
%gamma = 0.05;                  % percentage of people going to quarentnine (function of time)

TotalSuc  =  zeros(timePeriod,1);
TotalLat =  zeros(timePeriod,1);
TotalQuar =  zeros(timePeriod,1);
TotalInfect =  zeros(timePeriod,1);
TotalIso =  zeros(timePeriod,1);
TotalRecov  =  zeros(timePeriod,1);
TotalDeath  =  zeros(timePeriod,1);

[population,C] = NYSpopulation(N);

T1 = zeros(N);
T2 = zeros(N);
T3 = zeros(N);
T4 = zeros(N);

%Initiate all cells with susceptible
f = zeros(N);
P = zeros(N);
Tc = rand(N);
C = ones(N);

[index1,index2]= find(population==1);
if radius ~= 0
    for in1 = -radius:radius
        for in2 = -radius:radius
            if (index1(1)+in1 >= 1 && index1(1)+in1 <= N) && (index2(1)+in2 >= 1 && index2(1)+in2 <= N)
                %disp(index1(1)+in1)
                %disp(index2(1)+in2)
                if C(index1(1)+in1,index2(1)+in2) ~= -10
                    C(index1(1)+in1,index2(1)+in2) = 4;    %Perturb with 1 infection
                end
            end
        end
    end
else
    C(index1(1),index2(1)) = 4;
end

T2(find(C==4)) = 1;
f(find(C==4)) = rand;
P=sqrt(f);

TotalSuc(1,1)  =  sum(sum(C==1));
TotalLat(1,1) =  sum(sum(C==2));
TotalQuar(1,1)  =  sum(sum(C==3));
TotalInfect(1,1)  =  sum(sum(C==4));
TotalIso(1,1)  =  sum(sum(C==5));
TotalRecov(1,1)  =  sum(sum(C==6));
TotalDeath(1,1)  =  sum(sum(C==0));

h = pcolor(C);
set(h, 'EdgeColor', 'none');
%shading interp
colormap jet
caxis([-10 6])

t=1;
while t < timePeriod
    if t==15													% Varies with lockdown period %
        gamma=0.3;
    else
        gamma=0;
    end
    
    for i=2:N-1
        for j=2:N-1
            if C(i,j) == 1
                dpGamma = cumsum([gamma 1-gamma]);                 
                if t==15
                    cGamma = [2 3];
                else
                    cGamma = [3 2];
                end
                c1Gamma = cGamma(1+sum(dpGamma(end)*rand>dpGamma));
                if c1Gamma==2
                    p =sqrt(1-getRc(Tc(i,j),population(i,j)))*( P(i,j) + P(i-1,j) + P(i,j-1) + P(i-1,j-1) + P(i+1,j) + P(i,j+1) + P(i+1,j+1) + P(i-1,j+1) + P(i+1,j-1));
                    if (p*0.25) >= 0.01  % How to determine the thresold
                        C(i,j) = 2;
                        %f(i,j) = rand;
                        T1(i,j) = T1(i,j)+1;
                    else
                        C(i,j) = 1;
                    end
                else
                    C(i,j) = 3;
                    T4(i,j) = T4(i,j)+1;
                end
            elseif C(i,j) == 2
                dpBeta = cumsum([1-beta beta]);                 
                cBeta = [4 5];
                c1Beta = cBeta(1+sum(dpBeta(end)*rand>dpBeta));
                if c1Beta == 4 && T1(i,j)> t1
                    C(i,j) = 4;
                    f(i,j) = rand;
                    P = sqrt(f);
                    T2(i,j) = T2(i,j)+1;
                elseif c1Beta ~= 4 && T1(i,j)> t1
                    C(i,j) = 5;
                    T2(i,j) = T2(i,j)+1;
                elseif T1(i,j)<= t1
                    T1(i,j) = T1(i,j)+1;
                end
            elseif C(i,j) == 3 
                dpAlpha = cumsum([alpha 1-alpha]);                 
                cAlpha = [5 3];
                c1Alpha = cAlpha(1+sum(dpAlpha(end)*rand>dpAlpha));
                if c1Alpha == 5
                    C(i,j) = 5;
                    T2(i,j) = T2(i,j)+1;
                    T4(i,j) = 0;
                elseif c1Alpha == 3 && T4(i,j) > t4
                    C(i,j) = 1;
                    T4(i,j) = 0;
                elseif c1Alpha == 3 && T4(i,j) <= t4
                    C(i,j) = 3;
                    T4(i,j) = T4(i,j)+1;
                end
            elseif C(i,j) == 4 || C(i,j) == 5
                nfrac = (sum(sum(C==4))+sum(sum(C==5)))/(N*N);
                dplambda = cumsum([getLambda(lambda,nfrac,population(i,j)) 1-getLambda(lambda,nfrac,population(i,j))]);   % death rate %
                clambda = [0 6];
                c1lambda = clambda(1+sum(dplambda(end)*rand>dplambda));
                if c1lambda == 6 && T2(i,j)>t2
                    T3(i,j) = T3(i,j) + 1;
                    T2(i,j) = 0;
                    C(i,j) = 6;
                elseif c1lambda == 0
                    C(i,j) = 0;
                    T2(i,j) = 0;
                elseif T2(i,j)<=t2 && c1lambda == 6
                    T2(i,j) = T2(i,j)+1;
                    C(i,j)=C(i,j);
                end
            elseif C(i,j) == 6
                if T3(i,j) > t3
                    C(i,j) = 1;
                    T3(i,j) = 0;
                else 
                    T3(i,j) = T3(i,j)+ 1;
                end
            else 
                C(i,j)=C(i,j);
            end
        end
    end
    
    [C,f,P,T1,T2,T3] = randomWalk(C,f,P,T1,T2,T3,N,maxMov,population);
    t=t+1;
    TotalSuc(t,1)  =  sum(sum(C==1));
    TotalLat(t,1) =  sum(sum(C==2));
    TotalQuar(t,1)  =  sum(sum(C==3));
    TotalInfect(t,1)  =  sum(sum(C==4));
    TotalIso(t,1)  =  sum(sum(C==5));
    TotalRecov(t,1)  =  sum(sum(C==6));
    TotalDeath(t,1)  =  sum(sum(C==0));
    
    fprintf('Day: %i \n', t);
    fprintf('Total Susc: %i \n', TotalSuc(t,1));
    fprintf('Total Infected: %i \n', TotalInfect(t,1));
    fprintf('Total Recovered: %i \n', TotalRecov(t,1));
    fprintf('Total Death: %i \n', TotalDeath(t,1));
    
    h = pcolor(C);
    caxis([-10 6]);
    set(h, 'EdgeColor', 'none');
    day = t-1;
    title(['This is Day: ', day]);
    M(t-1) = getframe(gcf);
    
end

%time = visualize(TotalSuc,TotInfect,TotRecov,TotDeath);
%[DeathSmooth,conf] = visualize(TotalSuc,TotalLat,TotalQuar,TotalInfect,TotalIso,TotalRecov,TotalDeath,timePeriod);

end

function [Rc] = getRc(Tc,fpop)
% Function to calculate resistivity based on influential factors  %
% --------------------------------------------------------------- % 
% Description of arguments:										  %
% Tc      - Uniform random distribution of transmission	          %
% fpop    - Population density matrix		         			  %
% --------------------------------------------------------------- %


% Change values according to state
age = [0.126 0.625 0.249];                    %age group classes 0-10,10-50,50-up
ageF = [0.2 0.6 0.2];
fage = age*ageF';

health = [0.065 0.310078 0.6249];             %respiratory and cardiac
healthF = [0.05 0.15 0.8];
fhealth = health*healthF';

employment = [0.0557 0.0403 0.903];           %employment [heatlth care, Emmergency service, Other]
employmentF= [0.1 0.2 0.7];
femployment=employment*employmentF';

Rc = fhealth*fage*femployment*Tc*(1-fpop);

end

function [lambda] = getLambda(l,n,fpop)

% Change values according to state
age = [0.126 0.625 0.249];                    %age group classes 0-10,10-50,50-up
ageF = [0.3 0.1 0.6];
fage = age*ageF';

health = [0.065 0.310078 0.6249];             %respiratory and cardiac
healthF = [0.85 0.1 0.05];
fhealth = health*healthF';

employment = [0.0557 0.0403 0.903];           %employment [heatlth care, Emmergency service, Other]
employmentF= [0.1 0.2 0.7];
femployment=employment*employmentF';

lambda = l*exp(n)*fage*fhealth*femployment;

end


function [C,f,P,T1,T2,T3] = randomWalk(C,f,P,T1,T2,T3,N,maxMov,p)
% Function to calculate random moves of Cellular automata         %
% --------------------------------------------------------------- % 
% Description of arguments:										  %
% C          - Cell state matrix                                  %
% f          - Uniform distribution matrix		      			  %
% T1         - Latency period matrix			     			  %
% T2         - Immunization/Recovery period matrix				  %
% T3     	 - Loss of immunity period matrix					  %
% N          - Simulation time									  %
% maxMov     - Maximum randomWalk								  %
% p          - Initial radius of infection						  %
% --------------------------------------------------------------- %

% Global random move probability
globalP=0.001;													  

for i=1:N
    for j=1:N
        dplocal = cumsum([globalP 1-globalP]);
        local = [1 0];
        loc = local(1+sum(dplocal(end)*rand>dplocal));
        if loc == 0
            shuff=[randi([-maxMov maxMov],1,1),randi([-maxMov maxMov],1,1)];
            shuff(1) = round(shuff(1)*p(i,j));                   
            shuff(2) = round(shuff(2)*p(i,j));
        else
            shuff=[randi([-maxMov maxMov],1,1),randi([-maxMov maxMov],1,1)];
        end
        
        if (i+shuff(1)) <= 1 
            inew=2;
        elseif (i+shuff(1)) >= N
            inew=N-1;
        else
            inew = i + shuff(1);
        end
        
        if (j+shuff(2)) <= 1 
            jnew=2;
        elseif (j+shuff(2)) >= N
            jnew=N-1;
        else
            jnew = j + shuff(2);
        end
        if (C(i,j) ~= 0 && C(i,j) ~= 3 && C(i,j) ~= 5 && C(i,j) ~= -10 )
            if(C(inew,jnew) ~= 0 && C(inew,jnew) ~= 3 && C(inew,jnew) ~= 5 && C(inew,jnew) ~= -10)
                index1=sub2ind([N N],[i],[j]);
                index2=sub2ind([N N],[inew],[jnew]);
                C([index1 index2]) = C([index2 index1]);
                f([index1 index2]) = f([index2 index1]);
                P([index1 index2]) = P([index2 index1]);
                T1([index1 index2]) = T1([index2 index1]);
                T2([index1 index2]) = T2([index2 index1]);
                T3([index1 index2]) = T3([index2 index1]);
            else
                C(i,j) = C(i,j);
                f(i,j) = f(i,j);
                P(i,j) = P(i,j);
                T1(i,j) = T1(i,j);
                T2(i,j) = T2(i,j);
                T3(i,j) = T3(i,j);
            end
        else 
            C(i,j) = C(i,j);
            f(i,j) = f(i,j);
            P(i,j) = P(i,j);
            T1(i,j) = T1(i,j);
            T2(i,j) = T2(i,j);
            T3(i,j) = T3(i,j);
        end
    end
end

end

function [DeathSmooth,conf] = visualize(TotalSuc,TotalLat,TotalQuar,TotalInfect,TotalIso,TotalRecov,TotalDeath,timeperiod)
% Function to visualize output %

figure(2)
plot(TotalSuc,'LineWidth',4)
hold on;
plot(TotalLat,'LineWidth',4);
hold on;
plot(TotalQuar,'LineWidth',4);
hold on;
plot(TotalInfect+TotalIso,'LineWidth',4)
hold on;
%plot(TotalIso,'LineWidth',4)
%hold on;
plot(TotalRecov,'LineWidth',4)
hold on;
%plot(TotalDeath,'LineWidth',4)
xlim([0 timeperiod])
ylim([0 100000])
legend('Susciptible','Latent','Quarentine','Infected','Recovered','Death')
graph= size(TotalSuc);

figure(3)

for i=1:timeperiod
    if i==1
        Death(i) = TotalDeath(i);
    else
        Death(i) = abs(TotalDeath(i)-TotalDeath(i-1));
    end
end

bar(Death);
hold on;
DeathSmooth = smooth(Death);
plot(DeathSmooth);


figure(4)
for i=1:timeperiod
    if i==1
        totalconfirmed(i) = TotalInfect(i)+TotalIso(i)+TotalRecov(i)+TotalDeath(i);
    else
        totalconfirmed(i) = abs(TotalInfect(i)+TotalIso(i)+TotalRecov(i)+TotalDeath(i)-(TotalInfect(i-1)+TotalIso(i-1)+TotalRecov(i-1)+TotalDeath(i-1)));
    end
end

conf = smooth(totalconfirmed);

plot(conf,'LineWidth',4)

end

