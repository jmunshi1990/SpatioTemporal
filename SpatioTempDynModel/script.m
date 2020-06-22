% ------------------------------------------------------------------------- %
% This script implements Cellular Automata function "AutomataSpreadModel.m" %
% Due to randomness of the CA framework, we recommend sampling over several %
% simulations. Please modify the script according to required arguments.    %
% ------------------------------------------------------------------------- %

% Number of random samples
C = 10;

for c=1:C
    fprintf('Started:  %i -th  samples \n', c);
    [S(:,c),L(:,c),Q(:,c),I(:,c),J(:,c),R(:,c),D(:,c)] = AutomataSpreadModel(0.25,5.1,20.2,90,45,200,150,0);
    fprintf('Finished:  %i -th  samples \n', c);
end

fprintf('Averaging over %i samples... \n', C);

% Sampling
TotalSucLD120(:,1) = sum(S')/C;
TotalLatLD120(:,1) = sum(L')/C;
TotalQuarLD120(:,1) = sum(Q')/C;
TotalInfectLD120(:,1) = sum(I')/C;
TotalIsoLD120(:,1) = sum(J')/C;
TotalRecovLD120(:,1) = sum(R')/C;
TotalDeathLD120(:,1) = sum(D')/C;
%clear s i r d r S I R D
clear c C S L Q I J R D;