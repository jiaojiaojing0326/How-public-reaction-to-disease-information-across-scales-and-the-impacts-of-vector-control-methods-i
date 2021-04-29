% ==========================================================================================================================================================================================================================================

% Matlab code to generate all csv files under strength only situation in the manuscript titled "How public reaction to disease information across scales and the impacts of vector control methods influence disease prevalence and control efficacy"

% =========================================================================================================================================================================================================================================
 
function [S_out, I_out, R_out]  = Multi_Patch_StrengthOnly(folder, information, MaxTime, p, alpha, gamma, epsilon, N, seed, network)
% This function calculates and plots the temporal evolution of the model 
% The output is the values of human succeptible population and
% recovered population at the end of the simulation, 
% and  the peak of infected people for every patch, as a function of the
% environmental concern parameter (epsilon). 

% All the constants that I will not change
beta    = [0.15 0.30]/1000;       % transmision rate for [humans mosquitoes]
mu      = [2.35616e-05  1/13];    % natural death rate for [humans mosquitoes]
b       = 2.4657534e-05;          % human birth rate
r       = 0.037;                  % human recovery rate
w       = 0;                      % human death rate from disease
delta   = 5.468913e-05;           % composite rate 
nu      = 1/7;                    % maturation rate
eta     = 10;                     % egg laying rate - it says 10 on the table
m       = 10;                     % number of equations


if nargin == 0                          % if the number of input arguments is zero
    seed     = randi(intmax('int32'),1,1);
    MaxTime  = 300;                     % end of simulation
    epsilon  = [100 100] ;              % environment concern demotivates [mosquitos  larvae] pesticide
    gamma    = exp(-epsilon/80);        % motivation to control [mosquitos  larvae] due to Infected humans
    alpha    = 100 * gamma;             % motivation to control [mosquitos  larvae] due to Severe Outcomes (always 100 times larger than the regular infected)
    N        = 9;                      % number of patches
    p        = 0.01;                     % fraction of people that moves to other patches; this parameter determines the disease spreading from one location to others
    network  = 'SF';                    % type of network: SQ, SF, full
    information = 'local';             % 'global' or 'local' information
    folder = join([network, information, 'N', string(N),'p', string(p), 't', string(MaxTime),'K=800','strengthonly'], {'-'} );
    mkdir(folder)
end

rng(seed) % seed for pseudorandom number generator

filename = strcat(folder, '/2-a1-', '-a2-','-g1-', '-g2-', '-e1-', '-e2-');

for q = 1:50
  
%initial conditions in each patch
X0 = [];
Kl = [];
%delay = [];
for i = 1:N
       %  [R0 appears two times, the first one is for recovered people as a
       %  function of time in a given patch, the second one is reserved to
       %  accumulate the total number of recovered people since the begining
       %  [S0 Humans              ,I0, R0, D0, L0,              S0 mosq.,I0, Cm0, Cl0, R0]
    ini = [700 + randi(10,1,1);  0;  0;  0;  0; 1000 + randi(10,1,1); 0;   0;   0;   0];
    X0  = [X0; ini];
    Kl  = [Kl; 800 + randi(10,1,1)];  %carrying capacity of each patch: 500, 800 and 2000.
end


%%------------ Network and mobility --------------------------------------
% define the network structure and mobility conditions for different cases

if strcmp(network, 'SF')
    % Scale Free network: Molloy Reeds algorithm40
    mov = scale_free(N, 5, 4); %(total nodes, number of initial nodes, number of nodes a new added node is connected to)
    dd = distance(mov);
    %choose one random patch and infect one person
    J0 = randi(N-1,1,1);
end

X0((J0-1)*m+2) = 1;  % J0 is the first infected patch
mov1 = mov;
  
for jj=1:N
    mov(jj,:)= mov(jj,:)/sum(mov(jj,:));  %normalized by rows. (people leaving can't be more than 100%)
end 
mov = p * mov;    % p is the fraction of people that leave the patch
%%------------------------ END --------------------------------------------
                            
tspan = 0:MaxTime/600:MaxTime;  %500 points from 0 to MaxTime

%solving the equations
if strcmp(information, 'local')  
    options = odeset('RelTol',1e-3,'AbsTol',1e-3);
    [timeaxis,pop] = ode45(@SIR_eq_local, tspan, X0, options);
end

if strcmp(information, 'neighbor')
    options = odeset('RelTol',1e-3,'AbsTol',1e-3);
    [timeaxis,pop] = ode45( @SIR_eq_neighbor, tspan, X0, options);
end

if strcmp(information, 'global')
   options = odeset('RelTol',1e-3,'AbsTol',1e-3);
   [timeaxis,pop] = ode45(@SIR_eq_global, tspan, X0, options);  
end

% separate the results in different variables
Sh=pop(:,1:m:m*N);
Ih=pop(:,2:m:m*N);
Rh=pop(:,10:m:m*N);
Dh=pop(:,4:m:m*N);
Lm=pop(:,5:m:m*N);
Sm=pop(:,6:m:m*N);
Im=pop(:,7:m:m*N);
Cm=pop(:,8:m:m*N);
Cl=pop(:,9:m:m*N);

% % Normalize by the total number of people on each patch
N_total = zeros(length(Sh(:,1)),1);
M_total = zeros(length(Sh(:,1)),1);
for j = 1:N
    N_total(:) = Sh(:,j) + Ih(:,j) + Rh(:,j);
    M_total(:) = Sm(:,j) + Im(:,j);
   
    Sh(:,j) = Sh(:,j) ./ N_total;
    Ih(:,j) = Ih(:,j) ./ N_total;
    Rh(:,j) = Rh(:,j) ./ N_total;
    Dh(:,j) = Dh(:,j) ./ N_total;
    Lm(:,j) = Lm(:,j) ./ Kl(j);
    Sm(:,j) = Sm(:,j) ./ M_total;
    Im(:,j) = Im(:,j) ./ M_total;
end

% save all the variables as function of time in a text file for ploting later 
% timeaxis = pop.x.';
% save(strcat(filename, '-all_vs_t'),"timeaxis","*h","*m","Cl")

% plot all the variables as function of time and save them in a PDF file
%plot_all(strcat(filename,'.pdf'),pop.x,Sh,Ih,Rh,Dh,Lm,Sm,Im,Cm,Cl,N);

% S_out and R_out only save the last value of each of each patch
% I_out saves the maximum value of infected individuals
% these variables are the return value of this function
S_out = Sh(end,:) ;
I_out = Ih(end,:) ;
R_out = Rh(end,:) ;

%TEMPORAL SERIES, 9 FILES ONE FOR EACH VARIABLE, ONE COLUMN FOR EACH PATCH
dlmwrite(strcat(filename, '-SUSCEPTIBLE-HUMANS.csv'), [timeaxis Sh], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-INFECTED-HUMANS.csv'), [timeaxis Ih], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-RECOVERED-HUMANS.csv'), [timeaxis Rh], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-SEVERE-CASES-D.csv'), [timeaxis Dh], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-MOSQUITOES-LARVAE.csv'), [timeaxis Lm], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-SUSCEPTIBLE-MOSQUITOES.csv'), [timeaxis Sm], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-INFECTED-MOSQUITOES.csv'), [timeaxis Im], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-PESTICIDES-ADULTS.csv'), [timeaxis Cm], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-PESTICICDES-LARVAE.csv'), [timeaxis Cl], 'delimiter',',','-append')
dlmwrite(strcat(filename, '-DISTANCE-PATCH.csv'), dd(J0,:), 'delimiter',',','-append')


end


function dX = SIR_eq_local(t, X)

    %m = 10; %% number of equations for each patch
    
    % References:
    % X(1+m*j) = susceptible humans
    % X(2+m*j) = infected humans
    % X(3+m*j) = recovered humans
    % X(4+m*j) = Severe Outcomes
    % X(5+m*j) = mosquito larvae
    % X(6+m*j) = susceptible mosquitoes
    % X(7+m*j) = infected mosquitoes
    % X(8+m*j) = mosquito control
    % X(9+m*j) = larvae control
    % X(10+m*j) = total number of recovered people for 
    
    %X = pop(1:m*N);
    dX = zeros(m*N,1);

    for j = 0:N-1 
              
        if  X(8+j*m) <0
            X(8+j*m) = 0;
        end
        
        if  X(9+j*m) <0
            X(9+j*m) = 0;
        end

        flux = flux_of_people(X(1:m:end), j+1, mov, N);
        dX(1+j*m) = b * (X(1+m*j) + X(2+m*j) + X(3+m*j)) - beta(1) * X(7+m*j) * X(1+m*j) - mu(1) * X(1+m*j) + flux; %=dS^h/dt
        
        flux = flux_of_people(X(2:m:end), j+1, mov, N);
        dX(2+j*m) = beta(1) * X(7+m*j) * X(1+m*j) + X(2+m*j) * ( -r - mu(1) - w ) + flux;%=dI^h/dt
        
        flux = flux_of_people(X(3:m:end), j+1, mov, N);
        dX(3+j*m) = r * X(2+m*j) - mu(1) * X(3+m*j) + flux;%=dR^h/dt
        
        dX(4+j*m) = delta * X(2+m*j);%=dD^h/dt

        tmp = F(eta * (X(6+m*j) + X(7+m*j)), Kl(j+1)) - nu * X(5+m*j) - X(9+m*j) * X(5+m*j) ;  % = dL^m/dt    
        dX(5+j*m) = check_neg(X(5+m*j), tmp);
         
        tmp = nu * X(5+m*j) - beta(2) *  X(2+m*j) * X(6+m*j) - mu(2) * X(6+m*j) - X(8+m*j)* X(6+m*j); %=dS^m/dt
        dX(6+j*m) = check_neg(X(6+m*j), tmp);
        
        tmp = beta(2) * X(2+m*j) * X(6+m*j) - mu(2) * X(7+m*j) - X(8+m*j)*X(7+m*j); %=dI^m/dt 
        dX(7+j*m) = check_neg(X(7+m*j), tmp);
        
        tmp = alpha(1) * (X(4+m*j)) * binary(X(4+m*j)) + gamma(1) * (X(2+m*j)) * binary(X(2+m*j)) - (epsilon(1) * X(8+m*j)); % mosquito control
        dX(8+j*m) = check_neg_or_one(X(8+m*j), tmp);
        dX(8+j*m) = check_neg(X(8+m*j), tmp);
        
        tmp = (alpha(2) * (X(4+m*j)) * binary(X(4+m*j)) + gamma(2) * (X(2+m*j)) * binary(X(2+m*j)) - (epsilon(2) * X(9+m*j)))*(1- X(8 + j*m)/(0.00001 + X(8 + j*m) + X(9 + j*m))); % larvae control
        dX(9+j*m) = check_neg_or_one(X(9+m*j), tmp);
        dX(9+j*m) = check_neg(X(9+m*j), tmp);
        
        dX(10+j*m) = r * X(2+m*j);
        
        if  X(8+j*m) <0
            X(8+j*m) = 0;
        end
        
        if  X(9+j*m) <0
            X(9+j*m) = 0;
        end
    end
end

function dX = SIR_eq_neighbor(t, X)

    %m = 10; %% number of equations for each patch
    
    % References:
    % X(1+m*j) = susceptible humans
    % X(2+m*j) = infected humans
    % X(3+m*j) = recovered humans
    % X(4+m*j) = Severe Outcomes
    % X(5+m*j) = mosquito larvae
    % X(6+m*j) = susceptible mosquitoes
    % X(7+m*j) = infected mosquitoes
    % X(8+m*j) = mosquito control
    % X(9+m*j) = larvae control
    % X(10+m*j) = total number of recovered people for 
    
    %X = pop(1:m*N);
    dX = zeros(m*N,1);

    for j = 0:N-1 
        
        if  X(8+j*m) <0
            X(8+j*m) = 0;
        end
        
        if  X(9+j*m) <0
            X(9+j*m) = 0;
        end
        
        Ih_neighbor = adj_function(j, X, 2); %sum(X(2:m:end));
        Dh_neighbor = adj_function(j, X, 4);%sum(X(4:m:end));
                
        flux = flux_of_people(X(1:m:end), j+1, mov, N);
        dX(1+j*m) = b * (X(1+m*j) + X(2+m*j) + X(3+m*j)) - beta(1) * X(7+m*j) * X(1+m*j) - mu(1) * X(1+m*j) + flux; %=dS^h/dt
        
        flux = flux_of_people(X(2:m:end), j+1, mov, N);
        dX(2+j*m) = beta(1) * X(7+m*j) * X(1+m*j) + X(2+m*j) * ( -r - mu(1) - w ) + flux;%=dI^h/dt
        
        flux = flux_of_people(X(3:m:end), j+1, mov, N);
        dX(3+j*m) = r * X(2+m*j) - mu(1) * X(3+m*j) + flux;%=dR^h/dt
        
        dX(4+j*m) = delta * X(2+m*j);%=dD^h/dt

        tmp = F(eta * (X(6+m*j) + X(7+m*j)), Kl(j+1)) - nu * X(5+m*j) - X(9+m*j) * X(5+m*j) ;  % = dL^m/dt    
        dX(5+j*m) = check_neg(X(5+m*j), tmp);
         
        tmp = nu * X(5+m*j) - beta(2) *  X(2+m*j) * X(6+m*j) - mu(2) * X(6+m*j) - X(8+m*j)* X(6+m*j); %=dS^m/dt
        dX(6+j*m) = check_neg(X(6+m*j), tmp);
        
        tmp = beta(2) * X(2+m*j) * X(6+m*j) - mu(2) * X(7+m*j) - X(8+m*j)*X(7+m*j); %=dI^m/dt 
        dX(7+j*m) = check_neg(X(7+m*j), tmp);
        
        tmp = alpha(1) * (Dh_neighbor + X(4+m*j))* binary(X(4+m*j))+ gamma(1) * (Ih_neighbor + X(2+m*j))* binary(X(2+m*j)) - (epsilon(1) * X(8+m*j)); % mosquito control
        dX(8+j*m) = check_neg_or_one(X(8+m*j), tmp);
        dX(8+j*m) = check_neg(X(8+m*j), tmp);
        
        tmp = (alpha(2) * (Dh_neighbor + X(4+m*j))* binary(X(4+m*j)) + gamma(2) * (Ih_neighbor + X(2+m*j))* binary(X(2+m*j)) - (epsilon(2) * X(9+m*j)))*(1- X(8 + j*m)/(0.00001 + X(8 + j*m) + X(9 + j*m))); % larvae control
        dX(9+j*m) = check_neg_or_one(X(9+m*j), tmp);
        dX(9+j*m) = check_neg(X(9+m*j), tmp);

        dX(10+j*m) = r * X(2+m*j);
        
        if  X(8+j*m) <0
            X(8+j*m) = 0;
        end
        
        if  X(9+j*m) <0
            X(9+j*m) = 0;
        end
    end
end

function neighborcase = adj_function(jj,XX,kk)
  neighborcase = 0;
  for ii = 0: N-1
  neighborcase = neighborcase + mov1(jj+1,ii+1) * XX(kk+m*ii);  
  end

end


function dX = SIR_eq_global(t, X)

    %m = 10; %% number of equations for each patch
    
    % References:
    % X(1+m*j) = susceptible humans
    % X(2+m*j) = infected humans
    % X(3+m*j) = recovered humans
    % X(4+m*j) = Severe Outcomes
    % X(5+m*j) = mosquito larvae
    % X(6+m*j) = susceptible mosquitoes
    % X(7+m*j) = infected mosquitoes
    % X(8+m*j) = mosquito control
    % X(9+m*j) = larvae control
    % X(10+m*j) = total number of recovered people for 
    
    %X = pop(1:m*N);
    dX = zeros(m*N,1);


    for j = 0:N-1 
        
        if  X(8+j*m) <0
            X(8+j*m) = 0;
        end
        
        if  X(9+j*m) <0
            X(9+j*m) = 0;
        end
        
        Ih_tot = sum(X(2:m:end));
        Dh_tot = sum(X(4:m:end)); 
        
        flux = flux_of_people(X(1:m:end), j+1, mov, N);
        dX(1+j*m) = b * (X(1+m*j) + X(2+m*j) + X(3+m*j)) - beta(1) * X(7+m*j) * X(1+m*j) - mu(1) * X(1+m*j) + flux; %=dS^h/dt
        
        flux = flux_of_people(X(2:m:end), j+1, mov, N);
        dX(2+j*m) = beta(1) * X(7+m*j) * X(1+m*j) + X(2+m*j) * ( -r - mu(1) - w ) + flux;%=dI^h/dt
        
        flux = flux_of_people(X(3:m:end), j+1, mov, N);
        dX(3+j*m) = r * X(2+m*j) - mu(1) * X(3+m*j) + flux;%=dR^h/dt
        
        dX(4+j*m) = delta * X(2+m*j);%=dD^h/dt

        tmp = F(eta * (X(6+m*j) + X(7+m*j)), Kl(j+1)) - nu * X(5+m*j) - X(9+m*j) * X(5+m*j) ;  % = dL^m/dt    
        dX(5+j*m) = check_neg(X(5+m*j), tmp);
         
        tmp = nu * X(5+m*j) - beta(2) *  X(2+m*j) * X(6+m*j) - mu(2) * X(6+m*j) - X(8+m*j)* X(6+m*j); %=dS^m/dt
        dX(6+j*m) = check_neg(X(6+m*j), tmp);
        
        tmp = beta(2) * X(2+m*j) * X(6+m*j) - mu(2) * X(7+m*j) - X(8+m*j)*X(7+m*j); %=dI^m/dt 
        dX(7+j*m) = check_neg(X(7+m*j), tmp);
        
        tmp = alpha(1) * Dh_tot * binary(X(4+m*j)) + gamma(1) * Ih_tot * binary(X(2+m*j))- (epsilon(1) * X(8+m*j)) ; % mosquito control
        dX(8+j*m) = check_neg_or_one(X(8+m*j), tmp);
        dX(8+j*m) = check_neg(X(8+m*j), tmp);
        
        tmp = (alpha(2) * Dh_tot * binary(X(4+m*j)) + gamma(2) * Ih_tot * binary(X(2+m*j))- (epsilon(2) * X(9+m*j)))*(1- X(8 + j*m)/(0.00001 + X(8 + j*m) + X(9 + j*m))); % larvae control
        dX(9+j*m) = check_neg_or_one(X(9+m*j), tmp);
        dX(9+j*m) = check_neg(X(9+m*j), tmp);

        dX(10+j*m) = r * X(2+m*j);
        
        if  X(8+j*m) < 0
            X(8+j*m) = 0;
        end
        
        if  X(9+j*m) < 0
            X(9+j*m) = 0;
        end
    end
end

end 

function new = check_neg_or_one(old_value, change)
% check if a variable is about to be larger than one
    if  old_value + change > 1
        new = 0;
        return 
    else
        new = change;
        return 
    end
% check if a variable is about to be lower than zero
    if  old_value + change < 0
        new = 0;
        return 
    else
        new = change;
        return 
    end
end 

function new = check_neg(old_value, change)
% check if a variable is about to be lower than zero
    if  old_value + change < 0
        new = 0;
    else
        new = change;
    end
end 

function cc = F(x, Kl)
% function for carrying capacity in the eq. for larvaes:
    cc = x * ( 1 - ( x / Kl ));
    %if cc < 0
        %cc = 0;  
    %end
end

function flux = flux_of_people(X, current_patch, mov, N)
% calculate the term for the in and out flux of people in a patch
% mov(i,j) is the fraction of people that moves from i to j

    j = current_patch;

    in_flux_tmp = mov * X;
    in_flux = in_flux_tmp(j);
    
    out_flux = sum(mov(j,:)) * X(j) ;
    
    flux = in_flux - out_flux ;
end

function h = plot_all(txt,t,X1,X2,X3,X4,X5,X6,X7,X8,X9,N)
fig=figure('Name','SIR', 'units','inch','position',[0,0,14,8]) ;

fig.PaperType = 'uslegal';
fig.PaperOrientation = 'landscape';

subplot(3,3,1);
h=plot(t,X1(:,1:N));
xlabel 'Time';
ylabel 'Susceptible Humans';

subplot(3,3,2);
h=plot(t,X2(:,1:N));
xlabel 'Time';
ylabel 'Infected Humans';

title(txt)

subplot(3,3,3);
h=plot(t,X3(:,1:N));
xlabel 'Time';
ylabel 'Recovered Humans';

subplot(3,3,4) ;
h=plot(t,X4(:,1:N));
xlabel 'Time';
ylabel 'Severe Outcomes';

subplot(3,3,5) ;
h=plot(t,X8(:,1:N));
xlabel 'Time';
ylabel 'Mosquitoes Control';

subplot(3,3,6) ;
h=plot(t,X9(:,1:N));
xlabel 'Time';
ylabel 'Larvae Control';

subplot(3,3,7) ;
h=plot(t,X5(:,1:N));
xlabel 'Time';
ylabel 'Larvae';

subplot(3,3,8) ;
h=plot(t,X6(:,1:N));
xlabel 'Time';
ylabel 'Susceptible Mosquitoes';

subplot(3,3,9) ;
h=plot(t,X7(:,1:N));
xlabel 'Time';
ylabel 'Infected Mosquitoes';

%print(txt,'-depsc')
saveas(gcf, txt);  %txt is the path and filename 

end

function D=distance(A)
% calculate de distance matrix from a given adjancency matrix A

N = max(size(A));
D = NaN(N);
B = A;
k = 1;
%while any(isnan(D(:)))
while k <= N
    % Check for new walks, and assign distance
    D(B > 0 & isnan(D)) = k;
    % Iteration
    k = k + 1;
    B = B * A;
end
D(1:N+1:end) = 0; % diagonal elements are always 0
% Now D contains the distance matrix
end 

function B = binary (value)

    if value > 1.000e-323
       B = 1;
       return
    else
       B = 0;
       return
    end
end


