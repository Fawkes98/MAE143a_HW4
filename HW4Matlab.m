%% Problem 1a, HW #4
clear all; close all; clc;
% trying to solve equation of form " a*x + b*y = f " where we know a, b, & f
% insert givens in required form
a1 = RR_PolyConv([1 1], [1 -1], [1 3], [1 -3], [1 5], [1 -5]);
b1 = RR_PolyConv([1 2], [1 -2], [1 4], [1 -4]);
f1 = RR_PolyConv([1 1], [1 1], [1 3], [1 3], [1 5], [1 5])

%plug into Diophantine fxn
[x1a y1a] = RR_Diophantine(a1, b1, f1)

%Show residual
test1=RR_PolyAdd(RR_PolyConv(a1,x1a),RR_PolyConv(b1,y1a)); residual1=norm(RR_PolyAdd(f1,-test1))
fprintf('Note that the solution {x1,y1} is improper, but solves a*x1+b*y1=f1, with ~ zero residual\n\n')

poles1 = roots(x1a); zeros1 = roots(y1a); 
fprintf('The poles of the system D1(s) are at s =\n'); fprintf(' %3.3f\n', poles1(:,1));
fprintf('\nThe zeros of the system D1(s) are at s =\n'); fprintf(' %3.3f\n', zeros1(:,1));
fprintf('\n');

%% Problem 1b, HW#4
% spec new given
f2 = RR_PolyConv([1 1], [1 1], [1 3], [1 3], [1 5], [1 5], [1 50]);

%Diophantine 
[x1b y1b] = RR_Diophantine(a1, b1, f2); %first instance, improper

%while loop
p = 1; %count of p-order
    while length(x1b) <= length(y1b)
    f2 = RR_PolyConv(f2, [1 50]);
    p = p+1; %coordinated to numb of [1 50] terms multiplied into f2
    [x1b y1b] = RR_Diophantine(a1, b1, f2); %reevaluated
    end
[x1b y1b] = RR_Diophantine(a1, b1, f2) %filaized strictly proper system
fprintf('\nThe value of ''p'' to make D2(s) proper is p=%i\n\n', p);

poles1 = roots(x1a); zeros1 = roots(y1a); %calculate zeros and poles
fprintf('The poles of the system D2(s) are at s =\n'); fprintf(' %3.3f\n', poles1(:,1));
fprintf('\nThe zeros of the system D2(s) are at s =\n'); fprintf(' %3.3f\n', zeros1(:,1));

%% Problem 2, HW#4
clear all; clc;

%see PDF


%% Problem 3, HW#4
clear all; clc;

%set up symbolic funcitons, per RR pg9-5
d = .1 %value of d via Prob3
n = [1 2 4 8] %values of n via Prob3
syms k s
NUM3 = cell(length(n), 1); DEN3 = cell(length(n), 1); %placeholders, for prob 4
    for N = 1:length(n)
        %catagorize Pade fxns in seperable way for tf & rlocus
        ck = (factorial((2*n(N))-k)*factorial(n(N)))/(factorial(2*n(N))*factorial(k)*factorial(n(N)-k));
        FnNum = ((-1)^k)*ck*((d*s)^k); %numeratior as per Pade
        FnDen = ck*((d*s)^k); %denominator as per Pade
        %compute & plot variations of n
        Fn_N = symsum(FnNum, k, 0, n(N)); Fn_D = symsum(FnDen, k, 0, n(N)); %create summation expansions
        Nu = sym2poly(Fn_N); De = sym2poly(Fn_D); %grab coefficients for tf & root locus
        NUM3{N} = Nu; DEN3{N} = De; %stored for prob 4
        figure(N); hold on;
        rlocus(tf(Nu, De));
        title(['Root Locus plot of F_n(s), where n=', num2str(n(N))]);
    end

%% Problem 4a, HW#4
close all; clc;
figure(1); hold on;
rlocus(tf([1], [1 0]));
title('Root Locus plot of L_1(s) wrt K');

%% Problem 4b, HW#4
for i = 1:length(n) %for loop to multi s in DEN3 & plot
    DEN3{i}(1, end+1) = 0; %bump everything left one, essentially mulitplying s through denom
    figure(i+1); hold on;
    rlocus(tf(NUM3{i}, DEN3{i})); %plotting based on stored numerators & multi denoms
    title(['Root Locus plot of G_n(s), where n=', num2str(n(i))]);
end

%% Problem 4bb, HW#4
%follow same process as previous problems, truncate for efficiency
close all; clc;
%set up symbolic funcitons, per RR pg9-5
d = .2; %value of d via Prob3
% NUM4 = cell(length(n), 1); DEN4 = cell(length(n), 1); %placeholders, for next loop
    for N = 1:length(n)
        %catagorize Pade fxns in seperable way for tf & rlocus
        ck = (factorial((2*n(N))-k)*factorial(n(N)))/(factorial(2*n(N))*factorial(k)*factorial(n(N)-k));
        FnNum = ((-1)^k)*ck*((d*s)^k); %numeratior as per Pade
        FnDen = ck*((d*s)^k); %denominator as per Pade
        %compute & plot variations of n
        Fn_N = symsum(FnNum, k, 0, n(N)); Fn_D = symsum(FnDen, k, 0, n(N)); %create summation expansions
        Nu = sym2poly(Fn_N); De = sym2poly(Fn_D); %grab coefficients for tf & root locus
        De(end + 1) = 0; %symbolically 'multiply by s' by extension in denominator
        figure(N); hold on;
        rlocus(tf(Nu, De)); %plotting
        title(['Root Locus plot of G_n(s), with d = .2, where n=', num2str(n(i))], FontSize=10);
    end
















