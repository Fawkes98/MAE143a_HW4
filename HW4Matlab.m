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

%% Problem 1c, HW#4
