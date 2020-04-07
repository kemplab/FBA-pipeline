%% Smith Decomposition Demo
% Author: Nadia Figueroa, Date: February 2014
% E-mail: nadia.figueroafernandez@epfl.ch
% Linear System Theory IC - Doctoral School
% Learning Algorithms and Systems Laboratory
% École polytechnique fédérale de Lausanne

%%%%%%%%%%%%%%%%%%%%%% Smith Decomposition Algorithm %%%%%%%%%%%%%%%%%%%%%%
% The Smith normal form (also called Smith Canonical form or Invariant Factor 
% theorem) is a diagonal matrix D that contains the invariant factors of any A 
% matrix of size n × m over a field F (in the attached implementation it is 
% provided for the ring of integers Z and rings of polynomials F[x]).
% By applying elementary row and column operations on A we can obtain the 
% following

% D = |d1 0 ... 0 ... 0|= TAS
%     |0 d2 ... 0 ... 0|
%     |:  : ... : ... :|
%     |0  0 ...dr ... 0|
%     |:  : ... : ... :|
%     |0  0 ... 0 ... 0|

% where d1 , ..., dr ∈ F are monic, dj |dj+1 for 1 ≤ k ≤ r − 1. T is a 
% product of elementary row unimodular matrices, and S is a product of 
% elementary column unimodular matrices.

% There are three types elementary row and column operations:
% 1. Multiply i-th row (column) of In by constant c.
% 2. Swap i-th and j-th row (column) from In .
% 3. Add α(s) times row j of In to row i.
% For row operations we define the unimodular matrices as T1 , T2 and T3 .

% T1(i,c)=|1 0 0 0 0|  T2(i,j)=|1 0 0 0 0| T3(i,j,a(s))=|1  0   0 0 0|
%         |0 1 0 0 0|          |0 0 0 1 0|              |0  1   0 0 0|
%         |0 0 c 0 0|          |0 0 1 0 0|              |0  0   1 0 0|
%         |0 0 0 1 0|          |0 1 0 0 0|              |0 a(s) 0 1 0|
%         |0 0 0 0 1|          |0 0 0 0 1|              |0  0   0 0 1|

% The results can be verified with the following Mathematica SmithForm implementation:
% http://library.wolfram.com/infocenter/MathSource/7081/


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Following some sample integer matrices whose smith form you can find online 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
clc

%Instructions: Uncomment A matrix

% Example 1
% Source: http://math.stackexchange.com/questions/133076/computing-the-smith-normal-form
% A = [ -6 111  -36 6;5 -672 210 74;0 -255 81 24;-7 255 -81 -10];
% D = |1 0 0  0|
%     |0 3 0  0|
%     |0 0 21 0|
%     |0 0 0  0|

% Example 2
% Source: http://en.wikipedia.org/wiki/Smith_normal_form
A = [ 2 4 4; -6 6 12; 10 -4 -16]; 
% D = |2 0 0 |  
%     |0 6 0 | 
%     |0 0 12|

% Example 3
% Source: http://www.mathworks.ch/ch/help/symbolic/mupad_ref/linalg-smithform.html
% A = [9 -36 30; -36 192 -180;30 -180 180]; 
% D = |3 0  0 | 
%     |0 12 0 | 
%     |0 0  60|

% Example 4
% Source: http://www.maplesoft.com/support/help/Maple/view.aspx?path=linalg%28deprecated%29/ismith
% A = [13 5 7;17 31 39];
% D = |1 0 0 |  
%     |0 2 0 | 

[T,D,S] = smithFormInt(A)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Following some sample polynomial matrices whose smith form you can find online
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
clc
syms x;

%Instructions: Uncomment A and B (when present)

% Example 1
% Source: Random Lecture
% A = [1+x^2 x;x 1+x ]
% D = | 1      0     |
%     | 0 x^3 + x + 1|


% Example 2
% Source: Random Lecture modified
% A = [1+x^2 x 5;x 1+x x^3; 0 7 x] 
% D = | 1  0                                     0|
%     | 0  1                                     0|
%     | 0  0  x^5 - x^4/7 + x^3 - x^2/7 - (36*x)/7|

% Example 3
% Source: % http://math.stackexchange.com/questions/77063/how-do-i-get-this-matrix-in-smith-normal-form-and-is-smith-normal-form-unique
% B = [5 2 -8 -8; -6 -3 8 8; -3 -1 3 4; 3 1 -4 -5]
% A = x*eye(size(B)) - B
% D = |1  0    0       0     |
%     |0 x+1   0       0     |
%     |0  0   x+1      0     |
%     |0  0    0   (x+1)(x-3)|

% Example 4
% Source: Random Lecture 
% B = [2 0 0 0; -1 1 0 0; 0 -1 0 -1; 1 1 1 2]
% A = x*eye(size(B)) - B
% D = |1  0     0       0       |
%     |0  1     0       0       |
%     |0  0    x-1      0       |
%     |0  0     0   (x-2)(x-1)^2|

% Example 5
% Source: Mathematica SmithDemoV6 (link above)
A = [2, -2 + x - x^3, -1 + 2*x - x^2, -x + x^3; ...
    -2 + x - x^3, -x + (1/2)*x^2 + 2*x^3 - x^4 + (1/2)*x^6, -1 - (3/2)*x + 2*x^2 - x^4 + (1/2)*x^5, 2 - x^2/2 - x^3 + x^4 - (1/2)*x^6; ...
    -1 + 2*x - x^2, -1 - (3/2)*x + 2*x^2 - x^4 + (1/2)*x^5, 15/2 - (23/2)*x + 7*x^2 - (5/2)*x^3 + (1/2)*x^4, -1 + 2*x - (3/2)*x^2 + x^4 - (1/2)*x^5]
% D = |1  0       0      0|
%     |0 x-2      0      0|
%     |0  0   (x-2)(x-3) 0|

% Example 6
% Tricky Polynomial Matrix 
% A = [31 - 160*x, -8*x - 180, 8*x + 80; ... 
%     37 - 170*x, -20*x - 190, 7*x + 90]
% D = |1 0 0|
%     |0 1 0|

[T,D,S] = smithFormPoly(A);

T = expand(T)
D = expand(D)
S = expand(S)