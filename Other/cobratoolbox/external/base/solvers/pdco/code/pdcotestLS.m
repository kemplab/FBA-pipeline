function pdcotestLS( m,n,nc )

% m=50;  n=100;  nc = 10;  pdcotestLS( m,n,nc );
% Generates a random m by n Weighted Least Squares problem
% with the last nc rows as linear constraints, and runs it on pdco.m.
% (We need nc <= m.)
%
% The problem is treated as
%
%    minimize    c'x + 1/2 ||D1*x||^2 + 1/2 ||r||^2
%      x,r
%    subject to  A*x + D2*r = b,   bl <= x <= bu,   r unconstrained,
%
% where A is an m by n matrix, c is an n-vector,
% D1 and D2 are positive-definite diagonal matrices,
% and everything relevant is the same as in pdco.m.
%
% Here, D1 and D2 are defined by d1 and d2 in the private function toydata.
% To make the i-th row of (A,b) act like a constraint, d2(i) is "small"
% (1e-3 or 1e-4).  Ordinary LS rows have d2(i) = 1.
% Other rows could have intermediate values if desired.

%-----------------------------------------------------------------------
% 23 Sep 2003: Example LS test program for pdco.m,
%              derived from pdcotestLP.m.
%              Michael Saunders, SOL, Stanford University.
%-----------------------------------------------------------------------

  [A,b,bl,bu,c,d1,d2] = toydata( m,n,nc ); % Private function below
% D  = sum(A,1);   D(find(D==0)) = 1;
% D  = sparse( 1:n, 1:n, 1./D, n, n );
% A  = A*D;                                % Normalize cols of A

  options = pdcoSet;

  x0    = 0.5*ones(n,1);     % Initial x
  y0    = zeros(m,1);        % Initial y
  z0    = ones(n,1);         % Initial z
  xsize = 1;                 % Estimate of norm(x,inf) at solution
  zsize = 1;                 % Estimate of norm(z,inf) at solution

  options.mu0       = 1e-0;  % An absolute value
  options.LSQRatol1 = 1e-8;  % For LPs, LSQR must solve quite accurately
  options.wait      = 1;     % Allow options to be reviewed before solve

  [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco( c,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );

% keyboard                   % Allow review of x,y,z, etc.
%-----------------------------------------------------------------------
% End function pdcotestLS
%-----------------------------------------------------------------------


function [A,b,bl,bu,c,d1,d2] = toydata( m,n,nc )

%        [A,b,bl,bu,c,d1,d2] = toydata( m,n,nc );
%        defines an m by n matrix A, rhs vector b, and cost vector c
%        for use with pdco.m.  The last nc rows of (A,b) are
%        treated differently to give a small residual.

%-----------------------------------------------------------------------
% 23 Sep 2003: Adapted from version in pdcotestLP.m.
%-----------------------------------------------------------------------

  rand('state',0);
  density = 0.25;
  rc      = 1e-1;

  em      = ones(m,1);
  en      = ones(n,1);
  zn      = zeros(n,1);
  bigbnd  = 1e+30;

  A       = sprand(m,n,density,rc);
  x       = en;

  gamma   = 1e-3;       % Primal regularization.
  delta   = 1e-3;       % 1e-3 or 1e-4 for LP;  1 for Least squares.

  d1      = gamma;      % Scalar for this test.
  d2      = em;         % D2 = I
  d2(m-nc+1:m) = delta; % Satisfy last nc rows more accurately.
	     
  b       = full(A*x) + d2;  % Plausible b, as if dual vector y = ones(m,1).
  c       = rand(n,1);

  bl      = zn;         % x > 0
  bu      = 10*en;      % x < 10
% bl(1)   = - bigbnd;   % Test "free variable" (no bounds)
% bu(1)   =   bigbnd;

%-----------------------------------------------------------------------
% End private function toydata
%-----------------------------------------------------------------------
