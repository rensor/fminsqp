classdef fminsqp < handle
  
  properties(Constant)
    name = 'fminsqp';
  end

  properties
    options = [];
  end
  
  properties(SetAccess = public, Hidden = true)
    
    % Inputs
    fun = [];
    x0 = [];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = [];
    nDV = [];
    
    % Second order Quasi-Newton Methods
    Hinv = []; % Inverse of hessian DFP method
    H = []; % Hessian BFGS method
    iH = 0; % Hessian update counter
    
    % Merit function variables
    f0 = []; % Initial, unpeanalized objective function value
    nGnl = []; % number of non-linear inequality constraints
    aFac = []; % Scaling parameter for merit function
    lambda = []; % Lagrange variables for "Augmented Lagrange (AL) method"
    
    % Global convergence filter 
    filter = struct();
    
    % Iteration history
    history = struct();
    
    % Switch
    initialized = false;
    
    
  end
  
  methods
    
    % Construct
    function this = fminsqp(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,varargin)
      
      % initialize data structures
      if nargin < 1 || ~isa(fun,'function_handle')
        error([this.name,': fun is required to be a function handle'])
      else
        this.fun = fun;
      end
      
      if nargin < 2 || isempty(x0)
        this.x0 = [];
      else
        if isnumeric(x0)
          this.x0 = x0(:);
        else
          error([this.name,' x0 is required to be of type numeric'])
        end
      end
      
      if nargin < 3 || isempty(A)
        this.A = [];
      else
        if isnumeric(A)
          this.A = A;
        else
          error([this.name,' A is required to be of type numeric'])
        end
      end
      
      if nargin < 4 || isempty(b)
        this.b = [];
      else
        if isnumeric(b)
          this.b = b(:);
        else
          error([this.name,' b is required to be of type numeric'])
        end
      end
      
      if size(this.A,1) ~= size(this.b,1)
        error([this.name,' A and b must contain equal number of rows'])
      end
      
      if nargin < 5 || isempty(Aeq)
        this.Aeq = [];
      else
        if isnumeric(Aeq)
          this.Aeq = Aeq;
        else
          error([this.name,' Aeq is required to be of type numeric'])
        end
      end
      
      if nargin < 6 || isempty(beq)
        this.beq = [];
      else
        if isnumeric(beq)
          this.beq = beq(:);
        else
          error([this.name,' beq is required to be of type numeric'])
        end
      end
      
      if size(this.Aeq,1) ~= size(this.beq,1)
        error([this.name,' Aeq and beq must contain equal number of rows'])
      end
          
      if nargin < 7 || isempty(lb)
        this.lb = [];
      else
        if isnumeric(lb)
          this.lb = lb(:);
        else
          error([this.name,' lb (lower bound) is required to be of type numeric'])
        end
      end
      
      if nargin < 8 || isempty(ub)
        this.ub = [];
      else
        if isnumeric(ub)
          this.ub = ub(:);
        else
          error([this.name,' ub (upper bound) is required to be of type numeric'])
        end
      end
      
      if nargin < 9 || isempty(nonlcon)
        this.nonlcon = [];
      elseif ~isempty(nonlcon)
        if ~isa(nonlcon,'function_handle')
          error([this.name,' nonlcon is required to be a function handle'])
        else
          this.nonlcon = nonlcon;
        end
      end
      
      % check that all sizes match
      if ~isempty(this.x0)
        this.nDV = numel(this.x0);
      end
      
      if ~isempty(this.lb) && ~isempty(this.ub)
        if numel(this.lb) ~= numel(this.ub)
          error([this.name,' lb and ub not equal dimensions'])
        end
      end
      
      if ~isempty(this.lb) && ~isempty(this.nDV)
        if numel(this.lb) ~= this.nDV
          error([this.name,' x0 and lb not equal dimensions'])
        end
      end
      
      if ~isempty(this.ub) && ~isempty(this.nDV)
        if numel(this.ub) ~= this.nDV
          error([this.name,' x0 and ub not equal dimensions'])
        end 
      end
      
      if ~isempty(this.A) && ~isempty(this.nDV)
        if size(this.A,2) ~= this.nDV
          error([this.name,' Columns of A(',num2str(size(this.A,2)),') does not match number of design variables(',num2str(this.nDV),')'])
        end
      elseif ~isempty(this.A) && isempty(this.nDV)
        this.nDV = size(this.A,2);
      end
      
      if ~isempty(this.Aeq) && ~isempty(this.nDV)
        if size(this.Aeq,2) ~= this.nDV
          error([this.name,' Columns of Aeq(',num2str(size(this.A,2)),') does not match number of design variables(',num2str(this.nDV),')'])
        end
      elseif ~isempty(this.Aeq) && isempty(this.nDV)
        this.nDV = size(this.Aeq,2);
      end
      
      % initialize options structure
      this.options = fminsqp.sqpoptions(varargin);
      
      % Initialize global convergence filter
      this.filter = fminsqp.initializeGlobalConvergenceFilter(this.options);
      
      % We made it this far
      this.initialized = true;
    end
    
    % Main function
    function [x,fval,exitflag,output] = solve(this)
      % Assume the code failed
      exitflag = -1;
      
      if strcmpi(this.options.Display,'iter')
        fprintf('*********************************************************************************************************')
        fprintf('\n \t \t \t \t \t \t fminsqp optimizer with global convergence filter')
        fprintf('\n*********************************************************************************************************\n')
        fprintf('\t %10s \t\t %10s \t \t %10s \t \t   %10s \t \t   %10s \n','f(x)','Max inf', 'Norm dx', 'nFeval','IterNo');
      end
        
      % Allocate iteration history array
      % Store function values, maximum infeasibility from non-linear
      % constraints, norm of design variable change
      this.history.f = zeros(this.options.MaxIterations,1);
      this.history.xnorm = zeros(this.options.MaxIterations,1);
      
      x = this.x0(:);
      
      % Initialize the current box constraints with upper and lower bounds.
      xLcur = this.lb(:);
      xUcur = this.ub(:);
      
      % Initialize "old" design variables
      xold1 = x;
      xold2 = x;
      
      
      % evaluate initial non-linear in-equality constraints
      if ~isempty(this.nonlcon)
        [gnl, gnleq] = this.nonlcon(x);
        % Determine number of constraints
        this.nGnl = numel(gnl) + numel(gnleq)*2;
        % Allocate iteration history for maximum infeasibility
        this.history.maxInf = zeros(this.options.MaxIterations,1);
      else
        this.nGnl = 0;
        this.history.maxInf = [];
      end
      
      % evaluate objective function at initial point
      this.f0  = this.fun(x);
      
      % Get scaling factor for merit function
      this.aFac = max([abs(this.f0),1]);
      
      % Evaluate merit function
      f = this.getObj(x,false);
      % Store initial value for filter
      this.filter.initF = f;
      % Store "old" value for convergence check
      fOld = f;
      df = zeros(numel(x),1);
      deltax = zeros(numel(x),1);
      this.H = eye(this.nDV);
      
      % Set counters and switches
      nFeval = 1;
      iterNo = 0;
      optimize = true;
      % Main loop
      while optimize
        
        % update iteration counter
        iterNo = iterNo + 1;
        
        % Store previous gradients (minus 1)
        dfm1 = df;
        % evaluate gradients
        [~,df] = this.getObj(x,true);
        [g] = this.getConstraints(x,false);
        [~,dg] = this.getConstraints(x,true);
        
        if iterNo > 1
          this.updateHessian(df,dfm1,deltax)
        end
      
        if ~isempty(this.A)
          % Assemble linear and non-linear in-equality constraints
          Am = zeros(size(this.A,1)+this.nGnl,size(this.A,2)); % Allocate
          Am(1:size(this.A,1),1:size(this.A,2)) = this.A; % Load linear
          bm = zeros(size(this.b,1)+this.nGnl,1); % Allocate
          bm(1:size(this.b,1)) = this.b; % Load linear
          
          if ~isempty(this.nonlcon)
            Am(size(this.A,1)+1:end,:) = dg; % Load non-linear
            bm(size(this.b,1)+1:end) = dg*x-g; % Load non-linear
          end
        else % Only non-linear
          if ~isempty(this.nonlcon)
            Am = dg;
            bm = dg*x-g;
          else
            Am = [];
            bm = [];
          end
        end
        
        % Expand linear equality constraints to account for slag variables
        if ~isempty(this.Aeq)
          Ameq = zeros(size(this.Aeq,1),size(this.Aeq,2)+this.nGnl);
          Ameq(1:size(this.Aeq,1),1:size(this.Aeq,2)) = this.Aeq;
        else
          Ameq = [];
        end
        
        % update move-limits
        reduceSwitch = false;
        [xLcur, xUcur] = this.AdaptiveMoveLimit(x, xLcur, xUcur, this.lb, this.ub, this.options.MoveLimit ,this.options.MoveLimitReduce, this.options.MoveLimitExpand, xold1, xold2, reduceSwitch);
        
        % update old values
        xold2 = xold1;
        xold1 = x;
        
        % Set switches
        backtrack = true;
        
        % Inner loop
        while backtrack
          
          % Set lower and upper bounds for the lp problem
          xL = xLcur(:);
          xU = xUcur(:);
        
          % Call optimizer
          [xNew,exitflag] = this.qpSolver(df,Am,bm,Ameq,this.beq,xL,xU,x);
          
          % Determine design deltas
          deltax = xNew - x;
          deltanorm = norm(deltax);
          
          % Determine optimality Norm for convergence check
          optimalityNorm = norm(df);
          
          % evaluate constraints
          [g] = this.getConstraints(xNew,false);
          
          % Evaluate objective function at new point
          nFeval = nFeval + 1;
          
          [f] = this.getObj(xNew,false);
          
          % Determine delta for merit function
          deltaf = fOld - f;
          
          % Determine maximum infeasibility for current design
          this.filter.h = max([g;0]); 
          % Store objective function value for current design
          this.filter.f = f/this.filter.initF;
          
          % Evaluate current (h,f) point against the convergence filter
          [this.filter] = this.EvaluateCurrentDesignPointToFilter(this.filter);
          
          % Assume that we don't want to update the convergence filter
          AddToFilter = false;
          if (this.filter.PointAcceptedByFilter)
            % Determine Delta values for filter checks
            deltaQ = -df'*deltax-0.5*deltax'*this.H*deltax;
            
            % Check if we should add accept the current point
            if ( (deltaf<this.filter.sigma*deltaQ) && (deltaQ>0.0) )
              % Not accepted
              reduceSwitch = true;
            else
              % Accepted
              reduceSwitch = false;
              backtrack = false;
              % Check if we should but the new point (h,f) into the filter
              if(deltaQ <= 0.0)
                AddToFilter = true;
              end
            end
          else
            reduceSwitch = true;
          end
          
          if reduceSwitch
            [xLcur, xUcur] = this.AdaptiveMoveLimit(x, xLcur, xUcur,this.lb, this.ub, this.options.MoveLimit ,this.options.MoveLimitReduce, this.options.MoveLimitExpand, xold1, xold2, reduceSwitch);
          end
          
          % check for convergence
          if (optimalityNorm <= this.options.OptimalityTolerance) || (deltanorm <=this.options.StepTolerance) || (iterNo >= this.options.MaxIterations) || (nFeval >= this.options.MaxFunctionEvaluations)
            optimize = false;
            backtrack = false;
            exitflag = 1;
            fval = this.fun(x);
          end
          
        end % Inner loop
        
        % Does the new point(h,f) qualify to be added to the filter?
        if (AddToFilter)
          [ this.filter ] = this.UpdateFilter(this.filter, this.filter.h, this.filter.f);
        end
        
        % Update design variables
        x = xNew;

        % Update "old" design
        fOld = f;
        this.history.f(iterNo) = f;
        this.history.xnorm(iterNo) = deltanorm;
        maxInf = max([g;0]);
        this.history.maxInf(iterNo) = maxInf;
        this.history.nIter = iterNo;
        this.history.nFeval = nFeval;
        
        if strcmpi(this.options.Display,'iter')
            fprintf('\t %6.4e \t \t %6.4e \t \t %6.4e \t \t %10i \t \t %10i \n' ,f, maxInf, deltanorm, nFeval ,iterNo);
        end
        
      end % Main loop
      
      this.history.f(iterNo+1:end)=[];
      this.history.xnorm(iterNo+1:end)=[];
      this.history.maxInf(iterNo+1:end)=[];
      output.history = this.history;
      
    end % Solve function
    
    function postprocess(this)
      % Save current "default" window style
      defaultWindowStyle=get(0,'DefaultFigureWindowStyle');
      % Set new window style to docked
      set(0,'DefaultFigureWindowStyle','docked')
      
      % Make iteration vector
      ivec = 1:this.history.nIter;
      
      f1=figure();
      plot(ivec,this.history.f)
      title('Objective')
      xlabel('Iteration Number')
      ylabel('Objective value')
      
      figure();
      plot(ivec,this.history.xnorm)
      title('Design change norm')
      xlabel('Iteration Number')
      yl=ylabel('Norm dx');
      set(yl,'Interpreter','none')
      
      figure();
      plot(ivec,this.history.maxInf)
      title('Maximum infeasibility')
      xlabel('Iteration Number')
      ylabel('-')
      
      % Jump back to figure 1
      figure(f1)
      % Restore default window style
      set(0,'DefaultFigureWindowStyle',defaultWindowStyle)
    end
    
  end % methods
  
  
  methods (Hidden = true)
    
    function [fval,df] = getObj(this,x,doDSA)
      fval = [];
      df = [];
      
      if ~doDSA
        fval = this.fun(x);
      else
        [~,df] = this.fun(x);
      end
    end
    
    function [g,dg] = getConstraints(this,x,doDSA)
      dg = [];
      g = [];
      if ~isempty(this.nonlcon)
        if ~doDSA
            [gn, gneq] = this.nonlcon(x);
            g = [gn(:);gneq(:);-gneq(:)];
        else
            [~,~,dgnl,dgneq] = this.nonlcon(x);
            dg = zeros(this.nGnl,this.nDV);
            dg(:,1:this.nDV) = [dgnl';dgneq';-dgneq'];
        end
      end
    end
    
    function updateHessian(this,df,dfm1,deltax)
      if this.iH >= this.options.HessianRest
        this.H = eye(this.nDV);
        this.iH = 0;
      end
      
      s = deltax;
      y = df-dfm1;
      sHs = s'*this.H*s;
      sy = s'*y;
      if sy >=0.2*sHs
        theta = 1;
        r = y;
      else
        theta = 0.8*sHs/(sHs-sy);
        r = theta*y+(1-theta)*this.H*s;
      end
      this.iH = this.iH + 1;
      this.H = this.H - this.H*s*s'*this.H/(sHs)+r*r'/(s'*r);
      
    end
    
    function [x,exitflag] = qpSolver(this,df,A,b,Aeq,beq,lb,ub,x)
      % Here you can add your own solvers
      switch this.options.Solver
        
        case 'quadprog' % MathWorks solver
          quadprogOptions = optimoptions('quadprog','Display','off');
          [x,~, exitflag]= quadprog(this.H,df,A,b,Aeq,beq,lb,ub,x,quadprogOptions);
          
        case 'qp' % Octave lp solver
          if ~isempty(A)
            [x, ~, info] = qp(x,this.H,df,Aeq,beq,lb,ub,[],A,b);
          else
            [x, ~, info] = qp(x,this.H,df,Aeq,beq,lb,ub);
          end
          
          switch info.info
            case {2,3,6}
              exitflag = 0;
            otherwise
              exitflag = 1;
          end

        otherwise
          error([this.name,' Unknown QP solver specified: ',this.options.Solver])
      end
    end
    
  end
  
  methods (Static = true, Hidden = true)
    
    % options initialization
    function options = sqpoptions(input)
      % Here you can add new options if needed
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
      
      % Set parameters
      p.addParameter('Solver','qp',  @(x) checkEmpetyOrChar(x));
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxFunctionEvaluations',100,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MaxIterations',100,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('InfeasibilityPenalization',100,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('OptimalityTolerance',1e-6,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-10,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('HessianRest',100,  @(x) checkEmptyOrNumericPositive(x));
      
      % Move-limit parameters
      p.addParameter('MoveLimitMethod','adaptive',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MoveLimit',0.1,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MoveLimitExpand',1.1,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MoveLimitReduce',0.25,  @(x) checkEmptyOrNumericPositive(x));
      
      % Define global convergene filter parameters
      p.addParameter('MaxInfeasibility',1e-3,  @(x) checkEmptyOrNumericPositive(x));
      
      % pars input
      if nargin < 1 || isempty(input)
        parse(p);
      else
        parse(p,input{:});
      end
      
      % Output results to options structure
      options = p.Results;
      
    end
    
    function filter = initializeGlobalConvergenceFilter(options)
      % SLP Filter settings
      filter.SmallVal = 1.0e-6; % Used for pertubation from zero and one
      filter.gamma = filter.SmallVal;
      filter.beta = 1.0 - filter.SmallVal;
      filter.sigma = 2*filter.SmallVal;
      filter.delta = filter.SmallVal;
      filter.vals = zeros(options.MaxIterations,2); % Filter is initialized
      filter.nVals = 1; % Number of points in the filter
      filter.PointAcceptedByFilter = false; % Is the current filter point accepted 
      filter.h = 1.0e30; % Current infeasibility point, large initialization
      filter.f = 1.0e30; % Current objective point, large initialization
      filter.initF = 0.0; % Initial objective function value
      
      % Initial values in SLP filter
      if options.MaxInfeasibility < inf
        % Here MaxInfeasibility sets the maximum acceptable infeasibility
        % wrt., the global constraints (merit constrints)
        % Low values e.g., 1e-10 can potentially slow the convergence rate
        % and/or result in the optimizer getting trapped in a local minimum
        
        filter.vals(1,1) = options.MaxInfeasibility; % Maximum infeasibility
        filter.vals(1,2) = -inf;  % Related obj value
      else
        % Here, the convergence filter will steadily add new points to the
        % filter to ensure stable convergence
        
        filter.vals(1,1) = inf;  % Maximum infeasibility
        filter.vals(1,2) = inf;  % Related obj value
      end
    end
    
    % SLP Global Convergence filter function
    function [filter] = EvaluateCurrentDesignPointToFilter(filter)
    % Extract current point
    h = filter.h;
    f = filter.f;
    % Loop for all values in filter
      for ii = 1:filter.nVals
        hi = filter.vals(ii,1);
        fi = filter.vals(ii,2);
        if ( (h <= hi*filter.beta) || (f+filter.gamma*h) <= fi )
          filter.PointAcceptedByFilter = true;
        else
          filter.PointAcceptedByFilter = false;
          break
        end % Accepted or not
      end % Loop for filter values
    end
    
    % SLP Global Convergence filter function
    function [newFilter] = UpdateFilter(filter, hk, fk)
    % This function adds the pair (hk,fk) to the filter.
    % In this process, it determines wheater the new point dominates 
    % some of the points already in the filter. If so, it removes these values
    % and adds the new point to the filter.
      newFilter = filter;
      newFilter.nVals = 0;
      newFilter.vals(:,:) = 0;
      Update = true;
      ii = 0;
      if (filter.nVals >= 1)
        while (Update)
            ii = ii + 1;
            hi = filter.vals(ii,1);
            fi = filter.vals(ii,2);
          % Dominated or not?
          if ( (hk <= hi) && (fk <= (fi)) )
            
          else % Copy old data to new filter
            newFilter.nVals = newFilter.nVals + 1;
            newFilter.vals(newFilter.nVals,1) = hi;
            newFilter.vals(newFilter.nVals,2) = fi;
          end
          
          if (ii >= filter.nVals)
            Update = false;
          end
        end 
      end
      
      % Add new values to filter
      newFilter.nVals = newFilter.nVals + 1;
      newFilter.vals(newFilter.nVals,1) = hk;
      newFilter.vals(newFilter.nVals,2) = fk;
    end
    
    % Adaptive move-limit algorithm
    function [xLcur, xUcur] = AdaptiveMoveLimit(x, xLcur, xUcur, xLorg, xUorg, moveLimit ,reduceFac, expandFac, xold1, xold2, reduceSwitch)
      
      if (reduceSwitch)
        Expand = reduceFac;
        Reduction = reduceFac;
      else
        Reduction = reduceFac;
        Expand = expandFac;
      end
      
      nDv = numel(x);
      for dvNo = 1:nDv
        delta = (xUcur(dvNo)-xLcur(dvNo))/2; % This was the previous allowable change
        % Use the iteration history to determine whether we have oscillations
        % in the design variables
        if (abs(x(dvNo)-xold1(dvNo)) > 1.e-10)
          s1 = (xold1(dvNo)-xold2(dvNo)) / (x(dvNo)-xold1(dvNo));
          if (s1 < 0.0)
            delta = delta*Reduction;      % oscillation, slow increase
          else
            delta = delta*Expand;       % Stable, allow more rapid increase
          end
        else
          delta = delta*moveLimit;
        end
        dmax = (xUorg(dvNo)-xLorg(dvNo))*moveLimit;
        if (delta > dmax) 
          delta = dmax;
        end
        % Initial extimate of lower and upper bound on x(i)
        xLcur(dvNo) = x(dvNo) - delta;
        xUcur(dvNo) = x(dvNo) + delta;
        % Make sure we are within the feasible domain
        xLcur(dvNo) = max(xLcur(dvNo),xLorg(dvNo));
        xUcur(dvNo) = min(xUcur(dvNo),xUorg(dvNo));
        % Take care of extremely small design changes where the bounds may be interchanged 
        if (xLcur(dvNo) >= xUcur(dvNo)); xLcur(dvNo) = 0.9999999*xUcur(dvNo); end;
        if (xUcur(dvNo) <= xLcur(dvNo)); xUcur(dvNo) = 1.0000001*xLcur(dvNo); end;
      end
    end
    
    
  end
  
end