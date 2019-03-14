classdef examples < handle

  properties
    solver = ''
    algorithm = ''
  end
  

  methods 
    
    function this = examples(solver,algorithm)
      
      if nargin < 1 || isempty(solver)
        this.solver = 'qp';
      else
        this.solver = solver;
      end
      
    end
    
    function [x,fval,exitflag,output,ex] = getExample(this,No)
      
      switch No
        case 1
          % exmple 1
          %	Minimize	f(x1,x2) = -2x1 -x2
          %
          %	S.t.    g1(x1,x2) = x1^2 + x2^2 <= 25       (nonlinear inequality constraint)
          %         g2(x1,x2) = x1^2 - x2^2 <= 7      (nonlinear inequality constraint)
          %         0 <= x1 <= 10;  0 <= x2 <= 10     (box constraints)
          %
          %   Starting point: (x1,x2) = (1,1)         (a feasible point)
          x0 = [1;1];
          lb = [0;0];
          ub = [10;10];
          f =@(x) examples.f1(x);
          g =@(x) examples.g1(x);
          ex = fminsqp(f,x0,[],[],[],[],lb,ub,g,'display','iter');
          [x,fval,exitflag,output] = ex.solve;
          ex.postprocess;
          
        case 2
          % exmple 2  
          %	Minimize	f(x1,x2) = x1^4 - 2*x1*x1*x2 + x1*x1 + x1*x2*x2 - 2*x1 + 4
          %
          %	S.t. h(x1,x2) = x1^2 + x2^2 = 2            (nonlinear equality constraint)
          %        g(x1,x2) = 1/4*x1^2 +3/4*x2^2 <= 1  (nonlinear inequality constraint)
          %        0 <= x1 <= 4;  0 <= x2 <= 4         (box constraints)
          %
          %   Starting point: (x1,x2) = (sqrt(2),0)    (a feasible point)
          x0 = [sqrt(2);0];
          lb = [0;0];
          ub = [4;4];
          f =@(x) examples.f2(x);
          g =@(x) examples.g2(x);
          ex = fminsqp(f,x0,[],[],[],[],lb,ub,g,'display','iter');
          [x,fval,exitflag,output] = ex.solve;
          ex.postprocess;
          
        case 3
          % Example 12.3 from Arora (2012) Example 7.3, pp. 285-287 (Arora (2004) pp. 418-420):
          %
          %	Minimize	f(x1,x2) = (x1-10)^3 + (x2-20)^3
          %
          %	S.t.  g1(x1,x2) = -(x1-5)^2 -(x2-5)^2 <= -100     (nonlinear inequality constraint)
          %         g2(x1,x2) = (x1-6)^2 + (x2-5)^2 <= 82.81  (nonlinear inequality constraint)
          %         13 <= x1 <= 100;  0 <= x2 <= 100          (box constraints)
          %
          %   Starting point: (x1,x2) = (20.1,5.84)                (a infeasible point)
          %   Optimum x=[14.095;0.84296], f(x) = -6961.8
          x0 = [20.1;5.84];
          lb = [13;0];
          ub = [100;100];
          f =@(x) examples.f3(x);
          g =@(x) examples.g3(x);
          ex = fminsqp(f,x0,[],[],[],[],lb,ub,g,'display','iter');
          [x,fval,exitflag,output] = ex.solve;
          ex.postprocess;
          
        case 4
          %	Minimize	e^(x1*x2*x3*x4*x5)-1/2*(x1^3+x2^3+1)^2
          %
          %	S.t.    g1(x) = x1^2+x2^2+x3^2+x4^2+x5^2-10 = 0  (nonlinear equality constraint)
          %         g2(x) = x2*x3-5*x4*x5 = 0                (nonlinear equality constraint)
          %         g3(x) = x1^3+x2^3+1 = 0                  (nonlinear equality constraint)
          %         
          %
          %   Starting point: x = [-1.71;1.59;1.81;-0.763;-0.763]       (a feasible point) 
          %   Optimum x = [-1.8;1.7;1.9;-0.8;-0.8]
          x0 = [-1.71;1.59;1.81;-0.763;-0.763];
          lb = [-2;-2;-2;-2;-2];
          ub = [2;2;2;2;2];
          f =@(x) examples.f4(x);
          g =@(x) examples.g4(x);
          ex = fminsqp(f,x0,[],[],[],[],lb,ub,g,'display','iter');
          [x,fval,exitflag,output] = ex.solve;
          ex.postprocess;
        case 5
          %	Minimize	f(x1,x2) = x1^2+2*x2^2
          %	S.t.  g1(x1,x2) = x1 + x2 >= 1 
          x0 = [1;1];
          lb = [-10;-10];
          ub = [10;10];
          f =@(x) examples.f5(x);
          g =@(x) examples.g5(x);
          ex = fminsqp(f,x0,[],[],[],[],lb,ub,g,'display','iter');
          [x,fval,exitflag,output] = ex.solve;
          ex.postprocess;
          
        case 6
          %	Minimize	rosenberg
          x0 = [-1;3];
          lb = [-10;-10];
          ub = [10;10];
          f =@(x) examples.f6(x);
          g =[];
          ex = fminsqp(f,x0,[],[],[],[],lb,ub,g,'display','iter');
          [x,fval,exitflag,output] = ex.solve;
          ex.postprocess;
      end % switch
    end % function
  end % methods
  
  methods(Static = true, Hidden = false) 
      function [f,df] = f1(x)
        f = [];
        df = [];
        switch nargout
          case 1
            f = -2*x(1) - x(2);
          case 2
            df(1,1) = -2;
            df(2,1) = 1;
        end
      end
      
      function [g,geq,dg,dgeq] = g1(x)
        g = [];
        geq = [];
        dg = [];
        dgeq = [];
        
        switch nargout
          case {1,2}
            g(1,1)= x(1)^2 + x(2)^2 - 25;
            g(2,1) = x(1)^2 - x(2)^2 - 7;
          case {3,4}
            dg(1,1) = 2*x(1);
            dg(2,1) = 2*x(2);
            dg(1,2) = 2*x(1);
            dg(2,2) = -2*x(2);
        end
        
      end
      
      function [f,df] = f2(x)
        f = [];
        df = [];
        
        switch nargout
          case 1
            f = x(1)^4 - 2*x(1)^2*x(2) + x(1)^2 + x(1)*x(2)^2 - 2*x(1) + 4;
          case 2
            df(1,1) = 4*x(1)^3 -4*x(1)*x(2) + 2*x(1) + x(2)^2 -2;
            df(2,1) = -2*x(1)^2 + 2*x(1);
        end
      end
      
      function [g,geq,dg,dgeq] = g2(x)
        g = [];
        geq = [];
        dg = [];
        dgeq = [];
        
        switch nargout
          case {1,2}
            g(1)= 1/4*x(1)^2 + 3/4*x(2)^2-1;
            geq(1) = x(1)^2 + x(2)^2 - 2;
          case {3,4}
            dg(1,1) = 2/4*x(1);
            dg(2,1) = 6/4*x(2);
            dgeq(1,1) = 2*x(1);
            dgeq(2,1) = 2*x(2);
        end
        
      end
      
      function [f,df] = f3(x)
        f = [];
        df = [];
        
        switch nargout
          case 1
            f = (x(1)-10)^3 + (x(2)-20)^3;
          case 2
            df(1,1) =3*(x(1)-10)^2;
            df(2,1) = 3*(x(2)-20)^2;
        end
      end
      
      function [g,geq,dg,dgeq] = g3(x)
        g = [];
        geq = [];
        dg = [];
        dgeq = [];
        
        switch nargout
          case {1,2}
            g(1,1)= -(x(1)-5)^2/100 - (x(2)-5)^2/100 + 1;
            g(2,1) = (x(1)-6)^2/82.81 + (x(2)-5)^2/82.81 - 1;
          case {3,4}
            dg(1,1) = -2*(x(1)-5)/100;
            dg(2,1) = -2*(x(2)-5)/100;
            dg(1,2) = 2*(x(1)-6)/82.81;
            dg(2,2) = -2*(x(2)-5)/82.81;
            
        end
        
      end
      
      function [f,df] = f4(x)
        f = [];
        df = [];
        
        switch nargout
          case 1
            f = e^(x(1)*x(2)*x(3)*x(4)*x(5))-1/2*(x(1)^3+x(2)^3+1)^2;
          case 2
            df(1,1) = x(2)*x(3)*x(4)*x(5)*(e^(x(1)*x(2)*x(3)*x(4)*x(5)))-3*x(1)^2*(x(1)^3+x(2)^3+1);
            df(2,1) = x(1)*x(3)*x(4)*x(5)*(e^(x(1)*x(2)*x(3)*x(4)*x(5)))-3*x(2)^2*(x(1)^3+x(2)^3+1);
            df(3,1) = x(1)*x(2)*x(4)*x(5)*(e^(x(1)*x(2)*x(3)*x(4)*x(5)));
            df(4,1) = x(1)*x(2)*x(3)*x(5)*(e^(x(1)*x(2)*x(3)*x(4)*x(5)));
            df(5,1) = x(1)*x(2)*x(3)*x(4)*(e^(x(1)*x(2)*x(3)*x(4)*x(5)));
        end
      end
      
      function [g,geq,dg,dgeq] = g4(x)
        g = [];
        geq = [];
        dg = [];
        dgeq = [];
        
        switch nargout
          case {1,2}
            geq(1,1) =  x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2-10;
            geq(2,1) =  x(2)*x(3)-5*x(4)*x(5);
            geq(3,1) =  x(1)^3 + x(2)^3 + 1;
          case {3,4}
            dgeq(1,1) = 2*x(1);
            dgeq(2,1) = 2*x(2);
            dgeq(3,1) = 2*x(3);
            dgeq(4,1) = 2*x(4);
            dgeq(5,1) = 2*x(5);
            
            dgeq(1,2) = 0;
            dgeq(2,2) = x(3);
            dgeq(3,2) = x(2);
            dgeq(4,2) = -5*x(5);
            dgeq(5,2) = -5*x(4);
            
            dgeq(1,3) = 3*(x(1)^2);
            dgeq(2,3) = 3*(x(2)^2);
            dgeq(3,3) = 0;
            dgeq(4,3) = 0;
            dgeq(5,3) = 0;
        end
      end
      
      function [f,df] = f5(x)
        f = [];
        df = [];
        
        switch nargout
          case 1
            f = x(1)^2+2*x(2)^2;
          case 2
            df(1,1) = 2*x(1);
            df(2,1) = 4*x(2);
        end
      end
      
      function [g,geq,dg,dgeq] = g5(x)
        g = [];
        geq = [];
        dg = [];
        dgeq = [];
        
        switch nargout
          case {1,2}
            g(1,1)= -x(1)-x(2)+1;
          case {3,4}
            dg(1,1) = -1;
            dg(2,1) = -1;
            
        end
      end
      
      function [f,df] = f6(x)
        f = [];
        df = [];
        
        switch nargout
          case 1
            f = (1-x(1))^2+100*(x(2)-x(1)^2)^2;
          case 2
            df(1,1) = 400*x(1)^3 - 400*x(1)*x(2)+2*x(1)-2;
            df(2,1) = 200*(x(2)-x(1)^2);
        end
      end
  end % Methods
  
end % class