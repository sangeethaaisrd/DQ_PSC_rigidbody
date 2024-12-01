classdef quaternion
  
  properties
    s (1,1) double;  % Real (scalar) part
    v (3,1) double;  % Imaginary vector (3x1) part
  end
  
  methods
  
    function q = quaternion(s, v)
      % Constructor
      if nargin == 0
        q.s = 0;
        q.v = zeros(3,1);
      elseif nargin == 1
          if isscalar(s)
              q.s = s;
              q.v = zeros(3,1);
          else
              q.s = 0;
              q.v = s;
          end
      else
        q.s = s;
        q.v = v(:);
      end
    end
    
    function q_out = uplus(q)
      % Quaternion unary addition
      q_out = q;
    end

    function q_out = plus(q1, q2)
      % Quaternion addition
      q_out = quaternion();
      if isa(q1,'double')
          q1 = quaternion(q1,[0 0 0]);
      elseif isa(q2,'double')
          q2 = quaternion(q2,[0 0 0]);
      end
      q_out.s = q1.s + q2.s;
      q_out.v = q1.v + q2.v;
    end
    
    function q_out = uminus(q)
      % Quaternion negation
      q_out = quaternion();
      q_out.s = -q.s;
      q_out.v = -q.v;
    end    
    
    function q_out = minus(q1,q2)
      % Quaternion subtraction
      q_out = quaternion();
      if isa(q1,'double')
          q1 = quaternion(q1,[0 0 0]);
      elseif isa(q2,'double')
          q2 = quaternion(q2,[0 0 0]);
      end
      q_out.s = q1.s - q2.s;
      q_out.v = q1.v - q2.v;
    end
        
    function q_out = mtimes(q1, q2)
      % Quaternion multiplication
      if isa(q1,'double')
          q1 = quaternion(q1,[0 0 0]);
      elseif isa(q2,'double')
          q2 = quaternion(q2,[0 0 0]);
      end
      s1 = q1.s;
      v1 = q1.v;
      s2 = q2.s;
      v2 = q2.v;
      q_out = quaternion();
      q_out.s = s1 * s2 - dot(v1, v2);  % Scalar part of the product
      q_out.v = (s1 * v2 + s2 * v1 + cross(v1, v2));  % Vector part of the product
    end

    function q_out = conj(q)
        % Quaternion conjugate
        q_out = quaternion();
        q_out.s =  q.s;
        q_out.v = -q.v;
    end

    function q_out = cross(q1,q2)
        % Quaternion cross product: (q1q2-q2q1)/2
        q_out = quaternion();
        q_out.s = 0;
        q_out.v = cross(q1.v,q2.v);
    end

    function dotprdt = dot(q1,q2)
        % dot product of quaternions: (q1 q2* + q1* q2)/2
        dotprdt = dot([q1.s;q1.v],[q2.s;q2.v]);
    end

    function abs_q = norm(q)
        abs_q = norm([q.s;q.v]);
    end

    function q_out = inv(q)
      % Quaternion inverse
      q_out = conj(q);  
      qnorm2 = norm(q)^2;
      if qnorm2 ~= 0
        q_out.s = q_out.s/qnorm2;
        q_out.v = q_out.v/qnorm2;
      else
          error('Zero quaternion! Inverse does not exist.')
      end
    end
    
    function q_out = mrdivide(q1, q2)
      % Quaternion division
      if isa(q1,'double')
          q1 = quaternion(q1,[0 0 0]);
      elseif isa(q2,'double')
          q2 = quaternion(q2,[0 0 0]);
      end
      q_out = q1*inv(q2);  
    end
    
    function q_out = normalize(q)
        q_out = q/norm(q);
    end
    
  end
  
end

