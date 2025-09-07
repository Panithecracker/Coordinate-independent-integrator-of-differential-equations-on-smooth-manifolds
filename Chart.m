classdef Chart
    properties
        %requested:
        domain      %domain of coord function
        image       % image of coord function
        phi         % manifold -> coords
        phi_inv     % coords -> manifold 
        n           % dimension of embedding space
        k           % dimension of manifold
        name        % name of the chart (optional)
        %internally inferred:
        coords_sym  % Symbolic variables for the coordinates
        Jacobian            %Jacobian of inverse map
    end
    
    methods
        function obj = Chart(d,i, p, pinv, n, k, name)
            % Parameterized Constructor 
            obj.domain = d;
            obj.image = i;
            obj.phi = p;
            obj.phi_inv = pinv;
            obj.n = n;
            obj.k = k;
            
            % Define symbolic variables for the k-dimensional coordinate space
            % 'c' is used as a prefix for coordinate variables, e.g., c1, c2, ...
            obj.coords_sym = sym('c', [1 k]); 
            
            % Assign the symbolic inverse map
            % The user provides an anonymous function that accepts these symbolic variables
            % and returns the symbolic expression for phi_inv.
            phi_inv_sym = obj.phi_inv(obj.coords_sym');
            obj.Jacobian = jacobian(phi_inv_sym, obj.coords_sym);

            if nargin < 6
                obj.name = 'Unnamed';
            else
                obj.name = name;
            end
        end
        
        function in = inside_domain(obj,p)
            % Check if the point in Rn lies in the domain
            in = obj.domain(p);
        end
        
        function in = inside_image(obj,c)
            in = obj.image(c);
            % Check if the point in Rk lies in the image of phi
        end
        
        function c = Coords(obj,p) %finds coordinates on chart
            c = obj.phi(p);
        end
        
        function p = Manifold(obj,c) %evaluates phi_inv
            p = obj.phi_inv(c);
        end
        
        function Jval = JacobianAt(obj,c)
            % Evaluates Jacobian at coordinate c            
            Jval = double(subs(obj.Jacobian,obj.coords_sym,c'));
        end
        function A = Projector(obj,p) %useful to convert a vector field of Rn into a vector field on M
            %returns the "projector" onto tangent space at p of M with
            %respect to this chart basis
            %(k by n matrix )
            c = obj.phi(p);
            J = obj.JacobianAt(c);
            A = inv(J'*J)*J';
        end
    end
end


function vj = TranslateVelocity(A,i,j,c_i,vi) 
    %expresses a tangent vector given in chart i into chart j basis.
    coords_i = sym('c',[1 2]);
    transition = A(j).Coords(A(i).Manifold(coords_i'));
    J = jacobian(transition,coords_i); 
    vj = double(subs(J,coords_i,c_i'))*vi; 
end


 % That is, given a particle moving at coordinate velocity v_i away from c_i on chart i 
    % it finds the respective coordinate
    %velocity, vj from coordinate c_j in chart j.
    
    %the change of basis matrix is the jacobian of transition at c_i
function list = GetCharts(Atlas,p) %for a given point and atlas, it returns the list of indexes from charts on the atlas , that contain the point
    list = []; 
    for i = 1:size(Atlas,1)
        if Atlas(i,1).inside_domain(p) == 1
            list = [list i];
        end
    end
end
