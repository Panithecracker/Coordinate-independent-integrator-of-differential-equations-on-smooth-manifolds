function M = CoordinateAdaptive(A,M0,T,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input: 
%  A: Atlas
%  M0: Initial condition vector
%  T: Final time
%  h : time step

%% Output: Matrix of estimates ordered by columns, in steps of dt units of
%time from M0 included.function M = IntegratorUniform(A,M0,T,dt)
%input: 
% A: Atlas
% M0: Initial condition
% T: Final time
% dt : time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialization:
n = ceil(T/h);
M = zeros(3,n+1);
M(:,1) = M0;

for i=1:n 

    valid_charts = GetCharts(A,M(:,i)); %valid charts (indexes) for current point
    chart_index = valid_charts(1); %pick chart 
    C = A(chart_index).Coords(M(:,i)); %coords on current chart
    % 0-USE EXTERNAL OR INTERNAL VECTOR FIELD DESCRIPTION (USER CHOICE, JUST UNCOMMENT THE ONE NEEDED):
%     V = A(chart_index).Projector(M(:,i))*f(M(:,i)); %external description
    V =  V_chart(A,chart_index,C); %intrisically, using chart basis

    %% 1-INTEGRATE ON CHART :
    Cnext = C+V*h; %Euler e.g (any works)
    
    valid_charts = valid_charts(valid_charts~=chart_index); %update valid charts

    while (~A(chart_index).inside_image(C)) %outside of coord space, change charts
        if (~isempty(valid_charts)) %alternatives available
            new_index = valid_charts(1); %pick new chart
            %% 2-CHANGE TANGENT SPACE BASIS USING TRANSITION MAP AND INTEGRATE ON THE ADAPTED NEW CHART
            V = TranslateVelocity(A,chart_index,new_index,C,V);
            chart_index = new_index;
            C = A(chart_index).Coords(M(:,i));
            Cnext = C+V*h; %Euler (any other works)
            valid_charts = valid_charts(valid_charts~=chart_index); %update "stack" of valid charts
        else
            fprintf("Tested with all charts, decrease step size!\n") %charts consumed, notify and exit 
            return
        end
    end
    %% 3-COORDINATE ESTIMATED, NOW INVERT TO OBTAIN THE POINT ESTIMATE ON THE MANIFOLD
    M(:,i+1) = A(chart_index).Manifold(Cnext);
end
end



function vj = TranslateVelocity(A,i,j,c_i,vi) %obtains coordinate in chart j for point with coord i in chart i
    %expresses a tangent vector from basis of chart i with respect to basis
    %of chart j. In other words, it obtains the velocity of a particle in
    %the chart j from the coordinate c_i, based on chart c_j.
    coords_i = sym('c',[1 2]);
    transition = A(j).Coords(A(i).Manifold(coords_i'));
    J = jacobian(transition,coords_i); 
    vj = double(subs(J,coords_i,c_i'))*vi; %the change of basis matrix is the jacobian of transition at c_i
end

function list = GetCharts(A,p)  %returns indexes of charts from atlas A , which cover the input point p 
    list = []; 
    for i = 1:size(A,1)
        if A(i,1).inside_domain(p) == 1
            list = [list i];                                                               
        end
    end
end 

function v = f(x) %user can define the vector field in Rn here (external description) 
    h2 = 1/sqrt(2);
    h3 = 0;
    h4 = 1/sqrt(2);
    v = [h4*x(2)-h3*x(3);h2*x(3)-h4*x(1);h3*x(1)-h2*x(2)];
end

function v = V_chart(A,chart_index,coord) %user can define the vector field through the charts of atlas A here (intrinsic description) 
vect = f(A(chart_index).Manifold(coord));
if chart_index == 1 || chart_index == 2
v = [vect(1);vect(2)];
end
if chart_index == 3 || chart_index == 4
v = [vect(1);vect(3)];
end
if chart_index == 5 || chart_index == 6
v = [vect(2);vect(3)];
end
end 



