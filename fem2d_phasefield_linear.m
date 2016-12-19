function fem2d_phasefield_linear ( nx, ny, dt, cond, seed )

%% Input instructions:
%  nx : x descritization
%  ny : y descritization
%  dt : time descritization
%  cond : Type of BCs => 'N' for Neumann
%                        'D' for Dirichlet
%                        'PB' for Periodic
%  seed : 1 for seed in the center
%         2 for seed at the corner
%
%
%*****************************************************************************80

t=30;
epsilon = 10;
m = 0;
par_a = - 36/(epsilon^2);
par_b = 54/(epsilon^2) + 6 * m / epsilon;
par_c = -(18 / epsilon^2 + 6 * m / epsilon);

%%

phi = zeros(nx,ny);
if seed == 1
    phi(round(nx/3):floor(nx/3*2),round(ny/3):floor(ny/3*2)) = 1;
elseif seed == 2
    phi(1:round(nx/2), : ) = 1;
else
    phi(round(1):floor(nx/3),round(1):floor(ny/3)) = 1;
end
phi = reshape(phi,[],1);
nt = t/dt;

%%

element_order = 3;

xl = 0.0;
xr = .5;
yb = 0.0;
yt = .5;

%
%  NODE COORDINATES
%
%  Numbering of nodes is suggested by the following 5x10 example:
%
%    J=5 | K=41  K=42 ... K=50
%    ... |
%    J=2 | K=11  K=12 ... K=20
%    J=1 | K= 1  K= 2     K=10
%        +--------------------
%          I= 1  I= 2 ... I=10
%
node_num = nx * ny;

fprintf ( 1, '  Number of nodes =          %d\n', node_num );

node_xy = zeros(2,node_num);

k = 0;
for j = 1 : ny
    for i = 1 : nx
        
        k = k + 1;
        
        node_xy(1,k) = ( ( nx - i     ) * xl   ...
            + (      i - 1 ) * xr ) ...
            / ( nx     - 1 );
        
        node_xy(2,k) = ( ( ny - j     ) * yb   ...
            + (      j - 1 ) * yt ) ...
            / ( ny     - 1 );
        
    end
end
%
%  ELEMENT array
%
%  Organize the nodes into a grid of 3-node triangles.
%  Here is part of the diagram for a 5x10 example:
%
%    |  \ |  \ |  \ |
%    |   \|   \|   \|
%   21---22---23---24--
%    |\ 8 |\10 |\12 |
%    | \  | \  | \  |
%    |  \ |  \ |  \ |  \ |
%    |  7\|  9\| 11\|   \|
%   11---12---13---14---15---16---17---18---19---20
%    |\ 2 |\ 4 |\ 6 |\  8|                   |\ 18|
%    | \  | \  | \  | \  |                   | \  |
%    |  \ |  \ |  \ |  \ |      ...          |  \ |
%    |  1\|  3\|  5\| 7 \|                   |17 \|
%    1----2----3----4----5----6----7----8----9---10
%
element_num = 2 * ( nx - 1 ) * ( ny - 1 );

fprintf ( 1, '  Number of elements =       %d\n', element_num );

element_node = zeros ( element_order, element_num );

k = 0;

for j = 1 : ny - 1
    for i = 1 : nx - 1
        %
        %     (I,J+1)-
        %      |  \
        %      |   \
        %      |    \  |
        %    (I,J)---(I+1,J)
        %
        k = k + 1;
        element_node(1,k) = i     + ( j - 1 ) * nx;
        element_node(2,k) = i + 1 + ( j - 1 ) * nx;
        element_node(3,k) = i     +   j       * nx;
        %
        %    (I,J+1)--(I+1,J+1)
        %      |  \    |
        %          \   |
        %           \  |
        %         -(I+1,J)
        %
        k = k + 1;
        element_node(1,k) = i + 1 +   j       * nx;
        element_node(2,k) = i     +   j       * nx;
        element_node(3,k) = i + 1 + ( j - 1 ) * nx;
        
    end
end
%
%  ASSEMBLE THE SYSTEM
%
%  Assemble the coefficient matrix A and the right-hand side B of the
%  finite element equations, ignoring boundary conditions.
%
a = zeros(node_num, node_num);
c = a;

for e = 1 : element_num
    
    gi1 = element_node(1,e);
    gi2 = element_node(2,e);
    gi3 = element_node(3,e);
    
    area = 0.5 * ...
        ( node_xy(1,gi1) * ( node_xy(2,gi2) - node_xy(2,gi3) ) ...
        + node_xy(1,gi2) * ( node_xy(2,gi3) - node_xy(2,gi1) ) ...
        + node_xy(1,gi3) * ( node_xy(2,gi1) - node_xy(2,gi2) ) );
    %
    %  Consider each quadrature point.
    %  Here, we use the midside nodes as quadrature points.
    %
    for q1 = 1 : 3
        
        q2 = mod ( q1, 3 ) + 1;
        
        nq1 = element_node(q1,e);
        nq2 = element_node(q2,e);
        
        xq = 0.5 * ( node_xy(1,nq1) + node_xy(1,nq2) );
        yq = 0.5 * ( node_xy(2,nq1) + node_xy(2,nq2) );
        wq = 1.0 / 3.0;
        %
        %  Consider each test function in the element.
        %
        for ti1 = 1 : element_order
            
            ti2 = mod ( ti1,     3 ) + 1;
            ti3 = mod ( ti1 + 1, 3 ) + 1;
            
            nti1 = element_node(ti1,e);
            nti2 = element_node(ti2,e);
            nti3 = element_node(ti3,e);
            
            qi = 0.5 * ( ...
                ( node_xy(1,nti3) - node_xy(1,nti2) ) * ( yq - node_xy(2,nti2) ) ...
                - ( node_xy(2,nti3) - node_xy(2,nti2) ) * ( xq - node_xy(1,nti2) ) ) ...
                / area;
            dqidx = - 0.5 * ( node_xy(2,nti3) - node_xy(2,nti2) ) / area;
            dqidy =   0.5 * ( node_xy(1,nti3) - node_xy(1,nti2) ) / area;
            
            
            %
            %  Consider each basis function in the element.
            %
            for tj1 = 1 : element_order
                
                tj2 = mod ( tj1,     3 ) + 1;
                tj3 = mod ( tj1 + 1, 3 ) + 1;
                
                ntj1 = element_node(tj1,e);
                ntj2 = element_node(tj2,e);
                ntj3 = element_node(tj3,e);
                
                qj = 0.5 * ( ...
                    ( node_xy(1,ntj3) - node_xy(1,ntj2) ) * ( yq - node_xy(2,ntj2) ) ...
                    - ( node_xy(2,ntj3) - node_xy(2,ntj2) ) * ( xq - node_xy(1,ntj2) ) ) ...
                    / area;
                dqjdx = - 0.5 * ( node_xy(2,ntj3) - node_xy(2,ntj2) ) / area;
                dqjdy =   0.5 * ( node_xy(1,ntj3) - node_xy(1,ntj2) ) / area;
                
                c(nti1,ntj1) = c(nti1,ntj1) ...
                    + area * wq * ( qi * qj);
                a(nti1,ntj1) = a(nti1,ntj1) ...
                    + area * wq * ( dqidx * dqjdx + dqidy * dqjdy );
                
            end
            
        end
        
    end
    
end
% figure;
% surf(reshape(phi',nx,ny))
figure(1);
for it = 0:nt
    %   if rem (it * dt, .01) == 0
    imagesc(reshape(phi,nx,ny));
    title(['it=',num2str(it)]);
    colormap('gray');
    colorbar;
    pause(0.001)
    drawnow;
    %   end
    
    phi_dummy = reshape(phi,nx,ny);
    
    for j = 1 : ny
        
        for i = 1 : nx
            if (cond == 'D')
                % Dirichlet Boundary condition:
                if ( i == 1 || i == nx || j == 1 || j == ny )
                    
                    phi_dummy(i,j) = 0;
                    
                end
                
            elseif (cond == 'N')
                
                %Nuemann Boundary conditions:
                if ( i == 1 )
                    phi_dummy(i,j) = phi_dummy(i+1,j);
                elseif( i == nx )
                    phi_dummy(i,j) = phi_dummy(i-1,j);
                elseif( j == 1 )
                    phi_dummy(i,j) = phi_dummy(i,j+1);
                elseif(j == ny)
                    phi_dummy(i,j) = phi_dummy(i,j-1);
                end
                
            elseif (cond == 'PB')
                %Periodic Boundary condtions:
                if ( i == 1 )
                    phi_dummy(i,j) = phi_dummy(nx-1,j);
                elseif( i == nx )
                    phi_dummy(i,j) = phi_dummy(2,j);
                elseif( j == 1 )
                    phi_dummy(i,j) = phi_dummy(i,ny-1);
                elseif(j == ny)
                    phi_dummy(i,j) = phi_dummy(i,2);
                end
            elseif (cond == 'MX')
                if ( i == 1 )
                    phi_dummy(i,j) = 1;
                elseif( i == nx )
                    phi_dummy(i,j) = 0;
                elseif( j == 1 )
                    phi_dummy(i,j) = phi_dummy(i,j+1);
                elseif(j == ny)
                    phi_dummy(i,j) = phi_dummy(i,j-1);
                end
            else
                fprintf ( 1, 'Condition does not exist! \n' );
            end
        end
    end
    
    phi = reshape(phi_dummy,[],1);
    g_phi = calc_g ( node_num, node_xy,  ...
        element_order, element_num, element_node, phi ,par_a, par_b, par_c );
    b = -2 * a * phi + g_phi;
    
    dphidt = c \ b;
    nphi = dt * dphidt + phi;
    phi = nphi;
end
end
%
% calculate g(phi)=a phi^3+b phi^2+c phi
%
function g = calc_g ( node_num, node_xy, ...
    element_order, element_num, element_node, phi, par_a, par_b, par_c )
g = zeros(node_num,1);
for ie = 1 : element_num
    
    gi1 = element_node(1,ie);
    gi2 = element_node(2,ie);
    gi3 = element_node(3,ie);
    
    area = 0.5 * ...
        ( node_xy(1,gi1) * ( node_xy(2,gi2) - node_xy(2,gi3) ) ...
        + node_xy(1,gi2) * ( node_xy(2,gi3) - node_xy(2,gi1) ) ...
        + node_xy(1,gi3) * ( node_xy(2,gi1) - node_xy(2,gi2) ) );
    %
    %  Consider each quadrature point.
    %  Here, we use the midside nodes as quadrature points.
    %
    for q1 = 1 : 3
        
        q2 = mod ( q1, 3 ) + 1;
        
        nq1 = element_node(q1,ie);
        nq2 = element_node(q2,ie);
        
        xq = 0.5 * ( node_xy(1,nq1) + node_xy(1,nq2) );
        yq = 0.5 * ( node_xy(2,nq1) + node_xy(2,nq2) );
        wq = 1.0 / 3.0;
        
        phi_q = 0;
        %
        %  Consider each test function in the element.
        %
        for ti1 = 1 : element_order
            
            
            ti2 = mod ( ti1,     3 ) + 1;
            ti3 = mod ( ti1 + 1, 3 ) + 1;
            
            nti1 = element_node(ti1,ie);
            nti2 = element_node(ti2,ie);
            nti3 = element_node(ti3,ie);
            qi = 0.5 * ( ...
                ( node_xy(1,nti3) - node_xy(1,nti2) ) * ( yq - node_xy(2,nti2) ) ...
                - ( node_xy(2,nti3) - node_xy(2,nti2) ) * ( xq - node_xy(1,nti2) ) ) ...
                / area;
            
            phi_q = phi_q + qi * phi(nti1);
        end
        
        gvalue = par_a * phi_q.^3 + par_b * phi_q.^2 + par_c * phi_q;
        
        %
        %  Consider each basis function in the element.
        %
        for tj1 = 1 : element_order
            
            tj2 = mod ( tj1,     3 ) + 1;
            tj3 = mod ( tj1 + 1, 3 ) + 1;
            
            ntj1 = element_node(tj1,ie);
            ntj2 = element_node(tj2,ie);
            ntj3 = element_node(tj3,ie);
            
            qj = 0.5 * ( ...
                ( node_xy(1,ntj3) - node_xy(1,ntj2) ) * ( yq - node_xy(2,ntj2) ) ...
                - ( node_xy(2,ntj3) - node_xy(2,ntj2) ) * ( xq - node_xy(1,ntj2) ) ) ...
                / area;
            
            g(ntj1) = g(ntj1) + area * wq * gvalue * qj;
            
        end
        
    end
    
end
end
