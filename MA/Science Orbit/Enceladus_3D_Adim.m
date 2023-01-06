function Enceladus_3D_Adim(Rs, r)

% Enceladus_3D_Adim.m - Enceladus texture loaded in a plot, non-dimensional
% axes
%
% PROTOTYPE:
%   Enceladus_3D_Adim(Rs,r)
%
% DESCRIPTION:
%   Function to load the Sun modelled as a sphere inside a figure.
%
% INPUT:
%   Rs          [1x1]       Enceladus mean radius
%   r           [1x3]       location of the Enceladus
%
% OUTPUT:
%   []          [figure]    Figure open with the Sun picture loaded
%
% ------------------------------------------------------------------------


%%  Load the Earth image from a website
Enceladus_image = 'Enceladus_surface.jpg';

%% Figure

% Choose the color of the figure background
background_plot = 'w';

% Create the figure
figure('Color', background_plot);
hold on;
grid on;

% Set the axes scale equal
axis equal;

% Put the axes labels
xlabel('X [DU]');
ylabel('Y [DU]');
zlabel('Z [DU]');

% Set initial view
view(120,30);

%% Create Earth surface as a wireframe

% Define the number of panels to be used to model the sphere 
npanels = 180;  

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(r(1), r(2), r(3), Rs, Rs, Rs, npanels);

% Create the globe with the surf function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

%% Texturemap the globe

% Load Earth image for texture map
cdata = imread(Enceladus_image);

% Set the transparency of the globe: 1 = opaque, 0 = invisible
alpha = 1; 

% Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
% specify the image data using the 'CData' property with the data loaded 
% from the image. Finally, set the transparency and remove the edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

end
