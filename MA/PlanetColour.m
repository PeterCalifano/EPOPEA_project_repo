function [ColorHexCode] = PlanetColour(inputName)

% BODYCOLORS = {'Sun': '#ffcc00',...
%     'Mercury': '#8c8680',...
%     'Venus': '#e6db67',...
%     'Earth': '#2a7bd1',...
%     'Moon': '#999999',...
%     'Mars': '#cc653f',...
%     'Jupiter': '#bf8f5c',...
%     'Saturn': '#decf83',...
%     'Uranus': '#7ebec2',...
%     'Neptune': '#3b66d4'};



planets_names = {'SUN', 'MERCURY', 'VENUS', 'EARTH', 'MARS', 'JUPITER', 'SATURN', 'URANUS', 'NEPTUNE'};

planet_colours = {'#ffcc00', '#8c8680', '#e6db67', '#2a7bd1', '#cc653f', '#bf8f5c', '#decf83', '#7ebec2', '#3b66d4'};

try 
    if isnumeric(inputName)
        planets_numbers = 0:(length(planets_names)-1);
        code_number = (planets_numbers == inputName);

        inputName = planets_names{code_number};
    end
    BODYCOLOURS = containers.Map(planets_names, planet_colours);
    ColorHexCode = BODYCOLOURS(upper(inputName));
catch
    warning('Body not in list. Exception output: #000000')
    ColorHexCode = '#000000';
end

end