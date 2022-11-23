function [planets_id,planets,N] = sequence_selector(marker)

switch marker
    case 1 % E - VEJ - S
        planets = {'Earth','Venus','Earth','Jupiter','Saturn'};
        planets_id = [3,2,3,5,6];
        N = length(planets_id) - 2;
    case 2 % E - VEE - S
        planets = {'Earth','Venus','Earth','Earth','Saturn'};
        planets_id = [3,2,3,3,6];
        N = length(planets_id) - 2;
    case 3 % E - VEVE - S
        planets = {'Earth','Venus','Earth','Venus','Earth','Saturn'};
        planets_id = [3,2,3,2,3,6];
        N = length(planets_id) - 2;
    case 4 % E - VEEJ - S
    planets = {'Earth','Venus','Earth','Earth','Jupiter','Saturn'};
    planets_id = [3,2,3,3,5,6];
    N = length(planets_id) - 2;
    case 5 % E - VVE - S
    planets = {'Earth','Venus','Venus','Earth','Saturn'};
    planets_id = [3,2,2,3,6];
    N = length(planets_id) - 2;
    case 6 % E - J - S
    planets = {'Earth','Jupiter','Saturn'};
    planets_id = [3,5,6];
    N = length(planets_id) - 2;
end

end