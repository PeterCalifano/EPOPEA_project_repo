function [planets_id,planets] = sequence_selector(marker)

switch marker
    case 1 % E - VEJ - S
        planets = {'Earth','Venus','Earth','Jupiter','Saturn'};
        planets_id = [3,2,3,5,6];
    case 2 % E - VEE - S
        planets = {'Earth','Venus','Earth','Earth','Saturn'};
        planets_id = [3,2,3,3,6];
    case 3 % E - VEVE - S
        planets = {'Earth','Venus','Earth','Venus','Earth','Saturn'};
        planets_id = [3,2,3,2,3,6];
    case 4 % E - VEEJ - S
    planets = {'Earth','Venus','Earth','Earth','Jupiter','Saturn'};
    planets_id = [3,2,3,3,5,6];
end

end