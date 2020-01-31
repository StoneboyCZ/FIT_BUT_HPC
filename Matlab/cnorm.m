%% cnorm
% Custom norm function used for results
function n = cnorm(v1,v2)
    n = norm(v1 - v2)/(norm(v2)+1);    
end