function [probVec] = SOFTMAX(vec)
% Apply softmax to vector to make it probability vector
    vec2 = exp(vec);
    try
        probVec = rshp(vec2./sum(vec2,2));
    catch
        probVec = vec2./sum(vec2,2);
    end
end

