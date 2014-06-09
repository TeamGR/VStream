function [ cosSim ] = similarity( FeatureVec1 , FeatureVec2  )
%SIMILARITY Computes the cosine similarity between two feature vectors
%   

cosSim = dot(FeatureVec1, FeatureVec2)/(norm(FeatureVec1,2)*norm(FeatureVec2,2));

end

