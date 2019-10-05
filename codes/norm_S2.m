function S2 = norm_S2(E, k,r, alpha)

% normalization 
alpha_E = E*alpha; 

% compute KE
S2temp1 = S2_from_E(E, k,r); 

% compute KKE 
S2temp2 = S2_from_E(S2temp1, k ,r); 

S2 = S2temp2 + alpha_E; 


end