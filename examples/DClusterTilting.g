
Read("../DClusterTilting.g");

## EXAMPLES

A := NakayamaAlgebra([2,2,2], GF(3));;
In := AllIndecModulesOfLengthAtMost(A, 3);; 

M := DirectSumOfQPAModules([In[4], In[5], In[6]]);; 
Display(Concatenation("Is P = P1 + P2 + P3 3-cluster tilting: ", String(IsDClusterTilting(M,3,In))));

M2 := DirectSumOfQPAModules([In[1], In[4], In[5], In[6]]);; 
Display(Concatenation("Is P = S1 + P1 + P2 + P3 3-cluster tilting: ", String(IsDClusterTilting(M2,3,In))));

M3 := DirectSumOfQPAModules([In[1], In[4], In[5], In[6]]);; 
Display(Concatenation("Is P = S1 + P1 + P2 + P3 3-cluster tilting: ", String(IsDClusterTilting(M3,3,In))));

