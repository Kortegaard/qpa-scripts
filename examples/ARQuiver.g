
Read("../ARQuiver.g");


#TEST
Q := DynkinQuiver("A", 3, ["r","l"]);
A := PathAlgebra(GF(3), Q);
I := IndecInjectiveModules(A);

AR_Q := ConstructARQuiver([I[1], I[3], I[2]],4);

Display("Q := 1 --> 2 <-- 3");
Display(Concatenation("AR Quiver based from: ", String(DimensionVector(I[1])), ", ",String(DimensionVector(I[3]))));
Display("AR QUIVER: \n--------------------");
Display(AR_Q[1]);
Display("--------------------\n\n");

# Second Quiver
Q2 := Quiver(4,[[1,2],[2,3],[3,4],[4,1]]);
kQ2 := PathAlgebra(GF(3), Q2);
R := NthPowerOfArrowIdeal(kQ2,3);
A2 := kQ2 / R;
S2 := SimpleModules(A2);
AR_Q2 := ConstructARQuiver([S2[1]],10);
#AR_q3 := ConstructARQuiverNamed(ModuleList, l, x -> x);


Display("Q := 1 --> 2 --> 3 --> 4 --> 1");
Display("ideal := {arrow ideal}^3 = {ideal of paths of length 3}");
Display(Concatenation("AR Quiver based from: ", String(DimensionVector(S2[1]))));
Display("AR QUIVER: \n--------------------");
Display(AR_Q2[2]);
Display("--------------------");



