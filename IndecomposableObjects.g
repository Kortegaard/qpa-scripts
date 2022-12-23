Read("./ARQuiver.g");

# Using constructing the AR quiver to find indecomposble objects.
AllIndecomposableModules := function(ind, steps)
    local lv, ARQ,v;

    ARQ := ConstructARQuiver(ind, steps);
    lv  := [];
    for v in ARQ[2] do
        Append(lv, [v.module]);
    od;
    return lv;
end;
