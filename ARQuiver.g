LoadPackage("QPA");

# Constructing the AR quiver, from:
#   ModuleList: [], a list of modules, to construct the AR based on
#   l: Int, Length for the predecessor function
#   NamingFunction: rec(module, numVertices) -> String, a function to name the vertices, given the module, and a rising number
# 
# Returns: Quiver, A part of the AR quiver.
ConstructARQuiverNamed := function(ModuleList, l, NamingFunction)
    local arrows, vertices, M, Q, v, a, r, w, d, rads,
    AddVertexIfNotIn, AddPredecessorsOfModule;

    arrows := [];
    vertices := [];

    # Check if a module is already added (up to isomorphism), and if it is not it will add the module to vertices.
    #
    # Returns tuple: (name: string, did_add: bool) 
    AddVertexIfNotIn := function(N)
        local v, v_name;
        v_name := false;

        # Going through the vertices, checking if it module is already added
        for v in vertices do
            if IsomorphicModules(N, v.module) then
                v_name := v.name;
            fi;
        od;

        # Deciding a name for vertex, based on the module
        if not IsString(v_name) then
            v_name := NamingFunction(rec(module := N, numVertices := Length(vertices)));
            # Adding vertex
            Append(vertices, [rec(
                name:= v_name,
                module:=N
            )]);
            return [v_name, true];
        fi;

        return [v_name, false];
    end;

    # Running the PredecessorsOfModule function for a given module 'M', with length 'l',
    # and knit the returned values into the quiver.
    AddPredecessorsOfModule := function(M, l)
        local P, dict, layer_num, layer, entry_num, index, v, a, arr, vertices_added, indec;
        vertices_added := 0;
        P := PredecessorsOfModule(M, l);

        # For keeping track of indecomposables names of (layer, entry) pairs, 
        # to know arrows start and termination points (from output P).
        dict := NewDictionary(false, true);

        # going throguh each of the indecomposable modules, of the P, adding them as vertices (if not already added),
        # and keep track of layer and entry, to add arrows later.
        layer_num := 1;
        for layer in P[1] do
            entry_num := 1;
            for indec in layer do
                v := AddVertexIfNotIn(indec);
                if v[2] then # if the vertex was just added
                    vertices_added := vertices_added + 1;
                fi;
                AddDictionary(dict, [layer_num,entry_num], v[1]);
                entry_num := entry_num + 1;
            od;
            layer_num := layer_num + 1;
        od;

        # Running through arrow part of output P, adding them to total list of arrows, if not already exists.
        layer_num := 1;
        for layer in P[2] do
            for arr in layer do
                a := [LookupDictionary(dict, [layer_num + 1, arr[1]]), # FROM
                      LookupDictionary(dict, [layer_num,     arr[2]]) # TO
                ];
                if not a in arrows then
                    Append(arrows, [a]);
                fi;
            od;
            layer_num := layer_num + 1;
        od;
        return vertices_added;
    end;



    # Running through each module given as input, and knitting their predecessors into the quiver.
    for M in ModuleList do
        if not IsProjectiveModule(M) then
            AddPredecessorsOfModule(M, l);
        fi;
    od;


    # Adding arrows to projectives, by decomposing radicals.
    for v in vertices do
        if IsProjectiveModule(v.module) then
            rads := RadicalOfModule(v.module);
            if Dimension(rads) = 0 then continue; fi;
            d := DecomposeModule(rads);
            for r in d do;
                w := AddVertexIfNotIn(r);
                a := [w[1], v.name];
                if not a in arrows then
                    Append(arrows, [a]);
                fi;
            od;
        fi;
    od;
    
    Q := Quiver(List(vertices, x -> x.name), arrows);
    return [Q, vertices];

end;

ConstructARQuiver := function(ModuleList, l)
    return ConstructARQuiverNamed(ModuleList, l, x -> String(DimensionVector(x.module)));
    #return ConstructARQuiverNamed(ModuleList, l, x -> Concatenation("v",String(x.numVertices)));
end;


