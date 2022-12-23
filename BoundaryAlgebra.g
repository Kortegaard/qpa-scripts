LoadPackage("QPA");

# Given k,n this will return the correspoding boundary algebra B_{k,n}.
#
# Parameters:
#   k - Int
#   n - Int, n>k.
#
# Returns: 
#   Algebra

BoundaryAlgebra := function(k, n)
    local i,j, arrows, temp, temp2, ideal, kQ, Q, I, A, getX, getY;

    arrows := [];
    ideal  := [];

    # Generate arrows
    for i in  [1..n] do
        temp := i+1; if temp > n then temp :=1; fi;
        Append(arrows, [
            [i, temp, Concatenation("x",String(i))],
            [temp, i, Concatenation("y",String(i))]
         ]);
    od;
    Q := Quiver(n, arrows);
    kQ := PathAlgebra(GF(3), Q);

    # Return x,y correct generators given number. (0 indexed, i.e. return x(num+1))
    getX := function(num) return GeneratorsOfAlgebra(kQ)[n+num*2+1]; end;
    getY := function(num) return GeneratorsOfAlgebra(kQ)[n+num*2+2]; end;

    # Generate
    for i in  [1..n] do
        Append(ideal,[getX(i-1)*getY(i-1)]);

        temp  := getX(\mod(i,n));
        temp2 := getY(\mod(i-1,n));
        for j in  [1..k-1] do temp  := temp * getX(\mod(i+j,n)); od;
        for j in  [1..n-k-1] do temp2 := temp2 * getY(\mod(i-j-1,n)); od;
        Append(ideal,[temp - temp2]);
    od;

    I := Ideal(kQ, ideal);
    A := kQ/I;
    return A;
end;

