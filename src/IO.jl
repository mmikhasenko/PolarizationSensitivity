
function parse_values_from_datafile_name(filename, Np)
    (Np!=6) && error("the parsing might not work: Np!=6 not implemeted!")
    v = match(r"Kstar=(.*),(.*)_Delta=(.*),(.*)_Lambda=(.*),(.*).txt", filename)
    eval.(Meta.parse.([v[i] for i in 1:Np]))
end