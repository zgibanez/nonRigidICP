function R = getRotationMatrix(a,b,g)

ca = cos(a); cb = cos(b); cg = cos(g);
sa = sin(a); sb = sin(b); sg = sin(g);

R= [cg*cb, -sg*ca+cg*sb*sa, sg*sa+cg*sb*ca;
    sg*cb,  cg*ca+sg*sb*sa,  -cg*sa+sg*sb*ca;
    -sb,    cb*sa,           cb*ca];

end