function proj_ab = proj(a,b)
    proj_ab = (dot(a,b)/norm(a))*(a/norm(a));
end