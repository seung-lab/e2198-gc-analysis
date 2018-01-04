function patch_cut = cutedge(patch,cutX,cutY)
    patch_cut = patch(cutX+1:end-cutX,cutY+1:end-cutY);
end