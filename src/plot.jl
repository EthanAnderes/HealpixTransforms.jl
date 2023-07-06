
# alm plot 
# ------------------------


# TODO: add polarization
function alm_triangle_plot(
    Ifield::Xfield{<:ℍ0}; 
    imag_fun = x->abs2.(x),
    title = L"|a_{\ell m}|^2", 
    vmin  = nothing, vmax = nothing,
    xylabel_fontsize=7,
    title_fontsize=8, 
    ticklabel_fontsize=6, 
    ) 
    
    hpxℍ0 = fieldtransform(Ifield)
    lmax  = hpxℍ0.lmax
    # ls, ms = HT.lm(hpxℍ0.lmax) 
    # ls_mat = HT.alm2triangle(ls)
    # ms_mat = HT.alm2triangle(ms)
    l_ticks = m_ticks = 0:lmax

    alms_mat = alm2triangle(Ifield[!])
    nan_bool = isnan.(alms_mat)
    alms_mat[nan_bool] .= 0

    flm = alms_mat |> imag_fun
    flm[nan_bool] .= NaN

    fig, ax = subplots(1,dpi=147)
    img = ax.imshow(flm, cmap="viridis", 
        vmin=vmin, vmax=vmax, 
        # extent=[Hz_or_m(k[1]), Hz_or_m(k[end]), θ[end], θ[1]],
        # extent=[0, lmax, 0, lmax],
        extent=[0, lmax, lmax, 0],
        origin="upper"
    )
    ax.set_aspect("equal") 

    ax.set_xlabel(L"azmuthal frequency $m$",fontsize=xylabel_fontsize)
    ax.set_ylabel(L"\ell",fontsize=xylabel_fontsize)

    ax.tick_params(axis="both", labelsize=ticklabel_fontsize)

    ax.set_title(title, fontsize=title_fontsize)  
    fig.tight_layout()

    cbar1 = fig.colorbar(img, ax=ax, location="bottom", shrink = 0.5)
    cbar1.ax.tick_params(labelsize=ticklabel_fontsize)

    fig.tight_layout()

    return fig, ax
end


