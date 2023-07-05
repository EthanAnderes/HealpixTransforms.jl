
# alm plot 
# ------------------------


# TODO: add polarization
function alm_triangle_plot(
    Ifield::Xfield{<:ℍ0}; 
    imag_fun = x->abs2.(x),
    title = L"|a_{\ell m}|^2", 
    vmin  = nothing, vmax = nothing,
    ) 
    # alms_mat; # pass alms_mat = HT.alm2triangle(alms_hpx)
    # l_ticks, m_ticks,
    # logs = false, blur = 0, 
    # vmin=nothing, vmax=nothing, 
    # title=nothing, 
    # )

    
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


    # img1 = ax.imshow(flm, vmin=vmin, vmax=vmax, origin="upper")
    # ax.set_aspect("auto") # , adjustable="box")
    # ax.set_xlim(minimum(m_ticks), maximum(m_ticks))
    # ax.set_ylim(maximum(l_ticks), minimum(l_ticks))
    ax.set_xlabel(L"azmuthal frequency $m$",fontsize=7)
    ax.set_ylabel(L"\ell",fontsize=7)

    ax.tick_params(axis="both", labelsize=6)

    ax.set_title(title, fontsize=8)  
    fig.tight_layout()

    cbar1 = fig.colorbar(img, ax=ax, location="bottom", shrink = 0.5)
    cbar1.ax.tick_params(labelsize=6)

    fig.tight_layout()

    return fig, ax
end


