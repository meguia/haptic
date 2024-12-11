### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 376c4c9f-69c5-432d-8442-c350a4ed6931
begin
    import Pkg
    Pkg.add(; name = "Kaleido_jll", version = "0.1")
	Pkg.add("Symbolics")
	Pkg.add("PlotlyJS")
	Pkg.add("LinearAlgebra")
end

# ╔═╡ 504e5376-d2da-47da-8c02-cf970f74a784
Pkg.add("PlutoUI")

# ╔═╡ dea7bd10-b50a-11ef-1d3d-d165de9697d8
using  Symbolics, PlotlyJS, LinearAlgebra, PlutoUI

# ╔═╡ 7c62454e-efb7-4116-b37f-b74b4576bfa6
#parametric curve
@variables t,u,v;

# ╔═╡ eaae8960-8d27-452f-b9f5-c123fd0993ec
begin
	dt = Differential(t)
	du = Differential(u)
	dv = Differential(v)
end	

# ╔═╡ c7131a8b-8c36-4afc-baa0-7d4bc1fad7d8
md"""
x $(@bind fx Select([sin, cos])) $(@bind nx Select(collect(1:4))) $(string(fx)) ( $(nx) * t)\
y $(@bind fy Select([sin, cos])) $(@bind ny Select(collect(1:4))) $(string(fy)) ( $(ny) * t)\
z $(@bind fz Select([sin, cos])) $(@bind nz Select(collect(1:4))) $(string(fz)) ( $(nz) * t)\
N $(@bind N Slider(10:10:100,default=10;show_value=true))
"""

# ╔═╡ aea9c17c-7747-496e-a5d5-f513b1294f05
# parametric curve
f = [
	fx(nx*t),
	fy(ny*t),
	fz(nz*t)
]	

# ╔═╡ fa212423-eed7-49fd-a627-abb955eed195
md"""
# Parametric Surfaces
"""

# ╔═╡ af741731-1d0e-414d-a803-39823d5f7006
s = [
	(2+cos(u))*cos(v),
	(2+cos(u))*sin(v),
	sin(u)
]

# ╔═╡ ae8c48ca-3b8b-4331-b823-ba60c95ce1c6
function sym_evaluate(f,tlist::Array{Float64,1})
	# evaluate function f(t) at values tlist
	farr = [Symbolics.value(substitute(f, Dict(t=>tn))) for tn in tlist];
	x,y,z = [getindex.(farr,1),getindex.(farr,2),getindex.(farr,3)]
	return x, y, z
end

# ╔═╡ 122964c1-236b-4e15-925c-4b6be96868e0
function sym_evaluate(f,uvlist::Array{Array{Float64, 1}, 1})
	# evaluate function f(u,v) at values tlist
	farr = [Symbolics.value(substitute(f, Dict(u=>uv[1],v=>uv[2]))) for uv in uvlist];
	x,y,z = [getindex.(farr,1),getindex.(farr,2),getindex.(farr,3)]
	return x, y, z
end

# ╔═╡ 44980873-6c6d-48f9-81e0-ffb7744e2f79
function sym_evaluate(f,ulist,vlist)
	# evaluate function f(uv) at values ulist X vlist
	farr =[[Symbolics.value(substitute(f, Dict(u=>un,v=>vn))) for un in ulist] for vn in vlist];
	x,y,z = [getindex.(hcat(farr...),1),getindex.(hcat(farr...),2),getindex.(hcat(farr...),3)]
	return x, y, z
end

# ╔═╡ bc85882f-375c-4632-95d8-e6dbcbb80678
function arrow3d!(x, y, z,  u, v, w; as=0.1, lc=:black, lw=3.0, scale=1.0)
	arrows = [scatter3d(x=x,y=y,z=z,mode="markers",color=lc,showlegend = false, marker_size=3)]
    (as < 0) && (nv0 = -maximum(norm.(eachrow([u v w]))))
    for (x,y,z, u1,v1,w1) in zip(x,y,z, u,v,w)
		u,v,w = [u1,v1,w1]/scale
        nv = sqrt(u^2 + v^2 + w^2) # modulus of vector
        v1, v2 = -[u,v,w]/nv, nullspace(adjoint([u,v,w]))[:,1]
        v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
        v5 = v4 - 2*(v4'*v2)*v2
        (as < 0) && (nv = nv0) 
        v4, v5 = -as*nv*v4, -as*nv*v5
        l1 = scatter3d(x=[x,x+u], y=[y,y+v], z=[z,z+w], mode="lines",line_width=lw, line_color=lc, showlegend = false)
        l2 = scatter3d(x=[x+u,x+u-v5[1]], y=[y+v,y+v-v5[2]], z=[z+w,z+w-v5[3]],  mode="lines", line_width=lw, line_color=lc,showlegend = false)
        l3 = scatter3d(x=[x+u,x+u-v4[1]], y=[y+v,y+v-v4[2]], z=[z+w,z+w-v4[3]],  mode="lines", line_width=lw, line_color=lc,showlegend = false)
		append!(arrows,[l1,l2,l3])
    end
	return arrows
end

# ╔═╡ 06768653-5505-4a1a-bb4c-13d0cd47bbd0
function parametric_surface(f,umin,umax,vmin,vmax,Nu,Nv,uvlist; lc=:black, scale=3.0, lcu=:blue, lcv=:red, lcn=:green, lw=4.0)
	# returns a plot of a parametric surface given by f(u,v) in the range [umin,umax]X[vmin,vmax] with NuxNv points and # the tangent vectors evaluated at values (u,v) of the parameter uvlist
	x, y, z = sym_evaluate(f,umin:(umax-umin)/Nu:umax,vmin:(vmax-vmin)/Nv:vmax)
	x1, y1, z1 = sym_evaluate(f,uvlist)
	dfdu = [simplify(expand_derivatives(du(fn))) for fn in f]
	dfdv = [simplify(expand_derivatives(dv(fn))) for fn in f]
	dfn = cross(dfdv,dfdu)/norm(cross(dfdv,dfdu))
	dfdu /= norm(dfdu)
	dfdv /= norm(dfdv)
	dxu, dyu, dzu = sym_evaluate(dfdu,uvlist)
	dxv, dyv, dzv = sym_evaluate(dfdv,uvlist)
	dxn, dyn, dzn = sym_evaluate(dfn,uvlist)
	#s1 = surface(x=x,y=y,z=z,opacity=0.7,surfacecolor=@. x^2+y^2+z^2,showscale=false)
	s1 = surface(x=x,y=y,z=z,showscale=false, opacity=0.8,surfacecolor=@. x)
	v1 = arrow3d!(x1,y1,z1,dxu,dyu,dzu;lc=lcu,lw=lw,scale=scale)
	v2 = arrow3d!(x1,y1,z1,dxv,dyv,dzv;lc=lcv,lw=lw,scale=scale)
	v3 = arrow3d!(x1,y1,z1,dxn,dyn,dzn;lc=lcn,lw=lw,scale=scale)
	return(append!(v1,v2,v3,[s1]))
end

# ╔═╡ 2911049c-c9a1-49f3-ac39-dd97dd10151d
begin
	go2 = parametric_surface(s,0,2*pi,0,2*pi,30,30,[[pi/3,pi/3],[0,pi/3]]; scale=2.0,lw=8.0)
	layout2 = Layout(title=join(string.(s)," , "), autosize=false, width=800, height=800)
	p2 = Plot(go2,layout2)
	p2
end	

# ╔═╡ 0583a5cf-1996-468c-b27f-fbadefb88f76
function parametric_curve(f,tmin,tmax,N,tlist; lc=:blue, scale=3.0, lca=:black, lw=0.4)
	# returns a plot of a parametric curve given by f(t) between tmin and tmax with N points and the 
	# tangent vector evaluated at values of the parameter tlist
	x, y, z = sym_evaluate(f,collect(tmin:(tmax-tmin)/N:tmax))
	x1, y1, z1 = sym_evaluate(f,collect(tlist))
	dfdt = [simplify(expand_derivatives(dt(fn))) for fn in f]
	dx, dy, dz = sym_evaluate(dfdt,collect(tlist))
	s1 = scatter3d(x=x,y=y,z=z,mode="lines",line_color=lc,showlegend=false)
	v1 = arrow3d!(x1,y1,z1,dx,dy,dz;scale=scale)
	return(append!(v1,[s1]))
end

# ╔═╡ f7c9f82c-e26a-43d7-968c-528d9d1af68f
begin
	go1 = parametric_curve(f,0,2*pi,300,0:2*pi/N:2*pi;scale=5.0)
	layout = Layout(title=join(string.(f)," , "), autosize=false, width=800, height=800)
	p1 = Plot(go1,layout)
	p1
end	

# ╔═╡ 6e39c24f-a52b-432f-95be-b3dc47ef8cb1
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 95%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═376c4c9f-69c5-432d-8442-c350a4ed6931
# ╠═504e5376-d2da-47da-8c02-cf970f74a784
# ╠═dea7bd10-b50a-11ef-1d3d-d165de9697d8
# ╠═7c62454e-efb7-4116-b37f-b74b4576bfa6
# ╠═eaae8960-8d27-452f-b9f5-c123fd0993ec
# ╟─f7c9f82c-e26a-43d7-968c-528d9d1af68f
# ╟─c7131a8b-8c36-4afc-baa0-7d4bc1fad7d8
# ╟─aea9c17c-7747-496e-a5d5-f513b1294f05
# ╟─fa212423-eed7-49fd-a627-abb955eed195
# ╠═af741731-1d0e-414d-a803-39823d5f7006
# ╠═2911049c-c9a1-49f3-ac39-dd97dd10151d
# ╠═06768653-5505-4a1a-bb4c-13d0cd47bbd0
# ╟─ae8c48ca-3b8b-4331-b823-ba60c95ce1c6
# ╟─122964c1-236b-4e15-925c-4b6be96868e0
# ╟─44980873-6c6d-48f9-81e0-ffb7744e2f79
# ╟─0583a5cf-1996-468c-b27f-fbadefb88f76
# ╠═bc85882f-375c-4632-95d8-e6dbcbb80678
# ╠═6e39c24f-a52b-432f-95be-b3dc47ef8cb1
