using Symbolics
using PlotlyJS
using LinearAlgebra


function sym_evaluate(f,tlist::Array{Float64,1})
	# evaluate function f(t) at values tlist
	farr = [Symbolics.value(substitute(f, Dict(t=>tn))) for tn in tlist];
	x,y,z = [getindex.(farr,1),getindex.(farr,2),getindex.(farr,3)]
	return x, y, z
end


function sym_evaluate(f,uvlist::Array{Array{Float64, 1}, 1})
	# evaluate function f(u,v) at values tlist
	farr = [Symbolics.value(substitute(f, Dict(u=>uv[1],v=>uv[2]))) for uv in uvlist];
	x,y,z = [getindex.(farr,1),getindex.(farr,2),getindex.(farr,3)]
	return x, y, z
end

function sym_evaluate(f,ulist,vlist)
	# evaluate function f(uv) at values ulist X vlist
	farr =[[Symbolics.value(substitute(f, Dict(u=>un,v=>vn))) for un in ulist] for vn in vlist];
	x,y,z = [getindex.(hcat(farr...),1),getindex.(hcat(farr...),2),getindex.(hcat(farr...),3)]
	return x, y, z
end


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