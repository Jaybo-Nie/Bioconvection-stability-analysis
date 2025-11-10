using Gridap
using Gridap.Geometry
using Gridap.Adaptivity
using DataStructures
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Algebra
using Base
using Printf
using Random
using Dates
using GridapPETSc

options = "-ksp_error_if_not_converged true -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"

GridapPETSc.with(args=split(options)) do 
    base     = @__DIR__ 
    save_dir = joinpath(base, "Solvers_single")
    mkpath(save_dir)

    total_start = time()
    
    # Parameters 
    #### Model Parameter ####
    L, H = 16, 2
    n = 40
    g = 980.665
    γ = 5e-10
    ν = 0.01
    θ = 0.01
    W = 0.05
    con = 1e4
    #### Time Discretization ####
    T, dt = 10, 0.01
    t = 0
    NN = Int(T/dt)
    #### Iteration Parameter ####
    tol = 1e-4
    max_subit = 10
    br = 1
    pr = 100

    # Discrete Model
    domain = (0,L,0,H)
    nC = (8*n, n)
    model = CartesianDiscreteModel(domain, nC; isperiodic=(true, false))

    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "top_bottom",["tag_5","tag_6"])


    # Approximation Space 
    #### Reference Finite Elements ####
    order = 2
    reffeV = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeP = ReferenceFE(lagrangian, Float64, order - 1)
    reffeC = ReferenceFE(lagrangian, Float64, order - 1)

    #### Spaces ####
    V = TestFESpace(model, reffeV; conformity=:H1, dirichlet_tags="top_bottom",dirichlet_masks=(false,true))
    U = TrialFESpace(V, VectorValue(0,0))

    Q = TestFESpace(model, reffeP, conformity=:H1, constraint=:zeromean)
    P = TrialFESpace(Q)

    VV = MultiFieldFESpace([V, Q])
    UU = MultiFieldFESpace([U, P])

    S = TestFESpace(model,reffeC;conformity=:H1)
    C = TrialFESpace(S)

    # Numerical Integration
    degree = order*2
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    

    # Marcos 
    ε(w) = 0.5 * (∇(w) + (∇(w))')
    sl2(u) = sqrt(sum( ∫( u*u )*dΩ ))
    vl2(u) = sqrt(sum( ∫( u ⋅ u )*dΩ ))
    conv(u,∇u) = (∇u')⋅u
    # Materialize linear combination in same space (fast)
    combine_interp(space, a, b, α, β) = interpolate(α*a + β*b, space)
    # L2 projection onto C
    function l2_project_C(f, Cspace, Sspace, dΩ)
        a(u,v) = ∫(u*v) * dΩ
        l(v) = ∫(f * v) * dΩ
        op = AffineFEOperator(a, l, Cspace, Sspace)
        return solve(op)
    end
    # Add Perturbation
    function perturb(solution, space; strength=0.1, base_seed=123)
        rng = Random.MersenneTwister(base_seed)
        n_dofs = num_free_dofs(space)
        solution_dofs = get_free_dof_values(solution)

        raw_noise = 2 .* rand(rng, n_dofs) .- 1
        zero_sum_noise = raw_noise .- mean(raw_noise)
   
        perturbation_dofs = strength * zero_sum_noise
    
        return FEFunction(space, solution_dofs + perturbation_dofs)
    end


    # Initialization
    uh0 = zero(U)
    uh = perturb(uh0, U; strength=0.0001, base_seed=929)
    uhold = uh
    uholder = uh
    uhk = uh
    uhkold = uh

    ph = zero(P)
    phold = zero(P)
    pholder = zero(P) 
    phk = zero(P) 
    phkold = zero(P) 

    initial_fun(x) = con * W * H / θ / (exp(W * H / θ) - 1) * exp(W * x[2] / θ)
    ch0 = l2_project_C(initial_fun, C, S, dΩ)
    ch = perturb(ch0, C; strength=0.01, base_seed=678)
    chold = ch
    cholder = ch
    chk = ch
    chkold = ch


    # Weak Form 
    aNSb((u,p), (v,q)) = 
        ∫((2*(u ⋅ v)/dt  - (∇⋅v)*p + q*(∇⋅u) + 2*ν*ε(u) ⊙ ε(v)))*dΩ

    c(u,v,uold) = ∫(v⊙(conv∘(uold,∇(u))))dΩ

    aNS((u,p),(v,q),uold) = aNSb((u,p),(v,q)) + c(u,v,uold)

    lNS((v,q), uold, ckold) = 
        ∫(-g*(1 + γ*ckold)*(v ⋅ VectorValue(0.0, 1.0)) + 2*(uold ⋅ v)/dt)*dΩ

    aCP(c, s, ukold) = 
        ∫((2*c*s/dt + (ukold ⋅ ∇(c))*s + θ*∇(c)⋅∇(s) - (VectorValue(0.0, 1.0) ⋅ ∇(s))*W*c))*dΩ

    lCP(s, cold) = ∫((2*cold*s/dt))*dΩ


    io = open(joinpath(save_dir, "DataSingleSolver.txt"), "a")
    @printf(io, "%s %s %s %s %s\n", "t", "totalmassIntegral", "KE", "timeDerU", "timeDerC")
    iotime = open(joinpath(save_dir, "TimeSingleSolver.txt"), "a")
    @printf(iotime, "%s %s %s %s\n", "t", "NS", "CP" ,"StepTime")

    ls = PETScLinearSolver()

    # Iteration 
    for step in 1:NN

        iter_start = time()

        t=step * dt
        cnt=1
        
        uhk = combine_interp(U, uhold, uholder, 1.5, -0.5)
        phk = combine_interp(P, phold, pholder, (1+dt), -dt)
        chk = combine_interp(C, chold, cholder, 1.5, -0.5)

        total_ns_solve = 0.0
        total_cp_solve = 0.0

        while cnt <= max_subit
            uhkold = uhk
            chkold = chk
            phkold = phk

            biformNS(x,y) = aNS(x,y, uhold)
            liformNS(y) = lNS(y,uhold,chkold)
            opNS = AffineFEOperator(biformNS,liformNS,UU,VV)
            t1 = time()
            xh = solve(ls,opNS)
            uhk, phk = xh
            total_ns_solve += (time() - t1)

            biformCP(x,y) = aCP(x,y,uhkold)
            liformCP(y) = lCP(y,chold)
            opCP = AffineFEOperator(biformCP,liformCP,C,S)
            t1 = time()
            chk = solve(opCP)
            total_cp_solve += (time() - t1)

            errU = vl2(uhk - uhkold)
            errC = sl2(chk - chkold)
            l2U = vl2(uhk)
            l2C = sl2(chk)
            errUl2 = errU/max(l2U,1e-8)

            cnt += 1

			if errUl2 < tol && errC < tol  
                break  
            end

        end

        uh = combine_interp(U, uhk, uhold, 2.0, -1.0)
        ph = interpolate(phk, P)
        ch = combine_interp(C, chk, chold, 2.0, -1.0)

        totalmassIntegral = sum( ∫(ch)*dΩ ) /H /L
        KE = 0.5 * sum( ∫(uh ⋅ uh)*dΩ )
        timeDerU = sqrt( sum( ∫((uh-uhold) ⋅ (uh-uhold) / dt / dt)*dΩ ) / sum( ∫(uh ⋅ uh / dt / dt)*dΩ ) )
		timeDerC = sqrt( sum( ∫((ch-chold) * (ch-chold) / dt / dt)*dΩ ) / sum( ∫(ch*ch/dt/dt)*dΩ ) )

        uholder, uhold = uhold, uh
        pholder, phold = phold, ph
        cholder, chold = chold, ch

        iter_time = time() - iter_start
        @printf(iotime, " %g %g %g %g\n", t, total_ns_solve, total_cp_solve, iter_time)
        flush(iotime)

        @printf(io, " %g %g %g %g %g\n", t, totalmassIntegral, KE, timeDerU, timeDerC)
        flush(io)

        if(mod(step,pr) ==0)
            writevtk(Ω, joinpath(save_dir, "results"*string(br)), cellfields=["u"=>uh, "p"=>ph, "c"=>ch])
            br=br+1
        end

    end

    close(io)
    close(iotime)

    total_time = time() - total_start
    println("Total time = $total_time")

end