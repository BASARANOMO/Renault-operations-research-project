function create_graph(G_original, dims, us, fs, J_init, J_fin, e)
    U = dims.U
    F = dims.F
    J = J_fin - J_init + 1
    # J = dims.J
    cost_dispatch = dims.γ
    
    G = SimpleDiGraph((3 * U + 4 * F) * J + U + F + 1)
    capacity = zeros((3 * U + 4 * F) * J + U + F + 1, 
        (3 * U + 4 * F) * J + U + F + 1)
    cost_mat = zeros((3 * U + 4 * F) * J + U + F + 1, 
        (3 * U + 4 * F) * J + U + F + 1)
    
    for j = 1:J # for each day
        for u = 1:U # connect each u to each f stock
            for f = 1:F
                add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                    (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1)
                capacity[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                    (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1] += Inf
                cost_mat[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                    (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1] += G_original.d[u, U + f] * cost_dispatch
            end
        end
        
        for f = 1:F
            # connect fictive u (carton) to each f consommation
            add_edge!(G, nv(G), (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f)
            #add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f, nv(G)) # cycle
            capacity[nv(G), (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f] += Inf
            #capacity[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f, nv(G)] += Inf
            cost_mat[nv(G), (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f] += fs[f].cexc[e]
            # cost_mat[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f, nv(G)] += fs[f].cexc[e]
        end
    end
    
    for j = 1:J
        for u = 1:U
            # connect u day j to two additional vertices
            add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 2)
            add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 3)
            
            capacity[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 2] += us[u].r[e, j + J_init - 1]
            capacity[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 3] += Inf
            
            cost_mat[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 3] += us[u].cs[e]
            
            if j < J
                # connect two additional vertices to u day j+1
                add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 2, 
                    j * (3 * U + 4 * F) + 3 * (u - 1) + 1)
                add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 3, 
                    j * (3 * U + 4 * F) + 3 * (u - 1) + 1)
                
                capacity[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 2, 
                    j * (3 * U + 4 * F) + 3 * (u - 1) + 1] += us[u].r[e, j + J_init - 1]
                capacity[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 3, 
                    j * (3 * U + 4 * F) + 3 * (u - 1) + 1] += Inf
            end
        end
        
        for f = 1:F
            # connect f stock day j to two additional vertices
            add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 2)
            add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 3)
            
            capacity[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 2] += fs[f].r[e, j + J_init - 1]
            capacity[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 3] += Inf
            
            cost_mat[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1, 
                (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 3] += fs[f].cs[e]
            
            if j < J
                # connect each f stock to f consommation of the next day
                add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f, j * (3 * U + 4 * F) + 3 * U + 4 * f)
                capacity[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f, j * (3 * U + 4 * F) + 3 * U + 4 * f] += Inf
                
                # connect two additional vertices to f day j+1
                add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 2, 
                    j * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1)
                add_edge!(G, (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 3, 
                    j * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1)
                
                capacity[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 2, 
                    j * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1] += fs[f].r[e, j + J_init - 1]
                capacity[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 3, 
                    j * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1] += Inf
            end
        end
    end
    
    # connect initial stock to the stock of day 1 (for each u)
    for u = 1:U
        add_edge!(G, (3 * U + 4 * F) * J + u, 3 * (u - 1) + 1)
        
        capacity[(3 * U + 4 * F) * J + u, 3 * (u - 1) + 1] += Inf
    end
    
    # connect initial stock to the stock of day 1 (for each f) & to the f consommation of day 1
    for f = 1:F
        add_edge!(G, (3 * U + 4 * F) * J + U + f, 3 * U + 4 * (f - 1) + 1)  
        add_edge!(G, (3 * U + 4 * F) * J + U + f, 3 * U + 4 * f)
        
        capacity[(3 * U + 4 * F) * J + U + f, 3 * U + 4 * (f - 1) + 1] += Inf
        capacity[(3 * U + 4 * F) * J + U + f, 3 * U + 4 * f] += Inf
    end
    
    return G, capacity, cost_mat
end

# for a given e
function set_demand(dims, us, fs, J_init, J_fin, e) # b parameter
    U = dims.U
    F = dims.F
    J = J_fin - J_init + 1
    # J = dims.J
    demand = zeros((3 * U + 4 * F) * J + U + F + 1)
    
    for j = 1:J
        # for each u
        for u = 1:U
            demand[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1] -= us[u].b⁺[e, j + J_init - 1]
        end
        
        # for each f consommation
        for f = 1:F
            demand[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * f] += fs[f].b⁻[e, j + J_init - 1]
        end
    end
    
    # for initial stock
    for u = 1:U
        demand[J * (3 * U + 4 * F) + u] -= us[u].s0[e]
    end
    
    for f = 1:F
        demand[J * (3 * U + 4 * F) + U + f] -= fs[f].s0[e]
    end
    
    println(sum(demand))
    # additional vertices from day 1 to day J-1: b = 0
    # additional vertices of day J: b >= 0
    demand[(3 * U + 4 * F) * J + U + F + 1] += Inf # +Inf means that b is not fixed
    for u = 1:U
        demand[(J - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 2] -= Inf # -Inf means b>=0
        demand[(J - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 3] -= Inf
    end
    for f = 1:F
        demand[(J - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 2] -= Inf
        demand[(J - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 3] -= Inf
    end
    
    return demand
end

function min_cost_flow(g, node_demand, edge_capacity, edge_cost, optimizer)
    m = JuMP.Model(optimizer)
    vtxs = vertices(g)
    
    source_nodes = [v for v in vtxs if node_demand[v] < 0 && node_demand[v] != -Inf]
    sink_nodes = [v for v in vtxs if node_demand[v] > 0 && node_demand[v] != Inf]
    
    idx_dict = Dict()
    ridx_dict = Dict()
    i = 1
    for e in lg.edges(g)
        idx_dict[i] = [src(e), dst(e)]
        ridx_dict[(src(e), dst(e))] = i
        i += 1
    end
    
    @variable(m, 0 <= f[i = 1:ne(g)] <= edge_capacity[idx_dict[i][1], idx_dict[i][2]])
    @objective(m, Min, sum(f[i] * edge_cost[idx_dict[i][1], idx_dict[i][2]] for i = 1:ne(g)))
    # @variable(m, 0 <= f[i=vtxs, j=vtxs; (i,j) in lg.edges(g)] <= edge_capacity[i, j])
    # @objective(m, Min, sum(f[src(e),dst(e)] * edge_cost[src(e), dst(e)] for e in lg.edges(g)))

    for v in lg.vertices(g)
        if v in source_nodes
            @constraint(m,
                sum(f[ridx_dict[(v, vout)]] for vout in outneighbors(g, v)) - sum(f[ridx_dict[(vin, v)]] for vin in lg.inneighbors(g, v)) == -node_demand[v]
            )
        elseif v in sink_nodes
            @constraint(m,
                sum(f[ridx_dict[(vin, v)]] for vin in lg.inneighbors(g, v)) - sum(f[ridx_dict[(v, vout)]] for vout in outneighbors(g, v)) == node_demand[v]
            )
        else
            if node_demand[v] == -Inf
                @constraint(m, sum(f[ridx_dict[(vin, v)]] for vin in lg.inneighbors(g, v)) - sum(f[ridx_dict[(v, vout)]] for vout in outneighbors(g, v)) >= 0)
            elseif node_demand[v] != Inf
                @constraint(m,
                    sum(f[ridx_dict[(vin, v)]] for vin in lg.inneighbors(g, v)) == sum(f[ridx_dict[(v, vout)]] for vout in outneighbors(g, v))
                )
            end
        end
    end

    optimize!(m)
    ts = termination_status(m)
    result_flow = spzeros(nv(g), nv(g))
    if ts != MOI.OPTIMAL
        @warn "Problem does not have an optimal solution, status: $(ts)"
        return result_flow
    end
    for e in lg.edges(g)
        (i,j) = Tuple(e)
        result_flow[i,j] = JuMP.value(f[ridx_dict[(i,j)]])
    end
    return result_flow
    
end

function run_opt(g, dims, us, fs, J_init, J_fin, e, optimizer)
    g, capacity, cost_mat = create_graph(g, dims, us, fs, J_init, J_fin, e)
    demand = set_demand(dims, us, fs, J_init, J_fin, e)
    flow = min_cost_flow(g, demand, capacity, cost_mat, optimizer)
    return flow
end

function read_flow(flow, U, F, J_init, J_fin)
    # create a U * F * J matrix (for a given e) for dispatching
    # create a U * J matrix (for a given e) for stockage
    # create a F * J matrix (for a given e) for stockage
    J = J_fin - J_init + 1
    dispatch = zeros(U, F, J)
    stock_U = zeros(U, J)
    stock_F = zeros(F, J)
    for j = 1:J
        for u = 1:U
            for f = 1:F
                dispatch[u, f, j] += flow[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1,
                    (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1]
            end
        end
        
        for u = 1:U
            stock_U[u, j] = flow[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 2] +
            flow[(j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 1, (j - 1) * (3 * U + 4 * F) + 3 * (u - 1) + 3]
        end
        
        for f = 1:F
            stock_F[f, j] = flow[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1,
                (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 2] + 
            flow[(j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 1,
                (j - 1) * (3 * U + 4 * F) + 3 * U + 4 * (f - 1) + 3]
        end
    end
    
    return dispatch, stock_U, stock_F
end