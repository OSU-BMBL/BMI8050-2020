# spatial modularity optimization 
cluster_sm <- function(G,D,maxit = 10, c_init = NULL)
{
    n = igraph::vcount(G) # number of nodes
    C_mat <- matrix(0,nrow = maxit+1,ncol = n)
    if(!is.null(c_init))
    {
        ci = c_init
    }
    else
    {
        ci = 1:n # initial clustering
    }
    C_mat[1,] <- ci
    # optimize modularity
    for(i in 1:maxit)
    {
        # loop through nodes
        for(j in 1:n)
        {
            # get clusters of each neighbor
            c_adj = unique(ci[igraph::neighbors(G,j)])
            n_adj = length(c_adj) # number of unique neighboring clusters
            # Matrix of proposed communities
            Z_prop = matrix(ci,
                            nrow = n_adj, 
                            ncol = n,
                            byrow = TRUE)
            Z_prop[,j] = c_adj # new proposed labels for j
            m_prop = rep(0,n_adj) # empty modularity storage
            for(l in 1:n_adj)
            {
                #m_prop[l] = spatial_modularity(G,D,Z_prop[l,])
                m_prop[l] = igraph::modularity(G,Z_prop[l,])
            }
            c_choose = c_adj[BiocGenerics::which.max(m_prop)]
            ci[j] = c_choose
        }
        print(igraph::modularity(G,ci))
        C_mat[i+1,] = ci
    }
    return(C_mat)
}