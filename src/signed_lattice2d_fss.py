from LRGSG_package.LRGSG import *
os.chdir("../")

ssize_list = [2**k for k in range(4, 8)]
nrep = [range(2**(16-2*k)) for k in range(4, 8)]
pval_list = [0.01, 0.09, 0.098, 0.099, 0.100, 0.102, 0.104, 0.2]
print('list of sizes', ssize_list,
      '\nlist of flip probs', pval_list,
      '\nlist of replicas', nrep)

d_lmin = []

for i,L in enumerate(ssize_list):
    nedges = 2*L**2
    for p in pval_list:
        path = f"{datPath_lminl2d}N={L*L}/"
        savename = f"{path}p={p:.3g}.txt"
        if os.exists()
        lmin = []
        for nr in tqdm(nrep[i], desc=f"replicas for L={L}, p={p}"):
            G = nx.grid_2d_graph(L, L, periodic=True)
            ransample = random.sample(range(nedges), int(p*nedges))
            #
            all_weights = {e: 1 for e in G.edges()}
            neg_weights = {e: -1 for i,e in enumerate(G.edges()) if i in ransample}
            #
            nx.set_edge_attributes(G, values=all_weights, name='weight')
            nx.set_edge_attributes(G, values=neg_weights, name='weight')
            slapl = get_graph_lapl(G)
            try:
                lmin_tmp, _ = eigsh(slapl, k=1, which='SM', tol=1e-10)
            except ArpackNoConvergence:
                pass
            lmin.append(lmin_tmp)
        lminarr = np.array(lmin).flatten()
        d_lmin.append({'L': L, 'p':p, 'lmin': lminarr})
        if not os.path.isdir(path):
            os.makedirs(path)
        np.savetxt(savename, lminarr)