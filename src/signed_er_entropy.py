from LRGSG_package.LRGSG import *

no_p = 10
no_rw = 5
pow2_m = 7
pow2_M = 13

lsN = np.array([2**i for i in range(pow2_m, pow2_M)])
lsNavg = 2 * np.array([2**(pow2_M-i) for i in range(pow2_m, pow2_M)])
lsp = np.array([(1+i)/(no_p) for i in range(no_p-1)])
lsrw = np.array([(1+i)/(no_rw) for i in range(no_rw-1)])

plt.rcParams.update({'font.size': 30})

for N, no_avg in zip(lsN, lsNavg):
    for p in lsp:
        print(N, p)
        for r in lsrw:
            entropy_full = []
            Cspec_full = []
            fig, ax = plt.subplots(figsize=(12, 9))
            ax2 = ax.twinx()
            ax.set_ylabel(r'$1-S$')
            ax2.set_ylabel(r'$C\log N$')
            ax.set_xlabel(r'$\tau$')
            ax.set_xscale('log')
            for nvg in tqdm(range(no_avg)):
                
                G = nx.erdos_renyi_graph(N, p, seed=int(N*p), directed=False)
                for e in G.edges():
                    G.add_edge(e[0], e[1], weight=1)
                    if np.random.random() < r:
                        G[e[0]][e[1]]['weight'] = -1
                [Sm1, dS1, VarL1, t1] = entropy(G, is_signed=True)
                dS1 /= np.max(dS1)
                t11 = (t1[1:] + t1[:-1]) / 2.
                entropy_full.append([t1, Sm1])
                Cspec_full.append([t11, dS1])
                pSm1, = ax.plot(t1, Sm1, '-')
                ax2.plot(t11, dS1, '--', c=pSm1.get_color())
            fig.savefig(f"{pltPath_Sm1C}ErdosReny_N={N}_p={p:.3g}_r={r:.3g}_navg={no_avg}{ePDF}")
            plt.close('all')
            fig, ax = plt.subplots(figsize=(12, 9))
            ax2 = ax.twinx()
            ax.set_ylabel(r'$1-S$')
            ax2.set_ylabel(r'$C\log N$')
            ax.set_xlabel(r'$\tau$')
            ax.set_xscale('log')
            tSm1, avgdSm1 = averaged_Sm1(entropy_full)
            tCs, avgdCspec = averaged_Sm1(Cspec_full)
            pSm1, = ax.plot(tSm1, avgdSm1, '-')
            ax.plot(tCs, avgdCspec, '--', c=pSm1.get_color())
            fig.savefig(f"{pltPath_Sm1C}ErdosReny_N={N}_p={p:.3g}_r={r:.3g}_avgd{ePDF}")
            plt.close('all')
