from functions import *
from lh_brute import *
import time
if __name__ == '__main__':
    params = load_params(sys.argv[1])
    name = str(pathlib.Path().absolute()) + "/data/real/" + params['dname']
    df = pd.read_csv(name)  # sys.argv[0] is the name of the python script
    df = df[df['mode'] != 'unknown']  # ignore the unknown edges
    df = df.reset_index()  # make sure indexes pair with number of rows
    genelist, genedict = initialize_graph(df)
    genelist.sort(key=lambda g: g.out_m+g.in_m, reverse=True)
    for ind,g in enumerate(genelist):
        g.ind = ind
        genedict[g.name] = ind  # updating the indices
    n = len(genelist)
    sample_eps = 1 / n
    start_time = time.time()
    Z_est, grad, _ = compute_likelihood(genelist, 'BNC', df=df, genedict=genedict, params=params, T=int(100 * n * log10(n)),
                       S =int(100 * n ** 2 / sample_eps ** 2), eps=0.01)
    print(Z_est)
    print(grad)
    end_time = time.time()
    print(end_time-start_time)
