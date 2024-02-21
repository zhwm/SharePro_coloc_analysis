import pandas as pd

res = []
for lo in ['Locus1', 'Locus2', 'Locus3', 'Locus4', 'Locus5']:
    print(lo)
    for k in range(1, 6):
        for kc in range(1, k+1):
            ks = k - kc
            csnp = pd.read_csv('../sim/{}/{}_{}_0.05_0.01/csnp.txt'.format(lo, kc, ks), sep='\s+', header=None)
            for ite in range(1, 51):
                cs = pd.read_csv('../sim/{}/{}_{}_0.05_0.01/SH10/Q{}.z_G{}.z.cs'.format(lo, kc, ks, ite, ite), sep='\t')
                csco = cs.query('share>0.8')
                csset = [j for i in csco['cs'] for j in i.split('/')]
                csusie = pd.read_csv('../sim/{}/{}_{}_0.05_0.01/{}_coloc_susie.txt'.format(lo, kc, ks, ite), sep='\t')
                csuco = csusie.loc[csusie['PP.H4.abf']>0.8,]
                csuset = set(csuco['hit1']) | set(csuco['hit2'])
                ecaviar = pd.read_csv('../sim/{}/{}_{}_0.05_0.01/{}_col'.format(lo, kc, ks, ite), sep='\t')
                ecaco = ecaviar.query('CLPP>0.8')
                ecaset = [i for i in ecaco['SNP_ID']]
                res.append((lo, kc, ks, len([i for i in csnp.iloc[ite-1,ks:(kc+ks)] if i in csuset])/kc, len([i for i in csnp.iloc[ite-1,ks:(kc+ks)] if i in ecaset])/kc, len([i for i in csnp.iloc[ite-1,ks:(kc+ks)] if i in csset])/kc))

pd.DataFrame(res, columns=['Locus', 'KC', 'KS', 'COLOC+SuSiE', 'eCAVIAR', 'SharePro']).to_csv('../doc/sharepro_loc_sim_cs_power.csv', sep=',', index=False, header=True)
