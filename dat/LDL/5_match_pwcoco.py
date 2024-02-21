import gzip
import pandas as pd

def get_bim_dict(bimdir):
    # get bim
    bim_dict = {}
    bim_dup = []
    with open(bimdir) as bim:
        for line in bim:
            ll = line.strip().split()
            if ll[1] in bim_dict:
                bim_dup.append(ll[1])
            else:
                bim_dict[ll[1]]=(ll[1], ll[4], ll[5])
    # remove dups
    for drs in bim_dup:
        del bim_dict[drs]
    return bim_dict

def get_rsid_dict(rsdir):
    # convert to rsid
    rsid_dict = {}
    with open(rsdir) as rs:
        for line in rs:
            ll = line.strip().split()
            rsid_dict[ll[0]]=ll[1]
    return rsid_dict

def get_ss_dict(ssdir, rsid_dict, bim_dict):
    # match gwas
    ss_dict = {}
    with gzip.open(ssdir, 'rt') as ss:
        header = next(ss)
        for line in ss:
            ll = line.strip().split()
            rsid = rsid_dict.get(ll[0].replace(':', '.'))
            alt = ll[0].split(':')[2]
            ref = ll[0].split(':')[3]
            if ll[1]==ref:
               maf = round(float(ll[2]),4)
            else:
               maf = round(1-float(ll[2]),4)
            if rsid in bim_dict:
                if bim_dict.get(rsid)==(rsid,ref,alt):
                    ss_dict[rsid]= bim_dict[rsid] + (maf, round(float(ll[7]),4), round(float(ll[8]), 4), ll[10], ll[4])
                elif bim_dict.get(rsid)==(rsid,alt,ref):
                    ss_dict[rsid]= bim_dict[rsid] + (1-maf, round(-float(ll[7]),4), round(float(ll[8]), 4), ll[10], ll[4])
    return ss_dict

def get_qtl_dict(qtldir, bim_dict):
    # match qtl
    qtl_dict = {}
    with gzip.open(qtldir, 'rt') as ss:
        header = next(ss)
        for line in ss:
            ll = line.strip().split()
            rsid = ll[0]
            alt = ll[3].upper()
            ref = ll[2].upper()
            if rsid in bim_dict:
                if bim_dict.get(rsid)==(rsid,ref,alt):
                    qtl_dict[rsid]= bim_dict[rsid] + (round(float(ll[4]), 4), round(float(ll[8]), 4), round(float(ll[9]), 4), ll[10], ll[16])
                elif bim_dict.get(rsid)==(rsid,alt,ref):
                    qtl_dict[rsid]= bim_dict[rsid] + (round(1-float(ll[4]), 4), round(-float(ll[8]), 4), round(float(ll[9]), 4), ll[10], ll[16])
    return qtl_dict

def get_qtl_ukb_dict(qtldir, bim_dict, rsid_dict):
    # match qtl
    qtl_dict = {}
    with gzip.open(qtldir, 'rt') as ss:
        header = next(ss)
        for line in ss:
            ll = line.strip().split()
            rsid = rsid_dict.get(ll[2].replace(':imp:v1','').replace(':','.'))
            alt = ll[3]
            ref = ll[4]
            if rsid in bim_dict:
                if bim_dict.get(rsid)==(rsid,ref,alt):
                    qtl_dict[rsid]=bim_dict[rsid] + (round(float(ll[5]), 4), round(float(ll[9]),4), round(float(ll[10]),4), 10**(-float(ll[12])), ll[7])
                elif bim_dict.get(rsid)==(rsid,alt,ref):
                    qtl_dict[rsid]=bim_dict[rsid] + (round(1-float(ll[5]), 4), round(-float(ll[9]),4), round(float(ll[10]),4), 10**(-float(ll[12])), ll[7])
    return qtl_dict

def match_ss(ssdir, qtldir, rsdir, bimdir, sname, qname, ukb=False):
    rsid_dict = get_rsid_dict(rsdir)
    bim_dict = get_bim_dict(bimdir)
    ss_dict = get_ss_dict(ssdir, rsid_dict, bim_dict)
    if ukb:
        qtl_dict = get_qtl_ukb_dict(qtldir, bim_dict, rsid_dict)
    else:
        qtl_dict = get_qtl_dict(qtldir, bim_dict)
    # get common snps
    common_snps = set(ss_dict.keys()) & set(qtl_dict.keys())
    common_snps_lst = [i for i in bim_dict.keys() if i in common_snps]
    # write out
    pd.DataFrame([ss_dict[i] for i in common_snps_lst], columns=['SNP','A1', 'A2', 'EAF', 'BETA', 'SE', 'P', 'N']).to_csv(sname, sep='\t', header=True, index=False)
    pd.DataFrame([qtl_dict[i] for i in common_snps_lst], columns=['SNP','A1', 'A2', 'EAF', 'BETA', 'SE', 'P', 'N']).to_csv(qname, sep='\t', header=True, index=False)

# PCSK9 - LDL - Alirocumab, Evolucumab - 5231_79 - chr1:55,505,221-55,530,525
# APOB - LDL - Mipomersen - 2797_56 - chr2:21,224,301-21,266,945
# LPL - Tg - Lipoprotein lipase - P06858 - chr8:19,796,764-19,824,770
# ANGPTL3 - Tg - Evinacumab - 10391_1 - chr1:63,063,191-63,071,984
# APOC3 - Tg - Volanesorsen - 6461_54 - chr11:116,700,623-116,703,788

# SNP   A1      A2      EAF     BETA    SE      P       N

# match PCSK9 with LDL(30780)
match_ss('raw/30780_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 'raw/res_invn_X5231_79_Fenland_MA_auto_chrX_filtered_1pc.txt.gz', '/home/wmzh22/utils/ukb/idx/1.rsid', 'bed/PCSK9.bim', 'ss/LDL.PCSK9.LDL.pwcoco', 'ss/LDL.PCSK9.PCSK9.pwcoco')

# match APOB with LDL(30780)
match_ss('raw/30780_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 'raw/res_invn_X2797_56_Fenland_MA_auto_chrX_filtered_1pc.txt.gz', '/home/wmzh22/utils/ukb/idx/2.rsid', 'bed/APOB.bim', 'ss/LDL.APOB.LDL.pwcoco', 'ss/LDL.APOB.APOB.pwcoco')

# match LPL with Tg(30870)
match_ss('raw/30870_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 'raw/discovery_chr8_LPL:P06858:OID20188:v1:Cardiometabolic.gz', '/home/wmzh22/utils/ukb/idx/8.rsid', 'bed/LPL.bim', 'ss/Tg.LPL.Tg.pwcoco', 'ss/Tg.LPL.LPL.pwcoco', ukb=True)

# match ANGPTL3 with Tg(30870)
match_ss('raw/30870_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 'raw/res_invn_X10391_1_Fenland_MA_auto_chrX_filtered_1pc.txt.gz', '/home/wmzh22/utils/ukb/idx/1.rsid', 'bed/ANGPTL3.bim', 'ss/Tg.ANGPTL3.Tg.pwcoco', 'ss/Tg.ANGPTL3.ANGPTL3.pwcoco')

# match APOC3 with Tg(30870)
match_ss('raw/30870_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 'raw/res_invn_X6461_54_Fenland_MA_auto_chrX_filtered_1pc.txt.gz', '/home/wmzh22/utils/ukb/idx/11.rsid', 'bed/APOC3.bim', 'ss/Tg.APOC3.Tg.pwcoco', 'ss/Tg.APOC3.APOC3.pwcoco')
