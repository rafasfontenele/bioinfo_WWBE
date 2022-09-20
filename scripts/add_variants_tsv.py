#!/user/bin/env python
import pandas as pd
import sys
import vcf

GENE_MAP = {'ORF1a': [266, 13468],
            'ORF1b': [13468, 21555],
            'S': [21563, 25384],
            'ORF3a': [25393, 26220],
            'E': [26245, 26472],
            'M': [26523, 27191],
            'ORF6': [27202, 27387],
            'ORF7a': [27394, 27759],
            'ORF7b': [27756, 27887],
            'ORF8': [27894, 28259],
            'N': [28274, 29533],
            'ORF10': [29558, 29674]}


def read_lofreq(filename):
    lofreq_list  = []
    #lofreq_calls = pd.DataFrame(columns=["REGION", "POS", "REF", "ALT", "QUAL", "FILTER", "REF_DP", "REF_RV", "ALT_DP", "ALT_RV", "AF", "TOTAL_DP", "STRAND-BIAS", "EFFECT"])
    vcf_reader = vcf.Reader(filename=filename)
    for row in vcf_reader:
        if len(row.FILTER) == 0 and "EFF" not in row.INFO:
            lofreq_list.append({"REGION": row.CHROM,
                                                "POS": int(row.POS),
                                                "REF": str(row.REF),
                                                "ALT": str(row.ALT[0]),
                                                "QUAL": row.QUAL,
                                                "FILTER": "",
                                                "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1],
                                                "REF_RV": row.INFO["DP4"][1],
                                                "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                                "ALT_RV": row.INFO["DP4"][3],
                                                "AF": row.INFO["AF"],
                                                "TOTAL_DP": row.INFO["DP"],
                                                "STRAND-BIAS": row.INFO["SB"],
                                                "EFF": ""
                                                })
        if len(row.FILTER) > 0 and "EFF" in row.INFO:
            lofreq_list.append({"REGION": row.CHROM,
                                                "POS": int(row.POS),
                                                "REF": str(row.REF),
                                                "ALT": str(row.ALT[0]),
                                                "QUAL": row.QUAL,
                                                "FILTER": row.FILTER[0],
                                                "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1],
                                                "REF_RV": row.INFO["DP4"][1],
                                                "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                                "ALT_RV": row.INFO["DP4"][3],
                                                "AF": row.INFO["AF"],
                                                "TOTAL_DP": row.INFO["DP"],
                                                "STRAND-BIAS": row.INFO["SB"],
                                                "EFF": row.INFO["EFF"]
                                                })
        if len(row.FILTER) == 0 and "EFF" in row.INFO:
            lofreq_list.append({"REGION": row.CHROM,
                                                "POS": int(row.POS),
                                                "REF": str(row.REF),
                                                "ALT": str(row.ALT[0]),
                                                "QUAL": row.QUAL,
                                                "FILTER": "",
                                                "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1],
                                                "REF_RV": row.INFO["DP4"][1],
                                                "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                                "ALT_RV": row.INFO["DP4"][3],
                                                "AF": row.INFO["AF"],
                                                "TOTAL_DP": row.INFO["DP"],
                                                "STRAND-BIAS": row.INFO["SB"],
                                                "EFF": row.INFO["EFF"],
                                                })
        if len(row.FILTER) > 0 and "EFF" not in row.INFO:
            lofreq_list.append({"REGION": row.CHROM,
                                                "POS": int(row.POS),
                                                "REF": str(row.REF),
                                                "ALT": str(row.ALT[0]),
                                                "QUAL": row.QUAL,
                                                "FILTER": row.FILTER[0],
                                                "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1],
                                                "REF_RV": row.INFO["DP4"][1],
                                                "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                                "ALT_RV": row.INFO["DP4"][3],
                                                "AF": row.INFO["AF"],
                                                "TOTAL_DP": row.INFO["DP"],
                                                "STRAND-BIAS": row.INFO["SB"],
                                                "EFF": ""
                                                })
        lofreq_calls = pd.DataFrame.from_records(lofreq_list)

    if lofreq_calls.empty:
        return lofreq_calls
    lofreq_calls["Variant"] = lofreq_calls.apply(lambda row:
                                                 str(row["POS"]) + ":" +
                                                 str(row["ALT"][1:]) if len(row["ALT"]) > 1 else
                                                 (str(row["POS"]) + "-" +
                                                  str(row["POS"] + len(row["REF"][1:]) + 1) if len(row["REF"]) > 1 else
                                                  str(row["REF"]) + str(row["POS"]) +
                                                  str(row["ALT"])), axis=1)
    lofreq_calls["EFFECT"] = lofreq_calls.apply(lambda row:
                                                str(row["EFF"]).split("|")[1] if row["EFF"] else
                                                "SILENT", axis=1)
    lofreq_calls["AA_MUT"] = lofreq_calls.apply(lambda row:
                                                str(row["EFF"]).split("|")[3] if row["EFF"] else
                                                "", axis=1)
    lofreq_calls["GENE"] = lofreq_calls.apply(lambda row:
                                              str(row["EFF"]).split("|")[5] if row["EFF"] else
                                              "Non-coding", axis=1)
    lofreq_calls = lofreq_calls[lofreq_calls["FILTER"] != "AmpliconRemoval"]
    return lofreq_calls


def main():
    file = sys.argv[1]
    new_file = read_lofreq(file)
    new_file.to_csv(sys.argv[2], sep="\t", index=False)


if __name__ == "__main__":
    main()
