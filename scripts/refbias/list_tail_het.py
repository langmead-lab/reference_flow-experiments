import argparse
import pandas as pd

def main(fn_txt, fn_perc, fn_val, fn_out):
    f_out = open(fn_out, 'w')
    string = "File: " + fn_txt
    f_out.write(string)
    f_out.write("\n")

    df = pd.read_csv(fn_txt, sep='\t')
    if fn_val:
        string = "Reference Bias Threshold Value: " + str(fn_val)
        f_out.write(string)
        f_out.write("\n")
        df_pass = df[df['REFERENCE BIAS'] >= float(fn_val)]
        for element in df_pass['HET SITE']:
            f_out.write(str(element))
            f_out.write("\n")
    else:
        string = "Reference Bias Percentile: " + str(fn_perc)
        f_out.write(string)
        f_out.write("\n")
        rb = df['REFERENCE BIAS']
        thresh = (rb.quantile(q=[float(fn_perc)]))
        df_pass = df[df['REFERENCE BIAS'] >= float(thresh)]
        for element in df_pass['HET SITE']:
            f_out.write(str(element))
            f_out.write("\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--txt', help='reference bias text file')
    parser.add_argument('-o', '--out', help='output file')
    parser.add_argument('-p', '--perc', help='percentile specified by user')
    parser.add_argument('-v', '--val', help='value specified by user')

    args = parser.parse_args()
    fn_txt = args.txt
    fn_out = args.out

    fn_perc = args.perc
    fn_val = args.val


    print("fn_txt: ", fn_txt)
    print("fn_perc: ", fn_perc)
    print("fn_val: ", fn_val)
    print("fn_out: ", fn_out)

    main(fn_txt, fn_perc, fn_val, fn_out)
