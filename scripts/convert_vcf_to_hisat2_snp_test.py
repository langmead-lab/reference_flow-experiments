import unittest
from convert_vcf_to_hisat2_snp import convert_vcf_to_hisat2_snp_line

class TestConvertVcfToHisat2SnpLine(unittest.TestCase):
    def test_header(self):
        line = '##bcftools_concatVersion=1.9-206-g4694164+htslib-1.9-258-ga428aa2'
        output = convert_vcf_to_hisat2_snp_line(line)
        assert output == None
    
    def test_snp(self):
        line = 'chr1	51479	.	T	A	.	PASS	AC=531;AN=5096;DP=17461;AF=0.1;EAS_AF=0;EUR_AF=0.19;AFR_AF=0.02;AMR_AF=0.11;SAS_AF=0.23;VT=SNP;NS=2548'
        output = convert_vcf_to_hisat2_snp_line(line)
        assert output == 'chr1_51479\tsingle\tchr1\t51478\tA\n'

    def test_del(self):
        line= 'chr1	83911	.	AAGAG	A	.	PASS	AC=822;AN=5096;DP=17361;AF=0.16;EAS_AF=0.02;EUR_AF=0.14;AFR_AF=0.25;AMR_AF=0.23;SAS_AF=0.17;VT=INDEL;NS=2548'
        output = convert_vcf_to_hisat2_snp_line(line)
        assert output == 'chr1_83911\tdeletion\tchr1\t83910\t4\n'

    def test_ins(self):
        line = 'chr1	84222	.	C	CCT	.	PASS	AC=4995;AN=5096;DP=21311;AF=0.98;EAS_AF=1;EUR_AF=1;AFR_AF=0.93;AMR_AF=0.99;SAS_AF=1;VT=INDEL;NS=2548'
        output = convert_vcf_to_hisat2_snp_line(line)
        assert output == 'chr1_84222\tinsertion\tchr1\t84221\tCT\n'


if __name__ == '__main__':
    unittest.main()

