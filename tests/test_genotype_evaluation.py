import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import pysam

PATH = Path(__file__).resolve().parents[1]/'scripts/evaluate_genotypes.py'
SPEC = importlib.util.spec_from_file_location('genotype_evaluator', PATH)
EVALUATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EVALUATOR)


class GenotypeEvaluationTests(unittest.TestCase):
    def test_partial_alt_is_not_a_correct_heterozygote_even_if_gtmatch_is_zero(self):
        header=('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
                '##INFO=<ID=MatchId,Number=.,Type=String,Description="Match">\n'
                '##INFO=<ID=GTMatch,Number=1,Type=Integer,Description="AC difference">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
        def row(pos,key,gt):
            return f'chr1\t{pos}\t{key}\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100;MatchId={key},{key};GTMatch=0\tGT\t{gt}\n'
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for label,gts in [('tp-base',['0|1','0|1']),('tp-comp',['1|0','1/.']),('fp',[]),('fn',[])]:
                plain=root/(label+'.vcf');plain.write_text(header+''.join(row(500+i*500,str(i),gt) for i,gt in enumerate(gts)))
                pysam.tabix_compress(str(plain),str(plain)+'.gz',force=True)
                pysam.tabix_index(str(plain)+'.gz',preset='vcf',force=True)
            (root/'summary.json').write_text(json.dumps({'TP-base':2,'TP-comp':2,'FP':0,'FN':0,
                'precision':1,'recall':1,'f1':1,'gt_concordance':1}))
            result,_,_=EVALUATOR.evaluate(root)
        self.assertEqual(result['base']['concordant'],1)
        self.assertEqual(result['base']['partial_or_missing'],1)
        self.assertEqual(result['base']['concordance_fully_called'],1)
        self.assertEqual(result['base']['concordance_all_matched'],.5)
        self.assertEqual(result['genotype_aware']['f1'],.5)


if __name__=='__main__':unittest.main()
