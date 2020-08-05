from kipoi.data import SampleIterator
from mmsplice.junction_dataloader import JunctionPSI5VCFDataloader, \
    JunctionPSI3VCFDataloader
from mmsplice import MMSplice


class SpliceOutlierDataloader(SampleIterator):

    def __init__(self, ref_table, fasta_file, vcf_file,
                 spliceAI_vcf=None, **kwargs):
        self.dl5 = JunctionPSI5VCFDataloader(
            ref_table, fasta_file, vcf_file, **kwargs)
        self.dl3 = JunctionPSI3VCFDataloader(
            ref_table, fasta_file, vcf_file, **kwargs)

    def __next__(self):
        pass


class SpliceOutlier:

    def __init__(self):
        self.mmsplice = MMSplice()

    def predict_on_batch(batch):
        pass

    def predict(self):
        pass

    def predict_save(self):
        pass


class SplicingOutlierResult:

    def __init__(self, df):
        self.df = df
        # junction
        # sample
        # gene
        # splice_site
        # event
        # delta_logit_psi
        # delta_psi
        # ref
        # k or alpha
        # n or beta
        # p_value
        # z-score

    @property
    def junction(self):
        return self.df
        
    @property
    def splice_site(self):
        raise NotImplementedError()


    @property
    def gene(self):
        raise NotImplementedError()
