from mmsplice.junction_dataloader import JunctionPSI5VCFDataloader, JunctionPSI3VCFDataloader
from mmsplice import MMSplice


class SpliceOutlierPrediction:

    def __init__(self, ref_table, fasta_file, vcf_file):
        self.dl5 = JunctionPSI5VCFDataloader()
        self.dl3 = JunctionPSI3VCFDataloader()
        self.mmsplice = MMSplice()

    def predict_on_batch(batch):
        pass

    def predict(self):
        pass

    def predict_save(self):
        pass
