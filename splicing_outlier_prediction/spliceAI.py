"""
This file provides python interface for spliceAI
"""
from collections import namedtuple
from tqdm import tqdm
import numpy as np
import pandas as pd
from kipoiseq import Variant
from spliceai.utils import Annotator, get_delta_scores
from mmsplice.utils import df_batch_writer
from kipoiseq.extractors import MultiSampleVCF


class VariantDB:

    def __init__(self, path):
        import rocksdb
        self.db = rocksdb.DB(path, rocksdb.Options(
            create_if_missing=True), read_only=True)

    @staticmethod
    def _variant_to_byte(variant):
        return bytes(str(variant), 'utf-8')

    def _type(self, value):
        raise NotImplementedError()

    def _get(self, variant):
        if variant.startswith('chr'):
            variant = variant[3:]
        return self.db.get(self._variant_to_byte(variant))

    def __getitem__(self, variant):
        value = self._get(variant)
        if value:
            return self._type(value)
        else:
            raise KeyError('This variant "%s" is not in the db'
                           % str(variant))

    def __contains__(self, variant):
        return self._get(variant) is not None

    def get(self, variant, default=None):
        try:
            return self[variant]
        except KeyError:
            return default

    def items(self):
        it = self.db.iteritems()
        it.seek_to_first()
        for variant, value in it:
            yield variant.decode('utf-8'), self._type(value)


class SpliceAIDB(VariantDB):

    def _type(self, value):
        return list(self._parse(value))

    def _parse(self, value):
        for i in value.decode('utf-8').split(';'):
            results = i.split('|')
            scores = np.array(list(map(float, results[1:])))
            yield SpliceAI.Score(
                results[0], scores[:4].max(), *scores
            )


class SpliceAI:
    Score = namedtuple('Score', ['gene_name', 'delta_score',
                                 'acceptor_gain', 'acceptor_loss',
                                 'donor_gain', 'donor_loss',
                                 'acceptor_gain_position',
                                 'acceptor_loss_positiin',
                                 'donor_gain_position',
                                 'donor_loss_position'])
    Record = namedtuple('Record', ['chrom', 'pos', 'ref', 'alts'])

    def __init__(self, fasta=None, annotation=None, db_path=None,
                 dist=50, mask=1, samples=False, quality=False):
        """
        Args:
          fasta: fasta file path
          annotation: 'grch37' or 'grch38'
          dist: area of interest based on distance to variant
          mask: mask for 'N'
        """
        assert ((fasta is not None) and (annotation is not None)) \
            or (db_path is not None)
        self.db_only = fasta is None
        if not self.db_only:
            self.ann = Annotator(fasta, annotation)
        self.dist = dist
        self.mask = mask
        self.db = SpliceAIDB(db_path) if db_path else None
        self.samples = samples
        self.quality = quality

    @staticmethod
    def _to_record(variant):
        if type(variant) == str:
            variant = Variant.from_str(variant)
        return SpliceAI.Record(variant.chrom, variant.pos,
                               variant.ref, [variant.alt])

    @staticmethod
    def parse(output):
        results = output.split('|')
        results = [0 if i == '.' else i for i in results]
        scores = np.array(list(map(float, results[2:])))
        return SpliceAI.Score(
            results[1], scores[:4].max(), *scores
        )

    def predict(self, variant):
        record = self._to_record(variant)
        if self.db:
            try:
                return self.db[str(variant)]
            except KeyError:
                if self.db_only:
                    return []
        return [
            self.parse(i)
            for i in get_delta_scores(record, self.ann,
                                      self.dist, self.mask)
        ]

    def predict_df(self, variants, vcf=None):
        rows = self._predict_df(variants, vcf)
        columns = [
            'variant', 'gene_name', 'delta_score',
            'acceptor_gain', 'acceptor_loss',
            'donor_gain', 'donor_loss',
            'acceptor_gain_position',
            'acceptor_loss_positiin',
            'donor_gain_position',
            'donor_loss_position'
        ]
        if self.samples:
            columns.append('samples')
            if self.quality:
                columns.append('GQ')
                columns.append('DP_ALT')

        return pd.DataFrame(rows, columns=columns).set_index('variant')

    def _predict_df(self, variants, vcf=None):
        for v in variants:
            for score in self.predict(v):
                row = score._asdict()
                row['variant'] = str(v)

                if self.samples:
                    samples = vcf.get_samples(v)
                    row['samples'] = ';'.join(samples)

                    if self.quality:
                        GQ = (
                            v.source.gt_quals[vcf.sample_mapping[sample]]
                            for sample in samples
                        )
                        row['GQ'] = ';'.join(map(str, GQ))
                        dp_alt = (
                            v.source.gt_alt_depths[vcf.sample_mapping[sample]]
                            for sample in samples
                        )
                        row['DP_ALT'] = ';'.join(map(str, dp_alt))

                yield row

    def _predict_on_vcf(self, vcf_file, batch_size=100):
        vcf = MultiSampleVCF(vcf_file)
        for variants in vcf.batch_iter(batch_size):
            yield self.predict_df(variants, vcf).reset_index('variant')

    def predict_save(self, vcf_file, output_csv,
                     batch_size=100, progress=True):
        batches = self._predict_on_vcf(vcf_file, batch_size=batch_size)
        if progress:
            batches = iter(tqdm(batches))
        df_batch_writer(batches, output_csv)
