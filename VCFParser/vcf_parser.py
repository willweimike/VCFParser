import gzip
import pandas as pd
from typing import IO, Dict, Optional, List


class VcfParser:

    header: str
    fh: IO
    VCF_KEYS = [
        'CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    ]

    def __init__(self, vcf: str):
        if vcf.endswith('.gz'):
            self.fh = gzip.open(vcf, 'rt')  # rt: read text
        else:
            self.fh = open(vcf, 'r')
        self.set_header()

    # for context manager
    def __enter__(self):
        return self

    # for context manager
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    # for iterator
    def __iter__(self):
        return self

    # for iterator
    def __next__(self):
        r = self.next()
        if r is not None:
            return r
        else:
            raise StopIteration

    def set_header(self):
        num_header_lines = 0
        for line in self.fh:
            if line.startswith('#'):
                num_header_lines += 1

        self.fh.seek(0)

        temp = []
        i = 0
        for line in self.fh:
            i += 1
            temp.append(line)
            if i == num_header_lines:
                break
        self.header = ''.join(temp)

    def next(self) -> Optional[Dict[str, str]]:
        line = self.fh.readline()

        end_of_file = line == ''
        if end_of_file:
            return None

        temp_lst = line.split('\t')

        count = 0
        assemble_lst = []
        for i in temp_lst:
            count += 1
            assemble_lst.append(i)
            if count == 8:
                break

        data_dict = {}
        for i in range(len(self.VCF_KEYS)):
            data_dict[self.VCF_KEYS[i]] = assemble_lst[i]

        for element in data_dict['INFO'].split(';'):
            if '=' in element:
                pos = element.index('=')
                key = element[0:pos]
                val = element[pos+1:]
            else:
                key, val = element, None
            data_dict[key] = val

        data_dict.pop('INFO')

        return data_dict

    def close(self):
        self.fh.close()


def extract_vcf_to_dataframe(vcf: str, columns: List[str]) -> pd.DataFrame:

    collected_lst = []

    with VcfParser(vcf=vcf) as parser:
        for variant in parser:
            variant = {k: v for k, v in variant.items() if k in columns}
            collected_lst.append(variant)

    return pd.DataFrame(data=collected_lst, columns=columns)
