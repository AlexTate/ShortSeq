from .short_seq import ShortSeq, ShortSeqCounter, read_and_count_fastq
from .short_seq_var import ShortSeqVar, get_domain_var
from .short_seq_128 import ShortSeq128, get_domain_128
from .short_seq_64 import ShortSeq64, get_domain_64

MIN_VAR_NT, MAX_VAR_NT = get_domain_var()
MIN_128_NT, MAX_128_NT = get_domain_128()
MIN_64_NT,  MAX_64_NT  = get_domain_64()
