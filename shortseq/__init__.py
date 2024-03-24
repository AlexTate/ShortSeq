from .short_seq import (
    pack,
    from_str,
    from_bytes
)

from .short_seq_var import ShortSeqVar, get_domain_var
from .short_seq_192 import ShortSeq192, get_domain_192
from .short_seq_64 import ShortSeq64, get_domain_64
from .counter import ShortSeqCounter, read_and_count_fastq

MIN_VAR_NT, MAX_VAR_NT = get_domain_var()
MIN_192_NT, MAX_192_NT = get_domain_192()
MIN_64_NT,  MAX_64_NT  = get_domain_64()
