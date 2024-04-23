class TunableParameters:

    def __init__(self, usage, target_gc, repeated_nucleotides, redundant_encoding):

        self._usage = usage
        self._target_gc = target_gc
        self._repeated_nucleotides = repeated_nucleotides
        self._redundant_encoding = redundant_encoding

    def get_usage(self):
        return self._usage

    def get_target_gc(self):
        return self._target_gc

    def get_repeated_nucleotides(self):
        return self._repeated_nucleotides

    def get_redundant_encoding(self):
        return self._redundant_encoding
