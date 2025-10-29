from ..orca_tools.extract_from_orca_output import extract_from_orca_output
from ..logger import formatting


class get_charges(extract_from_orca_output):
    def read(self):
        self.result = {}

        self.types = ["mulliken", "loewdin", "hirshfeld", "mbis", "chelpg"]
        self.strings = [
            "MULLIKEN ATOMIC CHARGES",
            "LOEWDIN ATOMIC CHARGES",
            "HIRSHFELD ANALYSIS",
            "MBIS ANALYSIS",
            "CHELPG Charges",
        ]
        offsets = [2, 2, 7, 10, 2]

        for t, s, o in zip(self.types, self.strings, offsets):
            self.result[t] = self.parse_blocks(s, "", offset=o)
        return

    def write(self):
        if self.result is None:
            return

        methods = ["mulliken", "loewdin", "hirshfeld", "mbis", "chelpg"]
        self.logger.info(formatting().dashes_short)
        for i, m in enumerate(methods):
            if self.result[m] is not None:
                self.logger.info(self.strings[i].upper())
                for ll in self.result[m][-1]:
                    if "Sum" not in ll:
                        self.logger.info(ll)
                self.logger.info(formatting().dashes_short)
