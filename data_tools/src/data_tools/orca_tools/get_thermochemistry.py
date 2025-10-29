from ..orca_tools.extract_from_orca_output import extract_from_orca_output
from ..logger import formatting


class get_thermochemistry(extract_from_orca_output):
    """
    Extract the thermochemistry
    """

    def read(self):
        """
        Extract the thermochemistry
        """
        text_block = self.parse_blocks("THERMOCHEMISTRY", ["G-E(el)"], offset=0)

        if text_block is None:
            return

        self.result = []
        self.dict_of_fields = {
            "Temperature": "i",
            "Pressure": "i",
            # Summary of contributions to the inner energy U:
            "Electronic energy": "i",
            "Zero point energy": "i",
            "Thermal vibrational correction": "v",
            "Thermal rotational correction": "v",
            "Thermal translational correction": "v",
            "Total thermal energy": "i",
            # Summary of corrections to the electronic energy
            "Total thermal correction": "v",
            "Non-thermal (ZPE) correction": "v",
            # Enthalpy
            "Thermal Enthalpy correction": "v",
            "Total Enthalpy": "i",
            # The entropy contributions are T*S
            "Electronic entropy": "v",
            "Vibrational entropy": "v",
            "Rotational entropy": "v",
            "Translational entropy": "v",
            "Final entropy term": "i",
            # The Gibbs free energy is G = H - T*S
            "Final Gibbs free energy": "i",
            # "G-E(el)" : "v",
        }

        for block in text_block:
            thermochemistry = {}
            for ll in self.dict_of_fields.keys():
                x = self.parse_single_line(self.search_string(ll, block), block)[-1]
                if "Eh" in x:
                    x = x.split()
                    i = x.index("Eh")
                    thermochemistry[ll] = x[i - 1] + " " + x[i]
                else:
                    thermochemistry[ll] = x.split("...")[-1]
            self.result.append(thermochemistry)

    def write(self, **kwargs):
        if self.result is None:
            return

        self.logger.info(formatting().dashes)
        self.logger.info("Thermochemistry summary".upper())
        for thermochemistry in self.result:
            for k, v in self.dict_of_fields.items():
                if v == "i":
                    self.logger.info(self.fmt.format(k, thermochemistry[k]))
                else:
                    self.logger.verbose(self.fmt.format(k, thermochemistry[k]))
            if len(self.result) > 1:
                self.logger.info(formatting().dashes)
        # self.logger.info(formatting().dashes)
        return
