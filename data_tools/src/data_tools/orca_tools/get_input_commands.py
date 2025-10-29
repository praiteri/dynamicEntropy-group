from ..orca_tools.extract_from_orca_output import extract_from_orca_output


class get_input_commands(extract_from_orca_output):
    """
    Extract the input commands
    """

    def read(self):
        string = "INPUT FILE"
        escape_string = "END OF INPUT"
        self.read_blocks(string, escape_string, offset=3)

    def write(self, **kwargs):
        if self.result is None:
            return

        self.logger.info(self.mf.dashes)
        for ll in self.result[0][:-1]:
            self.logger.info(ll.split(">")[-1])
        self.logger.info(self.mf.dashes)
        return
