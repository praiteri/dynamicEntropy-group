from ..orca_tools.extract_from_orca_output import extract_from_orca_output


class template(extract_from_orca_output):
    """
    Extract the input commands
    """

    def read(self):
        return

    def write(self, **kwargs):
        if self.result is None:
            return
        self.logger.info(self.mf.dashes)
        self.logger.info(self.mf.dashes)
        return
